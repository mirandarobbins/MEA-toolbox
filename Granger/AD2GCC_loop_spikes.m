spikes = meta_spike_extract(data,params,4);
spikes = meta_MUA64plot(data,params,spikes);

tic
% granger.N=size(data.trimmed,1)%*size(data.trimmed,3);
granger=[];
granger.PVAL    =   0.001;       % probability threshold for Granger causality significance
granger.NLAGS   =   2;         % if -1, best model order is assessed automatically
granger.freqs   =   1:100;      % frequency range to analyze (spectral analysis only)
for idx=1%:size(spikes.PDF.PDF_trimmed{1},1)
    %%%%%For raw MEA data
    params.selected_rep=idx;
    % granger.Fs      =   2000;         % params.Fs;   % sampling frequency (for spectral analysis only)
    % granger.X       =   downsample(data.filtered_lfp(4000:20000,:,params.selected_rep),10)'; % One rep
    % granger.X       =   downsample(reshape(data.filtered_spikes(4000:20000,:,:),160010,64,1),10)'; % All reps
    
    % %%%%%For convolved NeuroExplore spike times
    % granger.Fs      =   1000;         % params.Fs;   % sampling frequency  (for spectral analysis only)
    % for chan_id=1:numel(spikes.PDF_trimmed)
    % %     granger.X(chan_id,:)= reshape(spikes.PDF_trimmed{chan_id},1,numel(spikes.PDF_trimmed{chan_id})); % all trials
    %     granger.X(chan_id,:)= spikes.PDF_trimmed{chan_id}(1,:); % selected rep
    % end; granger.X(isnan(granger.X))=0;
    
    %%%%%For convolved MUA spike times
    granger.Fs      =   1000;         % params.Fs;   % sampling frequency  (for spectral analysis only)
    for chan_id=1:numel(spikes.PDF.PDF_trimmed)
            granger.X(chan_id,:)= spikes.PDF.PDF_trimmed{chan_id}(idx,200:1000); % selected rep
%         granger.X(chan_id,:)= nansum(spikes.PDF.PDF_trimmed{chan_id},1); % summed spikes
        
    end; granger.X(isnan(granger.X))=0;
    
    granger.nvar    =   size(granger.X,1);       % no input channels
    granger.N       =   size(granger.X,2);       % number of data points
    granger.sfile   =[data.this_file,'_sweep_',num2str(params.selected_rep),'.net'];
    % detrend and demean data
    % disp('detrending and demeaning data...');
    granger.X = cca_detrend(granger.X);
    granger.X = cca_rm_temporalmean(granger.X);
    
    % check covariance stationarity
    disp('checking for covariance stationarity ...');
    granger.uroot = cca_check_cov_stat(granger.X,10);
    granger.inx = find(granger.uroot);
    if sum(granger.uroot) == 0,
        disp('OK, data is covariance stationary by ADF');
    else
        disp('WARNING, data is NOT covariance stationary by ADF');
        disp(strcat('unit roots found in channels: ',num2str(granger.inx')));
    end
    
    % check covariance stationarity again using KPSS test
    [granger.kh,granger.kpss] = cca_kpss(granger.X);
    granger.inx = find(granger.kh==0);
    if isempty(granger.inx),
        disp('OK, data is covariance stationary by KPSS');
    else
        disp('WARNING, data is NOT covariance stationary by KPSS');
        disp(strcat('unit roots found in channels: ',num2str(granger.inx)));
    end
    
    
    % find best model order
    if granger.NLAGS == -1,
        disp('finding best model order ...');
        [granger.bic,granger.aic] = cca_find_model_order(granger.X,2,12);
        disp(['best model order by Bayesian Information Criterion = ',num2str(granger.bic)]);
        disp(['best model order by Aikaike Information Criterion = ',num2str(granger.aic)]);
        granger.NLAGS = granger.aic;
    end
    
    
    
    %% -------------------------------------------------------------------------
    % analyze time-domain granger
    
    % find time-domain conditional Granger causalities [THIS IS THE KEY FUNCTION]
    disp('finding conditional Granger causalities ...');
    granger.ret = cca_granger_regress(granger.X,granger.NLAGS,1);   % STATFLAG = 1 i.e. compute stats
    
    % check that residuals are white
    granger.dwthresh = 0.05/granger.nvar;    % critical threshold, Bonferroni corrected
    granger.waut = zeros(1,granger.nvar);
    for ii=1:granger.nvar,
        if granger.ret.waut<granger.dwthresh,
            granger.waut(ii)=1;
        end
    end; clear ii
    granger.inx = find(granger.waut==1);
    
    if isempty(granger.inx),
        disp('All residuals are white by corrected Durbin-Watson test');
    else
        disp(['WARNING, autocorrelated residuals in variables: ',num2str(granger.inx)]);
    end
    
    
    % check model consistency, ie. proportion of correlation structure of the
    % data accounted for by the MVAR model
    if granger.ret.cons>=80,
        disp(['Model consistency is OK (>80%), value=',num2str(granger.ret.cons)]);
    else
        disp(['Model consistency is <80%, value=',num2str(granger.ret.cons)]);
    end
    
    
    % analyze adjusted r-square to check that model accounts for the data (2nd
    % check)
    granger.rss = granger.ret.rss_adj;
    granger.inx = find(granger.rss<0.3);
    if isempty(granger.inx)
        disp(['Adjusted r-square is OK: >0.3 of variance is accounted for by model, val=',num2str(mean(granger.rss))]);
    else
        disp(['WARNING, low (<0.3) adjusted r-square values for variables: ',num2str(granger.inx)]);
        disp(['corresponding values are ',num2str(rss(granger.inx))]);
        disp('try a different model order');
    end
    %%
    granger.PVAL    =   0.001;       % probability threshold for Granger causality significance
    % find significant Granger causality interactions (Bonferonni correction)
    [granger.PR,granger.q] = cca_findsignificance(granger.ret,granger.PVAL,1);
    disp(['testing significance at P < ',num2str(granger.PVAL), ', corrected P-val = ',num2str(granger.q)]);
    
    % extract the significant causal interactions only
    granger.GC = granger.ret.gc;
    granger.GC2 = granger.GC.*granger.PR;
    
    
    % calculate causal connectivity statistics
    disp('calculating causal connectivity statistics');
    granger.causd = cca_causaldensity(granger.GC,granger.PR);
    granger.causf = cca_causalflow(granger.GC,granger.PR);
    
    
    disp(['time-domain causal density = ',num2str(granger.causd.cd)]);
    disp(['time-domain causal density (weighted) = ',num2str(granger.causd.cdw)]);
    
    % create Pajek readable file
    cca_pajek(granger.PR,granger.GC,granger.sfile);
    
    
    % -------------------------------------------------------------------------
    % plot time-domain granger results
%     figure(1); clf reset;
%     granger.FSIZE = 8;
%     colormap(flipud(bone));
%     % plot raw time series
%     %   for i=2:granger.nvar,
%     %       granger.X(i,:) = granger.X(i,:)+(10*(i-1));
%     %   end; clear i
%     %   subplot(231);
%     set(gca,'FontSize',granger.FSIZE);
%     plot(granger.X');
%     axis('square');
%     set(gca,'Box','off');
%     xlabel('time');
%     set(gca,'YTick',[]);
%     xlim([0 granger.N]);
%     title('Causal Connectivity Toolbox v2.0');
    
    % plot granger causalities as matrix
%     figure(2); clf reset;
%     % subplot(232);
%     set(gca,'FontSize',granger.FSIZE);
%     imagesc(granger.GC2);
%     axis('square');
%     set(gca,'Box','off');
%     title(['Granger causality, p<',num2str(granger.PVAL)]);
%     xlabel('from');
%     ylabel('to');
%     set(gca,'XTick',[1:granger.N]);
%     set(gca,'XTickLabel',1:granger.N);
%     set(gca,'YTick',[1:granger.N]);
%     set(gca,'YTickLabel',1:granger.N);
    
    % plot granger causalities as a network
    granger.GC2(granger.GC2==0)=NaN;
    % figure(3); clf reset;
    % figure; cca_plotcausality(granger.GC2,[],1); camroll(270); set(gca,'pos',[0.1 0.1 0.8 0.8]);
    
    figure; cca_plotcausality64(granger.GC2,[],1); camroll(270); set(gca,'pos',[0.1 0.1 0.8 0.8]);
    title(['Granger Causality calculated over ', num2str(granger.NLAGS),' lags, showing links with P<',num2str(granger.PVAL),', trial no: ',num2str(params.selected_rep)])
end
toc
