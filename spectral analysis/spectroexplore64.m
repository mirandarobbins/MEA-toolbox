function [data params spectro fourier] = spectroexplore64(data,params)
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=1;
end
spectro.analysis_window=[params.search_win(1),40;params.search_win(2),60]; %i.e.... min time , min freq ; max time ; max freq
%% Calculate spectrogram on a swep-by-sweep basis for strongest channel
spectro.params.freqs=0:1:100;

t1=tic;
progbar = waitbar(0,'Initializing...',...
            'name','sweep analysis progress');%,...
%             'position',[640 600 275 50]);
for sweep_id=data.sweep_sort.successful_sweeps
    spectro.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
end
       
    for sweep_id=data.sweep_sort.successful_sweeps
        
    % update progress bar
    waitbar(sweep_id/params.last_sweep,progbar,...
    strcat(['Analysing sweep ' num2str(sweep_id)  '/' num2str(params.last_sweep) '...'] ))

    % Run the spectrogram
    [spectro.S(:,:,sweep_id),...
     spectro.F(:,:,sweep_id),...
     spectro.T(:,:,sweep_id),...
     spectro.P(:,:,sweep_id)]   =...
            spectrogram(data.filtered_lfp(:,spectro.StrongestChannel(2),sweep_id),...
            hamming(1000),...
            999,...
            spectro.params.freqs,...
            params.Fs,...
            'yaxis');
    end

    % finalise progress bar
waitbar(1,progbar,  ['Finished! Processed ',num2str(params.last_sweep),' Sweeps in ', num2str(toc(t1)), ' seconds.']);
clear t1
%%
    spectro.T=max(spectro.T,[],3);
    spectro.F=max(spectro.F,[],2);
    spectro.P    = 10*log10(abs(spectro.P)); %convert to dB
    temp=spectro.P(:,:,data.sweep_sort.successful_sweeps);
    spectro.mean = flipud(mean(temp,3)); 
    spectro.SD   = flipud(std(temp,0,3)); clear temp
%     temp=mean(spectro.mean(:,1),2);
    temp=mean(spectro.mean(:,params.baseline_win_samples(1):params.baseline_win_samples(2)),2);

    spectro.mean_norm=spectro.mean-repmat(temp,1,size(spectro.mean,2)); clear temp


spectro.mean_smooth=smooth2a(spectro.mean,20, 200); %2D smoothing
spectro.mean_norm_smooth=smooth2a(spectro.mean_norm,20, 200); %2D smoothing
spectro.SD_smooth=smooth2a(spectro.SD,20, 200); %2D smoothing
spectro.tb=(0:max(data.tb)-1)./params.Fs;
spectro.mean_correct=zeros(size(spectro.mean));

for idx=1:size(spectro.mean,2)
    spectro.mean_correct(:,idx)=spectro.mean(:,idx)./mean(spectro.mean(:,1),2);
end; clear idx
%% Plot
% x_temp=repmat(spectro.T,params.last_sweep,1);x_temp=x_temp';

if plotYN==1
h_specto=figure;
subplot(211); hold on
    plot(data.tb,sgolayfilt(squeeze(data.filtered_lfp(:,spectro.StrongestChannel(2),data.sweep_sort.successful_sweeps)),0,15),'Color',[0.8 0.8 0.8]); 
    plot(data.tb,mean(squeeze(data.filtered_lfp(:,spectro.StrongestChannel(2),data.sweep_sort.successful_sweeps)),2),'-k')
    axis([min(min(spectro.T))*1000 max(max(spectro.T))*1000 -0.2 0.2])
    xlabel('time (ms)')    
    ylabel('fEPSP potential (mV)')
    box off
ax1=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);

subplot(212); hold on
%     imagesc(repmat(spectro.T,numel(spectro.params.freqs),1)',...
%             repmat(spectro.params.freqs,numel(spectro.T),1),...
%             spectro.mean_norm');

imagesc(flipud(spectro.mean_norm_smooth));
        axis tight; axis xy
    colormap(hot); %caxis([-80 -40]) ;    
ax2a=gca;
pos=get(ax2a,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);
pos=get(gca,'pos');
    hc=colorbar('location','eastoutside',...
                'position',[pos(1)+pos(3)+0.01 pos(2) 0.03 pos(4)*0.5],...
                'YAxisLocation','right');
    ylabel(hc,'dBuV/Hz')
title('Power Spectral Density Estimate')
ylabel('Frequency (Hz)');
xlabel('Time (s)');  
ax2a = gca;
ax2b = axes('Position',get(ax2a,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k',...
           'Xtick',[],...
           'YTick',[],...
           'xlim',[min(min(spectro.T)), max(max(spectro.T))*1000],...
           'ylim',[min(min(spectro.F)), max(max(spectro.F))]);
rectangle('Position',[spectro.analysis_window(1,1),...
                      spectro.analysis_window(1,2),...
                      spectro.analysis_window(2,1)-spectro.analysis_window(1,1),...
                      spectro.analysis_window(2,2)-spectro.analysis_window(1,2)],...
                      'LineWidth',2,'LineStyle','-','EdgeColor',[0 0 0],'Parent',ax2b)

clear ax1 ax2a ax2b hc pos x_temp
end 
%% trial to trial aveage power and variability
for sweep_id=data.sweep_sort.successful_sweeps
spectro.average_window_power(sweep_id)=mean(mean(spectro.P(spectro.analysis_window(1,2):spectro.analysis_window(2,2),...
                                                 spectro.analysis_window(1):spectro.analysis_window(2),...
    sweep_id)));
end
clear sweep_id
spectro.average_window_power_mean=nanmean(spectro.average_window_power);
spectro.average_window_power_sem=nansem(spectro.average_window_power);

%% calculate FFT using windowed data
StrongestChannel=spectro.StrongestChannel(2);
fourier=[];
fourier.input    = squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),StrongestChannel,:));%data.sweep_sort.successful_sweeps);
fourier.baseline = squeeze(data.filtered_lfp(params.baseline_win_samples(1):params.baseline_win_samples(2),StrongestChannel,:));%data.sweep_sort.successful_sweeps);

fourier.T = 1/params.Fs;                               % Sample time
fourier.L = size(fourier.input,1);                     % Length of signal

fourier.T_baseline = 1/params.Fs;                      % Sample time
fourier.L_baseline = size(fourier.baseline,1);         % Length of signal

% fourier.NFFT =2^nextpow2(fourier.L); % Next power of 2 from length of y
fourier.NFFT =2^15; % Next power of 2 from length of y


fourier.f = params.Fs/2*linspace(0,1,fourier.NFFT/2+1); %frequency vector for plotting

fourier.degree=0;
fourier.frame=15;


t1=tic;
progbar = waitbar(0,'Initializing...',...
            'name','FFT progress')%,...
%             'position',[640 600 275 50]);

for sweep_id=1:size(fourier.input,2);%data.sweep_sort.successful_sweeps
    
    % update progress bar
    waitbar(sweep_id/params.last_sweep,progbar,...
    strcat(['Analysing sweep ' num2str(sweep_id)  '/' num2str(params.last_sweep) '...'] ))


fourier.fft(:,sweep_id) = fft(fourier.input(:,sweep_id),fourier.NFFT)/(fourier.L/2);
fourier.psd(:,sweep_id) = 10*log10(abs(fourier.fft(:,sweep_id))/(fourier.L/2));
fourier.fft_baseline(:,sweep_id) = fft(fourier.baseline(:,sweep_id),fourier.NFFT)/(fourier.L_baseline/2);
fourier.psd_baseline(:,sweep_id) = 10*log10(abs(fourier.fft_baseline(:,sweep_id))/(fourier.L_baseline/2));

fourier.fft_norm(:,sweep_id)=abs(fourier.fft(:,sweep_id))./(abs(fourier.fft_baseline(:,sweep_id))); % subtractive normalisation here?
% fourier.psd(:,sweep_id)=fourier.psd(:,sweep_id)./fourier.psd(1,sweep_id)

end
    % finalise progress bar
waitbar(sweep_id/params.last_sweep,progbar,  ['Finished: Processed ',num2str(params.last_sweep),' Sweeps in ', num2str(toc(t1)), ' seconds.']);
clear t1
%     spectro.mean=mean(10*log10(abs(spectro.P)),3);
%      spectro.SD=std(10*log10(abs(spectro.P)),0,3);

fourier.fft_mean_success = nanmean(fourier.fft(:,data.sweep_sort.successful_sweeps),2);
fourier.fft_SEM_success  = nansem(fourier.fft(:,data.sweep_sort.successful_sweeps),2);

fourier.fft_mean_failure = nanmean(fourier.fft(:,data.sweep_sort.failures),2);
fourier.fft_SEM_failure  = nansem(fourier.fft(:,data.sweep_sort.failures),2);

fourier.fft_norm_mean_success = nanmean(fourier.fft_norm(:,data.sweep_sort.successful_sweeps),2);
fourier.fft_norm_SEM_success  = nansem(fourier.fft_norm(:,data.sweep_sort.successful_sweeps),2);

fourier.fft_norm_mean_failure = nanmean(fourier.fft_norm(:,data.sweep_sort.failures),2);
fourier.fft_norm_SEM_failure  = nansem(fourier.fft_norm(:,data.sweep_sort.failures),2);

fourier.psd_mean_success = nanmean(fourier.psd(:,data.sweep_sort.successful_sweeps),2);
fourier.psd_SEM_success  = nansem(fourier.psd(:,data.sweep_sort.successful_sweeps),2);

fourier.psd_mean_failure = nanmean(fourier.psd(:,data.sweep_sort.failures),2);
fourier.psd_SEM_failure  = nansem(fourier.psd(:,data.sweep_sort.failures),2);
clear progbar sweep_id
%% ploting FFT
if plotYN==1
   figure; hold on;
for sweep_id=data.sweep_sort.successful_sweeps       
    plot(fourier.f,fourier.psd(1:fourier.NFFT/2+1,sweep_id),':b')
end
semilogx(fourier.f,fourier.psd_mean_success(1:fourier.NFFT/2+1),'b','LineWidth',1.5)
ciplot(fourier.psd_mean_success(1:fourier.NFFT/2+1)-fourier.psd_SEM_success(1:fourier.NFFT/2+1),...
       fourier.psd_mean_success(1:fourier.NFFT/2+1)+fourier.psd_SEM_success(1:fourier.NFFT/2+1),...
       fourier.f,'b')
if ~isempty(data.sweep_sort.failures)
for sweep_id=data.sweep_sort.failures       
    plot(fourier.f,fourier.psd(1:fourier.NFFT/2+1,sweep_id),':r')
end   
semilogx(fourier.f,fourier.psd_mean_failure(1:fourier.NFFT/2+1),'r','LineWidth',1.5)
ciplot(fourier.psd_mean_failure(1:fourier.NFFT/2+1)-fourier.psd_SEM_failure(1:fourier.NFFT/2+1),...
       fourier.psd_mean_failure(1:fourier.NFFT/2+1)+fourier.psd_SEM_failure(1:fourier.NFFT/2+1),...
       fourier.f,'r')
end
% set(gca,'XScale','log')
% set(gca,'YScale','log')
xlabel('Frequency (Hz)')
axis([1 100 -80 0])
title('Power spectral density')
ylabel('dBuV/Hz')

else
end
fourier=rmfield(fourier,'input');
fourier=rmfield(fourier,'baseline');
assignin('base', 'spectro', spectro);
assignin('base', 'fourier', fourier);

    
    
    
    
    
    