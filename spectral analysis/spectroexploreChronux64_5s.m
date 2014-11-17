function [data params spectro fourier] = spectroexploreChronux64_5s(data,params)
params.flags.plot_online=1;
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=0;
end
%% Calculate spectrogram on a swep-by-sweep basis for strongest channel
spectro=[];
spectro.analysis_window=[params.search_win(1),40;params.search_win(2),60]; %i.e.... min time , min freq ; max time ; max freq

spectro.params.fpass     =   [0 100];
spectro.params.W         =   10;   % bandwidth in Hz /10
spectro.params.T         =   0.2; % duration to calculate tapers over (in s)  /0.1
spectro.params.tapers    =   [spectro.params.W     spectro.params.T    0];
spectro.params.tapers    =   [0.1 2]; % /off <<<<< NORMAL
spectro.params.Fs        =   params.Fs;
spectro.params.pad       =   3;
spectro.params.err       =   [2 1];
spectro.params.trialave  =   0;
spectro.params.movingwin =   [0.1 0.01]; %<<<<< NORMAL
% spectro.params.movingwin =   [100 10];
for sweep_id=data.sweep_sort.successful_sweeps
    spectro.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
end; clear sweep_id
spectro.StrongestChannel(spectro.StrongestChannel==0)=NaN;
% spectro.channeltoanalyse=mode(spectro.StrongestChannel);
spectro.channeltoanalyse=20;

% Run chronux spectrogram
    
    [spectro.S,...
     spectro.T,...
     spectro.F,...
     spectro.Serr] =  mtspecgramc(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,:)),...
                                       spectro.params.movingwin,...
                                       spectro.params);
% for sweep_id=1:size(spectro.S,3)
%     for tb_id=1:size(spectro.S,1)
%     spectro.S(tb_id,:,sweep_id)=spectro.S(tb_id,:,sweep_id)-mean(spectro.S(1:5,:,sweep_id)) ;
%     end
% end
%     
    spectro.P    = 10*log10(abs(spectro.S)); %convert to dB

% means
    temp=spectro.P(:,:,data.sweep_sort.successful_sweeps);
    spectro.mean = max(temp,[],3); %mean(temp,3); 
%     spectro.mean =smooth2a(spectro.mean,1,1);
    spectro.SD   = flipud(std(temp,0,3)); clear temp
%     spectro.mean_smooth=smooth2a(spectro.mean,1, 20); %2D smoothing
%     spectro.SD_smooth=smooth2a(spectro.SD,1, 20); %2D smoothing
%     temp=mean(spectro.mean(:,params.baseline_win_samples(1):params.baseline_win_samples(2)),2);
%     spectro.mean_norm=spectro.mean-repmat(temp,1,size(spectro.mean,2)); clear temp
%     spectro.mean_norm_smooth=smooth2a(spectro.mean_norm,20, 200); %2D smoothing

clear temp 
% trial to trial aveage power and variability
t_temp=floor(spectro.T*1000);
% (1) Fixed window 
%     spectro.average_window_power=...
%         squeeze(mean(spectro.P(find(t_temp==params.search_win(1)):find(t_temp==params.search_win(2)),...
%                                                         :,data.sweep_sort.successful_sweeps)));
% 
%     spectro.average_window_power_mean   =   nanmean(spectro.average_window_power,2);
%     spectro.average_window_power_sem    =   nansem(spectro.average_window_power,2);
%     if ~isempty(data.sweep_sort.failures)
%             spectro.average_window_power_failures       =   squeeze(mean(spectro.P(find(t_temp==params.search_win(1)):find(t_temp==params.search_win(2)),...
%                                                                     :,data.sweep_sort.failures)));
%             spectro.average_window_power_failures_mean  =   nanmean(spectro.average_window_power_failures,2);
%             spectro.average_window_power_failures_sem   =   nansem(spectro.average_window_power_failures,2);
%     end

% (2) OR, centre around burst minima time...
%     spectro.average_window_power                =    zeros(numel(spectro.F),numel(data.sweep_sort.successful_sweeps));
%     spectro.average_window_power...
%         (spectro.average_window_power==0)       =    NaN;
%     spectro.power_ROI=[50 200];
%     for trial_id=data.sweep_sort.successful_sweeps
%         this_latency                                    =   round(data.burst_timing.latency{trial_id}(spectro.channeltoanalyse));
%         [min_difference, array_position]                =   min(abs(t_temp - this_latency)); clear min_difference
%         spectro.burst_centre(trial_id)                  =   spectro.T(array_position);
%         spectro.average_window_power(:,trial_id)        =   max(...
%                                                                     spectro.P(array_position-spectro.power_ROI(1)/5:array_position+spectro.power_ROI(2)/5,:,trial_id),...
%                                                                 [],1);
%         spectro.average_window_power... 
%             (spectro.average_window_power==0)           =   NaN;
% %         spectro.average_window_power(:,trial_id)        =   spectro.average_window_power(:,trial_id)-...
% %                                                                 spectro.average_window_power(5,trial_id);
%     end
%     spectro.average_window_power_mean   =   nanmean(spectro.average_window_power,2);
%     spectro.average_window_power_sem    =   nansem(spectro.average_window_power,2);    

clear t_temp array_position this_latency trial_id
% Plot

if plotYN==1
    figure;set(gcf,'position',[20 185 755 800])
    subplot(211); hold on
%         plot(data.tb,sgolayfilt(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,data.sweep_sort.successful_sweeps)),0,15),'Color',[0.8 0.8 0.8]); 
        plot(data.tb/1000,mean(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,data.sweep_sort.successful_sweeps)),2),'-k')
        axis([0 5 -0.2 0.2])
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
        x= repmat(spectro.T,size(spectro.mean,2),1);
        y=repmat(spectro.F,size(spectro.mean,1),1)';
    %     z=spectro.P(:,:,2)';
        z=spectro.mean';
%         z=smooth2a(z,1,1);
        ih=surf(x,y,z); view(2); axis tight;
        set(ih, 'edgecolor', 'none');
        set(ih, 'facecolor', 'interp'); 
        colormap(jet); 
        set(gca,'clim',[-80 -60])
        ax2a=gca;
        pos=get(ax2a,'pos');
        set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);
        pos=get(gca,'pos');
    %  set(gca,'YScale','log')
        hc=colorbar('location','eastoutside',...
                    'position',[pos(1)+pos(3)+0.01 pos(2) 0.03 pos(4)*0.5],...
                    'YAxisLocation','right');
        ylabel(hc,'dBuV^2/Hz')
        title('Power Spectral Density Estimate')
        ylabel('Frequency (Hz)');
        xlabel('Time (s)');  
        
%     subplot(313);hold on;
%         spectro.channeltoanalyse=spectro.channeltoanalyse;
%         % for sweep_id=data.sweep_sort.successful_sweeps       
%         %     plot(fourier.f,fourier.psd(:,spectro.channeltoanalyse,sweep_id),':b')
%         % end
%         % semilogx(fourier.f,fourier.psd_mean_success(:,spectro.channeltoanalyse),'b','LineWidth',1.5)
%         % ciplot(fourier.psd_mean_success(:,spectro.channeltoanalyse)-fourier.psd_SEM_success(:,spectro.channeltoanalyse),...
%         %        fourier.psd_mean_success(:,spectro.channeltoanalyse)+fourier.psd_SEM_success(:,spectro.channeltoanalyse),...
%         %        fourier.f,'b')
%         % if ~isempty(data.sweep_sort.failures)
%         %     for sweep_id=data.sweep_sort.failures       
%         %         plot(fourier.f,fourier.psd(:,spectro.channeltoanalyse,sweep_id),':r')
%         %     end   
%         %     semilogx(fourier.f,fourier.psd_mean_failure(:,spectro.channeltoanalyse),'r','LineWidth',1.5)
%         %     ciplot(fourier.psd_mean_failure(:,spectro.channeltoanalyse)-fourier.psd_SEM_failure(:,spectro.channeltoanalyse),...
%         %            fourier.psd_mean_failure(:,spectro.channeltoanalyse)+fourier.psd_SEM_failure(:,spectro.channeltoanalyse),...
%         %            fourier.f,'r')
%         % end
%         % plot(spectro.F,spectro.average_window_power)%,':b')
%         for trial_id=data.sweep_sort.successful_sweeps
%                 if data.sweep_sort.window_std(spectro.channeltoanalyse,trial_id)>params.detection_threshold
%                     plot(spectro.F,spectro.average_window_power(:,trial_id),'r','LineWidth',0.5);
%                 else
%                     plot(spectro.F,spectro.average_window_power(:,trial_id),':','LineWidth',0.5,'color',[0.3 0.3 0.3]);
%                 end
%         end
%         semilogx(spectro.F,spectro.average_window_power_mean,'b','LineWidth',1.5)
%         ciplot(spectro.average_window_power_mean+spectro.average_window_power_sem,...
%                spectro.average_window_power_mean-spectro.average_window_power_sem,...
%                spectro.F,'b')
%         pos=get(gca,'pos');
%         set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);
% 
%         % if ~isempty(data.sweep_sort.failures)
%         % plot(spectro.F,spectro.average_window_power_failures,':r')
%         % semilogx(spectro.F,spectro.average_window_power_failures_mean,'r','LineWidth',1.5)
%         % ciplot(spectro.average_window_power_failures_mean+spectro.average_window_power_failures_sem,...
%         %        spectro.average_window_power_failures_mean-spectro.average_window_power_failures_sem,...
%         %        spectro.F,'r')
%         % end
%         % set(gca,'XScale','log')
%         % set(gca,'YScale','log')
%         xlabel('Frequency (Hz)')
%         axis([1 spectro.params.fpass(2) -Inf 0])
%         title('Power spectral density')
%         ylabel('dBuV/Hz')
end
%% Run chronux multi-taper FFT
fourier=[];

for trial_id=1%:params.last_sweep% size(data.filtered_lfp,3)   
    fourier.input       = data.filtered_lfp;%(params.search_win_samples(1):params.search_win_samples(2),:,trial_id);%data.sweep_sort.successful_sweeps);
    fourier.T           = 1/params.Fs;                      % Sample time
    fourier.L           = size(fourier.input,1);            % Length of signal
%     fourier.baseline    = data.filtered_lfp(params.baseline_win_samples(1):params.baseline_win_samples(2),:,trial_id);%data.sweep_sort.successful_sweeps);
%     fourier.T_baseline  = 1/params.Fs;                      % Sample time
%     fourier.L_baseline  = size(fourier.baseline,1);         % Length of signal
    
    [fourier.fft(:,:,trial_id),fourier.f]               =    mtspectrumc(fourier.input,spectro.params);    
    fourier.psd(:,:,trial_id)                           =    10*log10(abs(fourier.fft(:,:,trial_id)))   ;%/(fourier.L/2));
%     [fourier.fft_baseline(:,:,trial_id),fourier.g]      =    mtspectrumc(fourier.baseline,spectro.params);
%     fourier.psd_baseline(:,:,trial_id)                  =    10*log10(abs(fourier.fft_baseline(:,:,trial_id)))  ;%/(fourier.L_baseline/2));
    
%     fourier.psd_norm(:,:,trial_id)   =  abs(fourier.psd(:,:,trial_id))./...
%                                         nanmean(abs(fourier.psd_baseline(:,:,trial_id)),3); %  normalisation here?
end; 
% fourier=rmfield(fourier,'g'); 
clear trial_id
fourier.fft_mean_success = nanmean(fourier.fft(:,:,data.sweep_sort.successful_sweeps),3);
fourier.fft_SEM_success  = nansem(fourier.fft(:,:,data.sweep_sort.successful_sweeps),3);
% fourier.fft_mean_failure = nanmean(fourier.fft(:,:,data.sweep_sort.failures),3);
% fourier.fft_SEM_failure  = nansem(fourier.fft(:,:,data.sweep_sort.failures),3);
fourier.psd_mean_success = nanmean(fourier.psd(:,:,data.sweep_sort.successful_sweeps),3);
fourier.psd_SEM_success  = nansem(fourier.psd(:,:,data.sweep_sort.successful_sweeps),3);
% fourier.psd_mean_failure = nanmean(fourier.psd(:,:,data.sweep_sort.failures),3);
% fourier.psd_SEM_failure  = nansem(fourier.psd(:,:,data.sweep_sort.failures),3);

if plotYN==1        
    figure; hold on
    plot(fourier.f,fourier.psd_mean_success(:,spectro.channeltoanalyse),'k')
%     ciplot(fourier.psd_mean_success(:,spectro.channeltoanalyse) - fourier.psd_SEM_success(:,spectro.channeltoanalyse),...
%            fourier.psd_mean_success(:,spectro.channeltoanalyse) + fourier.psd_SEM_success(:,spectro.channeltoanalyse),...
%            fourier.f,'k')
    xlabel('Frequency (Hz)')
%     axis([1 spectro.params.fpass(2) -Inf -40])
    title('Power spectral density')
    ylabel('dBuV/Hz')
end;
clear x y z ax1 ax2a ax2b ih hc pos

%%
fourier=rmfield(fourier,'input');
% fourier=rmfield(fourier,'baseline');
clear sweep_id 


assignin('base', 'spectro', spectro);
assignin('base', 'fourier', fourier);


    