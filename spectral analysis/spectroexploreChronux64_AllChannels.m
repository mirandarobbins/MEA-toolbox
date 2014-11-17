function [data params spectro] = spectroexploreChronux64_AllChannels(data,params);
params.flags.plot_online=1;
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=0;
end
%% Calculate spectrogram on a swep-by-sweep basis for all channels
spectro=[];
spectro.analysis_window=[params.search_win(1),40;params.search_win(2),60]; %i.e.... min time , min freq ; max time ; max freq

spectro.params.fpass=[0 200];
spectro.params.W=5;   % bandwidth in Hz
spectro.params.T=0.01; % duration to calculate tapers over (in s)
spectro.params.tapers=[spectro.params.W...
                       spectro.params.T...
                       0];
spectro.params.tapers=[1 1];
    
spectro.params.Fs=1000;
spectro.params.pad=2;
spectro.params.err=[2 1];
spectro.params.trialave=0;
spectro.params.movingwin=[0.05 0.01];


% chose modal most active channel
for sweep_id=data.sweep_sort.successful_sweeps
    spectro.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
end; clear sweep_id
spectro.StrongestChannel(spectro.StrongestChannel==0)=NaN;
       spectro.channeltoanalyse=mode(spectro.StrongestChannel);

% Run chronux spectrogram
for channel_id=1:64
   disp(strcat('analysing channel... ',num2str(channel_id)))
   [spectro.S{channel_id},...
     spectro.T,...
     spectro.F,...
     spectro.Serr{channel_id}] =  mtspecgramc(downsample(squeeze(data.filtered_lfp(:,channel_id,:)),20),...
                                       spectro.params.movingwin,...
                                       spectro.params);

% normalise by baseline power spectra    
for sweep_id=1:size(spectro.S{1},3)
    for tb_id=1:size(spectro.S{1},1)
    spectro.S{channel_id}(tb_id,:,sweep_id)=spectro.S{channel_id}(tb_id,:,sweep_id)-mean(spectro.S{channel_id}(1:8,:,sweep_id)) ;
    end    
end
   
% convert to dB
spectro.P{channel_id}    = 10*log10(abs(spectro.S{channel_id})); 
temp=spectro.P{channel_id}(:,:,data.sweep_sort.successful_sweeps);
spectro.mean{channel_id} = max(temp,[],3); %mean(temp,3); 
spectro.mean{channel_id} =smooth2a(spectro.mean{channel_id},1,1);
spectro.SD{channel_id}   = flipud(std(temp,0,3)); clear temp
end
% means
   
%     spectro.mean_smooth=smooth2a(spectro.mean,1, 20); %2D smoothing
%     spectro.SD_smooth=smooth2a(spectro.SD,1, 20); %2D smoothing
%     temp=mean(spectro.mean(:,params.baseline_win_samples(1):params.baseline_win_samples(2)),2);
%     spectro.mean_norm=spectro.mean-repmat(temp,1,size(spectro.mean,2)); clear temp
%     spectro.mean_norm_smooth=smooth2a(spectro.mean_norm,20, 200); %2D smoothing
% imagesc(flipud(squeeze(spectro.P{channel_id}(:,:,6))'))

clear temp 
%% trial to trial aveage power and variability
t_temp=floor(spectro.T*1000);
spectro.power_ROI=[20 20];
for channel_id=1:64
%centre around burst minima time...
spectro.average_window_power{channel_id}=zeros(numel(spectro.F),numel(data.sweep_sort.successful_sweeps));
spectro.average_window_power{channel_id}(spectro.average_window_power{channel_id}==0)=NaN;


for trial_id=data.sweep_sort.successful_sweeps
    this_latency=round(data.burst_timing.latency{trial_id}(channel_id));
     [~, array_position] = min(abs(t_temp - this_latency));
     spectro.burst_centre{channel_id}(trial_id)=spectro.T(array_position);

%%%%%%for window that moves to centre of burst (minima)
    time_limits=array_position-spectro.power_ROI(1)/10:array_position+spectro.power_ROI(2)/10;
%     spectro.average_window_power{channel_id}(:,trial_id)=max(spectro.P{channel_id}...
%                                                                 (time_limits,:,trial_id),[],1);
%%%%%%for fixed window
    spectro.average_window_power{channel_id}(:,trial_id)=max(spectro.P{channel_id}(18:36,:,trial_id),[],1);
%     spectro.average_window_power{channel_id}(:,trial_id)=mean(spectro.P{channel_id}(18:36,:,trial_id),1);

    spectro.average_window_power{channel_id}(spectro.average_window_power{channel_id}==0)=NaN;

%%%%%%normalization options
    %divide by baseline spectrum
        
% spectro.average_window_power{channel_id}(:,trial_id)=spectro.average_window_power{channel_id}(:,trial_id)-spectro.average_window_power{channel_id}(1,trial_id);
%     spectro.average_window_power{channel_id}(:,trial_id)=spectro.average_window_power{channel_id}(:,trial_id)./max(spectro.P{channel_id}(1:5,:,trial_id),[],1)';

end
spectro.average_window_power_mean{channel_id}=nanmean(spectro.average_window_power{channel_id},2);
spectro.average_window_power_sem{channel_id}=nansem(spectro.average_window_power{channel_id},2);
end
% clear t_temp array_position this_latency trial_id channel_id channel_to_plot
%% find peak gamma frequency band
gamma_range_coords=9:22; % this is approx 31:82Hz
for channel_id=1:64
    for sweep_id=data.sweep_sort.successful_sweeps
        
        gamma_band_temp=smooth(spectro.average_window_power{channel_id}(gamma_range_coords,sweep_id),2);
%         plot(gamma_band_temp)
        [gamma_amp gamma_freq ]=findpeaks(gamma_band_temp);
        if gamma_freq~=1; found_peak=1; else found_peak=0; end
        if isempty(gamma_freq)
            gamma_freq=NaN;
            gamma_amp=NaN;
        else
        gamma_freq=gamma_freq +(gamma_range_coords(1)-1);gamma_freq = round(spectro.F(gamma_freq));
        gamma_freq=gamma_freq(1);
        gamma_amp=gamma_amp(1);
        end
        if isempty(gamma_freq)
            gamma_freq =NaN;
            gamma_amp  =NaN;
        end
        spectro.window_power.gamma_centre_freq(channel_id,sweep_id)=gamma_freq;
        spectro.window_power.gamma_centre_freq_power(channel_id,sweep_id)=gamma_amp;
        spectro.window_power.gamma_centre_freq_foundYN(channel_id,sweep_id)=found_peak;

    end
    spectro.window_power.gamma_centre_freq(spectro.window_power.gamma_centre_freq==0)=NaN;
    spectro.window_power.gamma_centre_freq_power(spectro.window_power.gamma_centre_freq_power==0)=NaN;
    spectro.window_power.gamma_centre_freq_mean(channel_id)=nanmean(spectro.window_power.gamma_centre_freq(channel_id,:));
    spectro.window_power.gamma_centre_freq_SEM(channel_id)=nansem(spectro.window_power.gamma_centre_freq(channel_id,:));
end
% errorbar(spectro.window_power.gamma_centre_freq_mean,spectro.window_power.gamma_centre_freq_SEM)
clear found_peak gamma_amp gamma_band_temp gamma_freq gamma_range_coords
%% optional plots
if plotYN==1
figure;
for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0.01, 'Margin', 0.01);  

        hold on,%box off; axis off
    plot(spectro.F,spectro.average_window_power_mean{channel_id},'-k','LineWidth',2);
    ciplot(spectro.average_window_power_mean{channel_id}+spectro.average_window_power_sem{channel_id},...
       spectro.average_window_power_mean{channel_id}-spectro.average_window_power_sem{channel_id},...
       spectro.F,'k');
   
   for trial_id=data.sweep_sort.successful_sweeps
        if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        plot(spectro.F,spectro.average_window_power{channel_id}(:,trial_id),'LineWidth',0.5,'color',[0.3 0.3 0.3]);
        else
        plot(spectro.F,spectro.average_window_power{channel_id}(:,trial_id),'LineWidth',0.5,'color',[0.3 0.3 0.3]);
        end
   end
%        set(gca,'XScale','log')
%        set(gca,'YScale','log')
%        axis square
axis([0 100 -80 -30])
end
% plot individual trials highlighting peak gamma freqs
figure; hold on
for trial_id=data.sweep_sort.successful_sweeps
             plot(spectro.F,spectro.average_window_power{spectro.channeltoanalyse}(:,trial_id),'k','LineWidth',0.5);
             scatter(spectro.window_power.gamma_centre_freq(spectro.channeltoanalyse,trial_id),spectro.window_power.gamma_centre_freq_power(spectro.channeltoanalyse,trial_id),'og');
             
end
xlabel('Frequency (Hz)');ylabel('Power (dB)')
axis([0 100 -80 -30])

end    
%%
for trial_id=data.sweep_sort.successful_sweeps
    for channel_id=1:64
        spectro.window_power.twenty{trial_id}(channel_id)=spectro.average_window_power{channel_id}(6,trial_id);
        spectro.window_power.fifty{trial_id}(channel_id)=spectro.average_window_power{channel_id}(14,trial_id);
        spectro.window_power.onehunderd{trial_id}(channel_id)=spectro.average_window_power{channel_id}(27,trial_id);
        spectro.window_power.twohunderd{trial_id}(channel_id)=spectro.average_window_power{channel_id}(52,trial_id);
    end
   spectro.window_power.twenty{trial_id}=reshape(spectro.window_power.twenty{trial_id},8,8); 
   spectro.window_power.fifty{trial_id}=reshape(spectro.window_power.fifty{trial_id},8,8); 
   spectro.window_power.onehunderd{trial_id}=reshape(spectro.window_power.onehunderd{trial_id},8,8); 
   spectro.window_power.twohunderd{trial_id}=reshape(spectro.window_power.twohunderd{trial_id},8,8); 
end
assignin('base', 'spectro', spectro);
