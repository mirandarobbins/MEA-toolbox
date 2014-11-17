function [data params spikes spectro] = shared_coherence(data,spikes,params,spectro)
params.flags.plot_online=1;
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=0;
end
%% choices for strongest channel
spectro.coherence_radius=[];
channel_index=1:64;    
%%% 1 choose mode stongest channel   
for sweep_id=data.sweep_sort.successful_sweeps
    meta.StrongestChannel(sweep_id)=...
    find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
         min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
         );
end;
meta.StrongestChannel(meta.StrongestChannel==0)=NaN;
chosen_channel=meta.StrongestChannel;
% chosen_channel=repmat(mode(meta.StrongestChannel),1,numel(data.sweep_sort.successful_sweeps)); % mode?
    
% for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
%         this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
% %% 2 - choose biggest for each trial individually
% chosen_channel(sweep_id)=channel_index(...
%                           min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))==...
%                           min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
% 
% %%% 3 - choose earliest for each trial individually                        
% % chosen_channel(sweep_id)=round(median(channel_index(data.burst_timing.latency{this_sweep}==min(min(data.burst_timing.latency{this_sweep})))));              
% end;
clear sweep_id this_sweep
%% find co-active channels
channel_index=1:64;
spectro.coherence_radius.threshold=0.7;
spectro.coherence_radius.activeYN= zeros(64,numel(data.sweep_sort.successful_sweeps));
spectro.coherence_radius.active_channels=cell(numel(data.sweep_sort.successful_sweeps),1);
    for trial_id=data.sweep_sort.successful_sweeps
%         this_trial=data.sweep_sort.successful_sweeps(trial_id);
        amp_temp=min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,trial_id),[],1); amp_temp=amp_temp/min(amp_temp);
        % hide strongest channel:
        amp_temp(amp_temp==1)=0;

        % option 1: all that exceed threshold?
        spectro.coherence_radius.activeYN(gt(amp_temp,spectro.coherence_radius.threshold),trial_id)=1;
        spectro.coherence_radius.active_channels{trial_id}= channel_index(gt(amp_temp,spectro.coherence_radius.threshold));
        
        % option 2: just 2nd largest?
%        spectro.coherence_radius.active_channels{trial_id}  =   channel_index(amp_temp==max(amp_temp));
        
    end

%     % option (2) cont.: use modal next strongest?
%     for trial_id=1:numel(data.sweep_sort.successful_sweeps)
%         spectro.coherence_radius.active_channels{trial_id} = mode(cell2mat(spectro.coherence_radius.active_channels));
%     end

clear trial_id this_trial amp_temp channel_id
%% extract stats

%%%%% choose peak gamma frequency, extract peak coherence and phase between co-active channels at this freq
for trial_id=data.sweep_sort.successful_sweeps
    this_trial=(trial_id);
    this_channel=chosen_channel(trial_id);
    this_freq=spectro.window_power.gamma_centre_freq(this_channel,this_trial); % throws NaN if couldn't successfully tag peak gamma... fudge to 50Hz (fix ASAP!)
        if isnan(this_freq); 
            this_freq = 50; 
        end
    this_freq2=spectro.coherence_LFPLFP.f(ceil(spectro.coherence_LFPLFP.f)==this_freq);
    this_freq_idx=1:numel(spectro.coherence_LFPLFP.f); 
    this_freq_idx=this_freq_idx(spectro.coherence_LFPLFP.f==this_freq2(1));
    
    for channel_id=1:numel(spectro.coherence_radius.active_channels{trial_id})
        active_channel=spectro.coherence_radius.active_channels{trial_id}(channel_id);
        % COMPLETE coherence/phase spectra across WHOLE ACTIVE CHANNEL RADIUS...
        % LFPLFP 
        spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence{this_trial}(:,channel_id)       =   spectro.coherence_LFPLFP.C{active_channel}(:,this_trial);
        spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase{this_trial}(:,channel_id)           =   spectro.coherence_LFPLFP.phi{active_channel}(:,this_trial);
        % LFPspike 
        spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence{this_trial}(:,channel_id)   =   spectro.coherence_LFPspike.C{active_channel}(:,this_trial);
        spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase{this_trial}(:,channel_id)       =   spectro.coherence_LFPspike.phi{active_channel}(:,this_trial);
        
        % Coherence/phase at PEAK GAMMA FREQ only across active channel radius...
        % LFPLFP 
        spectro.coherence_radius.LFPLFP.coherence_at_peak_gamma{this_trial}(channel_id)                 =   spectro.coherence_LFPLFP.C{active_channel}(this_freq_idx,this_trial);
        spectro.coherence_radius.LFPLFP.phase_at_peak_gamma{this_trial}(channel_id)                     =   spectro.coherence_LFPLFP.phi{active_channel}(this_freq_idx,this_trial);
        % LFPspike (coherence/phase) at peak gamma freq across all active channels
        spectro.coherence_radius.LFPspike.coherence_at_peak_gamma{this_trial}(channel_id)               =   spectro.coherence_LFPspike.C{active_channel}(this_freq_idx,this_trial);
        spectro.coherence_radius.LFPspike.phase_at_peak_gamma{this_trial}(channel_id)                   =   spectro.coherence_LFPspike.phi{active_channel}(this_freq_idx,this_trial);
    end
end
%%%%% Complete coherence/phase spectra across WHOLE ACTIVE CHANNEL RADIUS...
    % LFPLFP coherence stats
    spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_cat    =   cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence,'UniformOutput',0)');
    spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean   =   nanmean(spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_cat);
    spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_SEM    =   nansem(spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_cat);
    % LFPLFP phase stats
    spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_cat    =   cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase,'UniformOutput',0)');
    spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean   =   nanmean(spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_cat);
    spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_SEM    =   nansem(spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_cat);
    % LFPspike coherence stats
    spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_cat    =   cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence,'UniformOutput',0)');
    spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean   =   nanmean(spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_cat);
    spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_SEM    =   nansem(spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_cat);
    % LFPspike phase stats
    spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_cat    =   cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase,'UniformOutput',0)');
    spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean   =   nanmean(spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_cat);
    spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_SEM    =   nansem(spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_cat);
   
%%%%% Coherence/phase at PEAK GAMMA FREQ across active channel radius...
    % LFPLFP coherence stats
    spectro.coherence_radius.LFPLFP.coherence_at_peak_gamma_cat     =  cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPLFP.coherence_at_peak_gamma,'UniformOutput',0)');
    spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean       =  nanmean(spectro.coherence_radius.LFPLFP.coherence_at_peak_gamma_cat);
    spectro.coherence_radius.LFPLFP.peak_gamma_coherence_SEM        =  nansem(spectro.coherence_radius.LFPLFP.coherence_at_peak_gamma_cat);
    % LFPLFP phase stats
    spectro.coherence_radius.LFPLFP.phase_at_peak_gamma_cat         =  cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPLFP.phase_at_peak_gamma,'UniformOutput',0)');
    spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean           =  nanmean(spectro.coherence_radius.LFPLFP.phase_at_peak_gamma_cat);
    spectro.coherence_radius.LFPLFP.peak_gamma_phase_SEM            =  nansem(spectro.coherence_radius.LFPLFP.phase_at_peak_gamma_cat);
    
    % LFPspike coherence stats
    spectro.coherence_radius.LFPspike.coherence_at_peak_gamma_cat   =  cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPspike.coherence_at_peak_gamma,'UniformOutput',0)');
    spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean     =  nanmean(spectro.coherence_radius.LFPspike.coherence_at_peak_gamma_cat);
    spectro.coherence_radius.LFPspike.peak_gamma_coherence_SEM      =  nansem(spectro.coherence_radius.LFPspike.coherence_at_peak_gamma_cat);
    % LFPspike phase stats
    spectro.coherence_radius.LFPspike.phase_at_peak_gamma_cat       =  cell2mat(cellfun(@transpose,spectro.coherence_radius.LFPspike.phase_at_peak_gamma,'UniformOutput',0)');
    spectro.coherence_radius.LFPspike.peak_gamma_phase_mean         =  nanmean(spectro.coherence_radius.LFPspike.phase_at_peak_gamma_cat);
    spectro.coherence_radius.LFPspike.peak_gamma_phase_SEM          =  nansem(spectro.coherence_radius.LFPspike.phase_at_peak_gamma_cat);
clear trial_id this_trial this_channel this_freq this_freq2 this_freq_idx channel_id active_channel

%%
%% plots
if plotYN==1
% figure
% subplot(2,1,1)
% scatter(spectro.coherence_radius.LFPLFP_peak_gamma_phase_cat,spectro.coherence_radius.LFPLFP_peak_gamma_coherence_cat)
% subplot(2,1,2)
% scatter(spectro.coherence_radius.LFPspike_peak_gamma_phase_cat,spectro.coherence_radius.LFPspike_peak_gamma_coherence_cat)
% figure
% wind_rose(spectro.coherence_radius.LFPLFP_peak_gamma_phase_cat+180,spectro.coherence_radius.LFPLFP_peak_gamma_coherence_cat)
end
assignin('base', 'spectro', spectro);

