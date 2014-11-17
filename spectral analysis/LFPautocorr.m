function [spectro] = LFPautocorr(data, params)
params.flags.plot_online=0;
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=0;
end
%% runs autocorrelation of LFP in window
spectro.corr=[];
spectro.corr.maxlags= 50;
for channel_id=1:64
    for trial_id=1:params.last_sweep
        
        temp=downsample(data.filtered_spikes(params.search_win_samples(1):params.search_win_samples(2),channel_id,trial_id),20); %
%         temp=data.filtered_spikes(params.search_win_samples(1):params.search_win_samples(2)-1,channel_id,trial_id); %

        temp=temp-mean(temp); temp=temp/std(temp);
        spectro.corr.LFP_autocorr{channel_id}(trial_id,:)=xcov(temp,temp,spectro.corr.maxlags,'unbiased');
    end
    spectro.corr.LFP_autocorr_mean(channel_id,:)=nanmean(spectro.corr.LFP_autocorr{channel_id}(data.sweep_sort.successful_sweeps,:));
    spectro.corr.LFP_autocorr_sem(channel_id,:)=nansem(spectro.corr.LFP_autocorr{channel_id}(data.sweep_sort.successful_sweeps,:));
end
%% choices for strongest channel
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
chosen_channel=repmat(mode(meta.StrongestChannel),1,numel(meta.StrongestChannel)); % mode?
    
% for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
%     this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
% % 2 - choose biggest for each trial individually
% chosen_channel(sweep_id)=channel_index(...
%                           min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))==...
%                           min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
% 
% %%% 3 - choose earliest for each trial individually                        
% % chosen_channel(sweep_id)=round(median(channel_index(data.burst_timing.latency{this_sweep}==min(min(data.burst_timing.latency{this_sweep})))));              
% end;
clear sweep_id this_sweep
%% extract centre channel only
spectro.corr.centre.lags=-spectro.corr.maxlags:spectro.corr.maxlags;
spectro.corr.centre.LFP_autocorr=zeros(params.last_sweep,numel(spectro.corr.centre.lags));
spectro.corr.centre.LFP_autocorr(spectro.corr.centre.LFP_autocorr==0)=NaN;
for sweep_id=data.sweep_sort.successful_sweeps
    spectro.corr.centre.LFP_autocorr(sweep_id,:) = spectro.corr.LFP_autocorr{chosen_channel(sweep_id)}(sweep_id,:);
end
spectro.corr.centre.LFP_autocorr_mean   =   nanmean(spectro.corr.centre.LFP_autocorr,1);
spectro.corr.centre.LFP_autocorr_SEM    =   nansem(spectro.corr.centre.LFP_autocorr,1);
%% plot
if plotYN==1
    figure
    for channel_id=1:64
            subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); 
    hold on
           plot(-spectro.corr.maxlags:spectro.corr.maxlags,spectro.corr.LFP_autocorr{channel_id}(:,:),'-r','LineWidth',1);
            plot(-spectro.corr.maxlags:spectro.corr.maxlags,spectro.corr.LFP_autocorr_mean(channel_id,:),'-k','LineWidth',2);
            ciplot(spectro.corr.LFP_autocorr_mean(channel_id,:)+spectro.corr.LFP_autocorr_sem(channel_id,:),...
                spectro.corr.LFP_autocorr_mean(channel_id,:)-spectro.corr.LFP_autocorr_sem(channel_id,:),...
                -spectro.corr.maxlags:spectro.corr.maxlags,'k');

    axis([-1*spectro.corr.maxlags spectro.corr.maxlags 0 1])
    %        set(gca,'YScale','log')
    end
    figure; hold on
        plot(spectro.corr.centre.lags,spectro.corr.centre.LFP_autocorr_mean,'-k','LineWidth',2);
        ciplot(spectro.corr.centre.LFP_autocorr_mean - spectro.corr.centre.LFP_autocorr_SEM,...
               spectro.corr.centre.LFP_autocorr_mean + spectro.corr.centre.LFP_autocorr_SEM,...
               spectro.corr.centre.lags,'k');
end
assignin('base', 'spectro', spectro);
