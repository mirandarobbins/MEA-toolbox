function [data params spikes spectro] = MEAChannel_coherence(data,spikes,params,spectro)
params.flags.plot_online=1;
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=0;
end
%% params for LFP-spike coherence
spectro.params=[];
spectro.params.Fs=1000;
spectro.params.fpass=[0 500];
spectro.params.W=10;                % bandwidth in Hz (default 10)
spectro.params.T=1;                 % duration to calculate tapers over (in s)
spectro.params.tapers=[spectro.params.W,    spectro.params.T,   0];
spectro.params.tapers=[20 4]; % [20 4]
spectro.params.pad=1;
spectro.params.err=[2 1];
spectro.params.trialave=0;
% spectro.params.movingwin=[spectro.params.T spectro.params.T/10];
% spectro.params.win=1;
spectro.params.fscorr=0;

%%%%%%% convert spike times to seconds, pad array
% temp_spike_times2=cell(1,64);
% temp_spike_times=cell(64,params.last_sweep);
%     for channel_id=1:64
%         temp_spike_times2{channel_id}=zeros(1000,params.last_sweep);
%     end
% for channel_id=1:64
%     for trial_id = 1:params.last_sweep
%         temp_spike_times{channel_id,trial_id}=spikes.spiketimes{channel_id,trial_id}./100; % convert to seconds
%         temp=temp_spike_times{channel_id,trial_id};
% %         temp=temp*20000; % convert back to samples (@ 20KHz)
%         temp_spike_times2{channel_id}(1:numel(temp),trial_id)=temp;
%         temp_spike_times2{channel_id}(temp_spike_times2{channel_id}==0)=NaN;
% %         temp_spike_times2{channel_id}=round(temp_spike_times2{channel_id}*1000);
%         if ~isempty(temp_spike_times2{channel_id})
%         end
%     end
% end
% spectro.coherence_LFPspike.input_spikes=temp_spike_times2;
% clear temp_spike_times temp_spike_times2 temp

%LFP-spike coherence
spectro.coherence_LFPspike=[];
for channel_id=1:64
% for trial_id = 1:params.last_sweep
%     inputa = squeeze(data.filtered_lfp(:,channel_id,trial_id));
%     inputb = spectro.coherence_LFPspike.input_spikes{channel_id}(:,trial_id);
%     inputb(isnan(inputb))=[];

    inputa=downsample(squeeze(data.filtered_lfp(:,channel_id,:)),20);
    inputb=spikes.PDF.timestamps{channel_id}(:,:)';    
    inputa=inputa(200:1000,:);
    inputb=inputb(200:1000,:);
    
%     inputa=inputa(params.search_win(1):params.search_win(2),:);
%     inputb=inputb(params.search_win(1):params.search_win(2),:);
%     inputb(isnan(inputb))=0;

%     [spectro.coherence_LFPspike.C{channel_id}(:,trial_id),...
%      spectro.coherence_LFPspike.phi{channel_id}(:,trial_id),...
%      spectro.coherence_LFPspike.S12{channel_id}(:,trial_id),...
%      spectro.coherence_LFPspike.S1{channel_id}(:,trial_id),...
%      spectro.coherence_LFPspike.S2{channel_id}(:,trial_id),...
%      spectro.coherence_LFPspike.f,...
%      spectro.coherence_LFPspike.zerosp{channel_id}(:,trial_id),...
%      spectro.coherence_LFPspike.confC{channel_id}(trial_id),...
%      spectro.coherence_LFPspike.phistd,...
%      spectro.coherence_LFPspike.Cerr]=coherencycpb(...
%                                           inputa,...
%                                           inputb,...
%                                           spectro.params,...
%                                           spectro.params.fscorr);
%   
 
 
    [spectro.coherence_LFPspike.C{channel_id},...
     spectro.coherence_LFPspike.phi{channel_id},...
     spectro.coherence_LFPspike.S12{channel_id},...
     spectro.coherence_LFPspike.S1{channel_id},...
     spectro.coherence_LFPspike.S2{channel_id},...
     spectro.coherence_LFPspike.f,...
     spectro.coherence_LFPspike.zerosp{channel_id},...
     spectro.coherence_LFPspike.confC{channel_id},...
     spectro.coherence_LFPspike.phistd{channel_id},...
     spectro.coherence_LFPspike.Cerr{channel_id}]=coherencycpb(...
                                          inputa,...
                                          inputb,...
                                          spectro.params,...
                                          spectro.params.fscorr);

 
% end

     spectro.coherence_LFPspike.S1{channel_id}   =10*log10(abs(spectro.coherence_LFPspike.S1{channel_id})); %convert to dB

     spectro.coherence_LFPspike.S2{channel_id}   =10*log10(abs(spectro.coherence_LFPspike.S2{channel_id})); %convert to dB
%      for trial_id=1:max(data.sweep_sort.successful_sweeps)
%      spectro.coherence_LFPspike.S1{channel_id}(:,trial_id)   =-1*spectro.coherence_LFPspike.S1{channel_id}(:,trial_id)./...
%                                                                  spectro.coherence_LFPspike.S1{channel_id}(1,trial_id)+1;
%      spectro.coherence_LFPspike.S2{channel_id}(:,trial_id)   =spectro.coherence_LFPspike.S2{channel_id}(:,trial_id)./...
%                                                                  spectro.coherence_LFPspike.S2{channel_id}(1,trial_id)-1;
%      end
     
    spectro.coherence_LFPspike.C{channel_id}  (:,var(spectro.coherence_LFPspike.C{channel_id})==0)      = NaN;
    spectro.coherence_LFPspike.phi{channel_id}(:,var(spectro.coherence_LFPspike.phi{channel_id})==0)	= NaN;
    spectro.coherence_LFPspike.S12{channel_id}(:,var(spectro.coherence_LFPspike.S12{channel_id})==0)    = NaN;
    spectro.coherence_LFPspike.S1{channel_id}(:,var(spectro.coherence_LFPspike.S1{channel_id})==0)      = NaN;
    spectro.coherence_LFPspike.S2{channel_id}(:,isnan(var(spectro.coherence_LFPspike.S2{channel_id})))  = NaN;
    spectro.coherence_LFPspike.S1{channel_id}(:,isnan(spectro.coherence_LFPspike.S2{channel_id}(1,:)))  = NaN; % try this?

    
    spectro.coherence_LFPspike.C_mean(:,channel_id)=nanmean(spectro.coherence_LFPspike.C{channel_id},2);
    spectro.coherence_LFPspike.C_sem(:,channel_id,:)=nansem(spectro.coherence_LFPspike.C{channel_id},2);
    
    spectro.coherence_LFPspike.phi_mean(:,channel_id)=nanmean(abs(spectro.coherence_LFPspike.phi{channel_id}),2);
    spectro.coherence_LFPspike.phi_sem(:,channel_id)=nansem(abs(spectro.coherence_LFPspike.phi{channel_id}),2);
    
    spectro.coherence_LFPspike.S12_mean(:,channel_id)=nanmean(abs(spectro.coherence_LFPspike.S12{channel_id}),2);
    spectro.coherence_LFPspike.S12_sem(:,channel_id)=nansem(abs(spectro.coherence_LFPspike.S12{channel_id}),2);
    
    spectro.coherence_LFPspike.S1_mean(:,channel_id)=nanmean((spectro.coherence_LFPspike.S1{channel_id}),2);
    spectro.coherence_LFPspike.S1_sem(:,channel_id)=nansem((spectro.coherence_LFPspike.S1{channel_id}),2);
    
    spectro.coherence_LFPspike.S2_mean(:,channel_id)=nanmean(spectro.coherence_LFPspike.S2{channel_id},2);
    spectro.coherence_LFPspike.S2_sem(:,channel_id)=nansem(spectro.coherence_LFPspike.S2{channel_id},2);    
end
clear channel_id inputa inputb
%plot results
if plotYN==1
%%%%% Average spread of coherence_spikeLFP
    figure('name','mean LFP-spike coherence')        
    for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
        hold on,%box off; axis off
        % for trial_id=data.sweep_sort.successful_sweeps
        %     if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        %        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.C{channel_id}(:,trial_id),'r','LineWidth',0.5);
        %     else
        %        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.C{channel_id}(:,trial_id),'k','LineWidth',0.5);
        %     end
        % end
        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.C_mean(:,channel_id),'-k','LineWidth',1);
        ciplot(spectro.coherence_LFPspike.C_mean(:,channel_id)+spectro.coherence_LFPspike.C_sem(:,channel_id),...
            spectro.coherence_LFPspike.C_mean(:,channel_id)-spectro.coherence_LFPspike.C_sem(:,channel_id),...
            spectro.coherence_LFPspike.f,'k');
        % set(gca,'XScale','log')
        % set(gca,'YScale','log')
        % axis square
        axis([spectro.params.fpass(1) spectro.params.fpass(2) 0 1]) 
        if channel_id==spectro.channeltoanalyse
            set(gca,'color',[0 1 0])
        end
    end
%%%%% Plot average individual spectra - spikes and LFP    
    figure('name','mean LFP/spike power specta')      
    for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
        hold on, %=box off; axis off
        % for trial_id=data.sweep_sort.successful_sweeps
        %     if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        %        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.S2{channel_id}(:,trial_id),'r','LineWidth',0.5);
        %     else
        %        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.S2{channel_id}(:,trial_id),'k','LineWidth',0.5);
        %     end
        % end
        % LFP in red
        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.S1_mean(:,channel_id),'-r','LineWidth',1);
        ciplot(spectro.coherence_LFPspike.S1_mean(:,channel_id)+spectro.coherence_LFPspike.S1_sem(:,channel_id),...
            spectro.coherence_LFPspike.S1_mean(:,channel_id)-spectro.coherence_LFPspike.S1_sem(:,channel_id),...
            spectro.coherence_LFPspike.f,'r');
        % spikes in blue
        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.S2_mean(:,channel_id),'-b','LineWidth',1);
        ciplot(spectro.coherence_LFPspike.S2_mean(:,channel_id)+spectro.coherence_LFPspike.S2_sem(:,channel_id),...
            spectro.coherence_LFPspike.S2_mean(:,channel_id)-spectro.coherence_LFPspike.S2_sem(:,channel_id),...
            spectro.coherence_LFPspike.f,'b');
        % set(gca,'XScale','log')
        % set(gca,'YScale','log')
        % axis square
        axis([spectro.params.fpass(1) spectro.params.fpass(2) -100 0]) 
        if channel_id==spectro.channeltoanalyse
            set(gca,'color',[0 1 0])
        end
    end
%%%%% Plot mean phase    
    figure('name','mean LFP-spike phase')            
    for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
        hold on,%box off; axis off
        % for trial_id=data.sweep_sort.successful_sweeps
        %     if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        %        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.phi{channel_id}(:,trial_id)*180/pi,'r','LineWidth',0.5);
        %     else
        %        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.phi{channel_id}(:,trial_id)*180/pi,'k','LineWidth',0.5);
        %     end
        % end
        plot(spectro.coherence_LFPspike.f,spectro.coherence_LFPspike.phi_mean(:,channel_id)*180/pi,'-b','LineWidth',2);
        ciplot(spectro.coherence_LFPspike.phi_mean(:,channel_id)*180/pi+spectro.coherence_LFPspike.phi_sem(:,channel_id)*180/pi,...
            spectro.coherence_LFPspike.phi_mean(:,channel_id)*180/pi-spectro.coherence_LFPspike.phi_sem(:,channel_id)*180/pi,...
            spectro.coherence_LFPspike.f,'b');
        % set(gca,'XScale','log')
        % set(gca,'YScale','log')
        % axis square
        axis([spectro.params.fpass(1) spectro.params.fpass(2) 0 360]) 
        if channel_id==spectro.channeltoanalyse
            set(gca,'color',[0 1 0])
        end
    end
end

%% choose strongest channel to run against
for sweep_id=data.sweep_sort.successful_sweeps
    spectro.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
end; clear sweep_id
spectro.StrongestChannel(spectro.StrongestChannel==0)=NaN;
% spectro.channeltoanalyse=mode(spectro.StrongestChannel);% mode?
spectro.channeltoanalyse=spectro.StrongestChannel;
%% LFP-LFP coherence
spectro.params=[];
    spectro.params.Fs=1000;
    spectro.params.fpass=[0 500];   % spectro.params.Fs/2];
    spectro.params.W=10;            % bandwidth in Hz
    spectro.params.T=1;             % duration to calculate tapers over (in s)
    spectro.params.tapers=[spectro.params.W spectro.params.T,   2*spectro.params.T*spectro.params.W];
    spectro.params.tapers=[3 5];%[ 3 5]
    spectro.params.pad=1;
    % spectro.params.err=[2 1];
    spectro.params.trialave=0;
    
% compute segmented LFP coherence_LFPLFP between channels
spectro.coherence_LFPLFP=[];
for channel_id=1:64
    for trial_id = data.sweep_sort.successful_sweeps
        disp(strcat('processing Trial',num2str(trial_id),' Channel',num2str(channel_id)))
%         inputa = squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,spectro.channeltoanalyse,trial_id));
%         inputb = squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,channel_id,trial_id));

%         inputa=downsample(squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,spectro.channeltoanalyse,trial_id)),1);
%         inputb=downsample(squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,channel_id,trial_id)),1);
          inputa=downsample(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse(trial_id),trial_id)),20);
          inputb=downsample(squeeze(data.filtered_lfp(:,channel_id,trial_id)),20);
          inputa=zscore(inputa(200:1000));
          inputb=zscore(inputb(200:1000)); 
    [spectro.coherence_LFPLFP.C{channel_id}(:,trial_id),...
     spectro.coherence_LFPLFP.phi{channel_id}(:,trial_id),...
     spectro.coherence_LFPLFP.S12{channel_id}(:,trial_id),...
     spectro.coherence_LFPLFP.S1(:,trial_id),...
     spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id),...
     spectro.coherence_LFPLFP.f]=coherencyc(...
                                          inputa,...
                                          inputb,...
                                          spectro.params);
     
%      spectro.coherence_LFPLFP.S1(:,trial_id)=spectro.coherence_LFPLFP.S1(:,trial_id)./spectro.coherence_LFPLFP.S1(1,trial_id)*-1; % normalize to DC  
     spectro.coherence_LFPLFP.S1(:,trial_id)   =10*log10(abs(spectro.coherence_LFPLFP.S1(:,trial_id))); % convert to dB
    
%      spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id)=spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id)./spectro.coherence_LFPLFP.S2{channel_id}(1,trial_id)   *-1; % normalize to DC
     spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id) = 10*log10(abs(spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id))); % convert to dB
    end
    spectro.coherence_LFPLFP.C{channel_id}  (:,var(spectro.coherence_LFPLFP.C{channel_id})==0)      =   NaN;
    spectro.coherence_LFPLFP.phi{channel_id}(:,var(spectro.coherence_LFPLFP.phi{channel_id})==0)    =   NaN;
    spectro.coherence_LFPLFP.S12{channel_id}(:,var(spectro.coherence_LFPLFP.S12{channel_id})==0)    =   NaN;
    spectro.coherence_LFPLFP.S1             (:,var(spectro.coherence_LFPLFP.S1)==0)                 =   NaN;
    spectro.coherence_LFPLFP.S2{channel_id} (:,var(spectro.coherence_LFPLFP.S2{channel_id})==0)     =   NaN;

    
    spectro.coherence_LFPLFP.C_mean(:,channel_id)   = nanmean(spectro.coherence_LFPLFP.C{channel_id},2);
    spectro.coherence_LFPLFP.C_sem(:,channel_id,:)  = nansem(spectro.coherence_LFPLFP.C{channel_id},2);
    
    spectro.coherence_LFPLFP.phi_mean(:,channel_id) = nanmean(abs(spectro.coherence_LFPLFP.phi{channel_id}),2);
    spectro.coherence_LFPLFP.phi_sem(:,channel_id)  = nansem(abs(spectro.coherence_LFPLFP.phi{channel_id}),2);
    
    spectro.coherence_LFPLFP.S12_mean(:,channel_id) = nanmean(abs(spectro.coherence_LFPLFP.S12{channel_id}),2);
    spectro.coherence_LFPLFP.S12_sem(:,channel_id)  = nansem(abs(spectro.coherence_LFPLFP.S12{channel_id}),2);
    
    spectro.coherence_LFPLFP.S1_mean                = nanmean(abs(spectro.coherence_LFPLFP.S1),2);
    spectro.coherence_LFPLFP.S1_sem                 = nansem(abs(spectro.coherence_LFPLFP.S1),2);
    
    spectro.coherence_LFPLFP.S2_mean(:,channel_id)  = nanmean(spectro.coherence_LFPLFP.S2{channel_id},2);
    spectro.coherence_LFPLFP.S2_sem(:,channel_id)   = nansem(spectro.coherence_LFPLFP.S2{channel_id},2);    
end
% Construct coherency phase histograms

% Plot results
if plotYN==1
%%%%%%  Average spread of coherence_LFPLFP
    figure('name','mean LFP-LFP coherence')
    for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
        hold on,%box off; axis off
        % for trial_id=data.sweep_sort.successful_sweeps
        %     if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        %        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.C{channel_id}(:,trial_id),'r','LineWidth',0.5);
        %     else
        %        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.C{channel_id}(:,trial_id),'k','LineWidth',0.5);
        %     end
        % end
        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.C_mean(:,channel_id),'-k','LineWidth',1);
        ciplot(spectro.coherence_LFPLFP.C_mean(:,channel_id)+spectro.coherence_LFPLFP.C_sem(:,channel_id),...
            spectro.coherence_LFPLFP.C_mean(:,channel_id)-spectro.coherence_LFPLFP.C_sem(:,channel_id),...
            spectro.coherence_LFPLFP.f,'k');
        % set(gca,'XScale','log')
        % set(gca,'YScale','log')
        % axis square
        axis([spectro.params.fpass(1) spectro.params.fpass(2) 0 1]) 
        if channel_id==spectro.channeltoanalyse
            set(gca,'color',[0 1 0])
        end
    end
    
%%%%%%  Plot average individual spectra
    figure('name','mean LFP power spectra')
    for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
        hold on,%box off; axis off
        % for trial_id=data.sweep_sort.successful_sweeps
        %     if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        %        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id),'r','LineWidth',0.5);
        %     else
        %        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id),'k','LineWidth',0.5);
        %     end
        % end
        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.S2_mean(:,channel_id),'-k','LineWidth',1);
        ciplot(spectro.coherence_LFPLFP.S2_mean(:,channel_id)+spectro.coherence_LFPLFP.S2_sem(:,channel_id),...
            spectro.coherence_LFPLFP.S2_mean(:,channel_id)-spectro.coherence_LFPLFP.S2_sem(:,channel_id),...
            spectro.coherence_LFPLFP.f,'k');
        % set(gca,'XScale','log')
        % set(gca,'YScale','log')
        % axis square
        if channel_id==spectro.channeltoanalyse
            set(gca,'color',[0 1 0])
        end
        axis([spectro.params.fpass(1) spectro.params.fpass(2) -100 0]) 
    end
    
%%%%%%  Plot mean phase
    figure('name','mean LFP-LFP phase')
    for channel_id=1:64
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
        hold on,%box off; axis off
        % for trial_id=data.sweep_sort.successful_sweeps
        %     if data.sweep_sort.window_std(channel_id,trial_id)>params.detection_threshold
        %        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id),'r','LineWidth',0.5);
        %     else
        %        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id),'k','LineWidth',0.5);
        %     end
        % end
        plot(spectro.coherence_LFPLFP.f,spectro.coherence_LFPLFP.phi_mean(:,channel_id)*180/pi,'-b','LineWidth',2);
        ciplot(spectro.coherence_LFPLFP.phi_mean(:,channel_id)*180/pi+spectro.coherence_LFPLFP.phi_sem(:,channel_id)*180/pi,...
            spectro.coherence_LFPLFP.phi_mean(:,channel_id)*180/pi-spectro.coherence_LFPLFP.phi_sem(:,channel_id)*180/pi,...
            spectro.coherence_LFPLFP.f,'b');
        % set(gca,'XScale','log')
        % set(gca,'YScale','log')
        % axis square
        if channel_id==spectro.channeltoanalyse
            set(gca,'color',[0 1 0])
        end
        axis([spectro.params.fpass(1) spectro.params.fpass(2) 0 360]) 
    end
clear inputa inputb
end
%%
assignin('base', 'spectro', spectro);

