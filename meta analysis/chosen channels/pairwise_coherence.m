function [data params spikes spectro] = pairwise_coherence(data,spikes,params,L4,L23)
%% LFP-spike coherence
% params for LFP-spike coherence
    spectro.params=[];
    spectro.params.Fs=1000;
    spectro.params.fpass=[0 500];
    spectro.params.W=10;                % bandwidth in Hz (default 10)
    spectro.params.T=1;                 % duration to calculate tapers over (in s)
    spectro.params.tapers=[spectro.params.W,    spectro.params.T,   0];
    spectro.params.tapers=[3 5]; % [20 4]
    spectro.params.pad=1;
    spectro.params.err=[2 1];
    spectro.params.trialave=0;
    % spectro.params.movingwin=[spectro.params.T spectro.params.T/10];
    % spectro.params.win=1;
    spectro.params.fscorr=0;

spectro.coherence_LFPspike=[];
for channel_id=[L4,L23]
    
    
    inputa=downsample(squeeze(data.filtered_lfp(:,channel_id,:)),20);
    inputb=spikes.PDF.timestamps{channel_id}(:,:)';
    inputa=inputa(200:400,:);
    inputb=inputb(200:400,:);
    
    
    
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
%% choose L4 channel
for sweep_id=data.sweep_sort.successful_sweeps
    spectro.StrongestChannel(sweep_id)=L4;
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
    spectro.params.tapers=[3,5];%[ 3 5]
    spectro.params.pad=1;
    % spectro.params.err=[2 1];
    spectro.params.trialave=0;
    % compute segmented LFP coherence_LFPLFP between channels
    spectro.coherence_LFPLFP=[];
for channel_id=L23
    figure; hold on
    for trial_id = data.sweep_sort.successful_sweeps
        disp(strcat('processing Trial',num2str(trial_id),' Channel',num2str(channel_id)))
        %         inputa = squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,spectro.channeltoanalyse,trial_id));
        %         inputb = squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,channel_id,trial_id));
        
        %         inputa=downsample(squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,spectro.channeltoanalyse,trial_id)),1);
        %         inputb=downsample(squeeze(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2)-1,channel_id,trial_id)),1);
        inputa=downsample(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse(trial_id),trial_id)),20);
        inputb=downsample(squeeze(data.filtered_lfp(:,channel_id,trial_id)),20);
        inputa=zscore(inputa(200:400));
        inputb=zscore(inputb(200:400));
        spectro.coherence_LFPLFP.input1(:,trial_id)=inputa;
        spectro.coherence_LFPLFP.input2(:,trial_id)=inputb;
        plot(inputa,'b')
        plot(inputb,'color',[0.3 0.3 0.3]);
        [spectro.coherence_LFPLFP.C{channel_id}(:,trial_id),...
            spectro.coherence_LFPLFP.phi{channel_id}(:,trial_id),...
            spectro.coherence_LFPLFP.S12{channel_id}(:,trial_id),...
            spectro.coherence_LFPLFP.S1(:,trial_id),...
            spectro.coherence_LFPLFP.S2(:,trial_id),...
            spectro.coherence_LFPLFP.f]=coherencyc(...
            inputa,...
            inputb,...
            spectro.params);
        
        %      spectro.coherence_LFPLFP.S1(:,trial_id)=spectro.coherence_LFPLFP.S1(:,trial_id)./spectro.coherence_LFPLFP.S1(1,trial_id)*-1; % normalize to DC
        spectro.coherence_LFPLFP.S1(:,trial_id)   =10*log10(abs(spectro.coherence_LFPLFP.S1(:,trial_id))); % convert to dB
        
        %      spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id)=spectro.coherence_LFPLFP.S2{channel_id}(:,trial_id)./spectro.coherence_LFPLFP.S2{channel_id}(1,trial_id)   *-1; % normalize to DC
        spectro.coherence_LFPLFP.S2(:,trial_id)   =10*log10(abs(spectro.coherence_LFPLFP.S2(:,trial_id))); % convert to dB
    end
    spectro.coherence_LFPLFP.C{channel_id}  (:,var(spectro.coherence_LFPLFP.C{channel_id})==0)      =   NaN;
    spectro.coherence_LFPLFP.phi{channel_id}(:,var(spectro.coherence_LFPLFP.phi{channel_id})==0)    =   NaN;
    spectro.coherence_LFPLFP.S12{channel_id}(:,var(spectro.coherence_LFPLFP.S12{channel_id})==0)    =   NaN;
    spectro.coherence_LFPLFP.S1             (:,var(spectro.coherence_LFPLFP.S1)==0)                 =   NaN;
    spectro.coherence_LFPLFP.S2             (:,var(spectro.coherence_LFPLFP.S2)==0)                 =   NaN;
    
    
    spectro.coherence_LFPLFP.C_mean(:,channel_id)   = nanmean(spectro.coherence_LFPLFP.C{channel_id},2);
    spectro.coherence_LFPLFP.C_sem(:,channel_id,:)  = nansem(spectro.coherence_LFPLFP.C{channel_id},2);
    
    spectro.coherence_LFPLFP.phi_mean(:,channel_id) = nanmean(abs(spectro.coherence_LFPLFP.phi{channel_id}),2);
    spectro.coherence_LFPLFP.phi_sem(:,channel_id)  = nansem(abs(spectro.coherence_LFPLFP.phi{channel_id}),2);
    
    spectro.coherence_LFPLFP.S12_mean(:,channel_id) = nanmean(abs(spectro.coherence_LFPLFP.S12{channel_id}),2);
    spectro.coherence_LFPLFP.S12_sem(:,channel_id)  = nansem(abs(spectro.coherence_LFPLFP.S12{channel_id}),2);
    
    spectro.coherence_LFPLFP.S1_mean                = nanmean(abs(spectro.coherence_LFPLFP.S1),2);
    spectro.coherence_LFPLFP.S1_sem                 = nansem(abs(spectro.coherence_LFPLFP.S1),2);
    
    spectro.coherence_LFPLFP.S2_mean                = nanmean(abs(spectro.coherence_LFPLFP.S2),2);
    spectro.coherence_LFPLFP.S2_sem                 = nansem(abs(spectro.coherence_LFPLFP.S2),2);
end
% Construct coherency phase histograms

%%
assignin('base', 'spectro', spectro);

