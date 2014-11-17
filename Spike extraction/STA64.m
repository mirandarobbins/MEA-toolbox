%% Spike-triggered average LFP
spikes.STA=[];
%%%%%%% convert spike times to seconds, pad array
temp_spike_times2=cell(1,64);
temp_spike_times=cell(64,params.last_sweep);
    for channel_id=1:64
        temp_spike_times2{channel_id}=zeros(1000,params.last_sweep);
    end
for channel_id=1:64
    for trial_id = 1:params.last_sweep
        temp_spike_times{channel_id,trial_id}=spikes.spiketimes{channel_id,trial_id}./100; % convert to seconds
        temp=temp_spike_times{channel_id,trial_id};
%         temp=temp*20000; % convert back to samples (@ 20KHz)
        temp_spike_times2{channel_id}(1:numel(temp),trial_id)=temp;
        temp_spike_times2{channel_id}(temp_spike_times2{channel_id}==0)=NaN;
%         temp_spike_times2{channel_id}=round(temp_spike_times2{channel_id}*1000);
        if ~isempty(temp_spike_times2{channel_id})
        end
%         a=temp_spike_times2{channel_id}(:,trial_id);
%         spikes.STA.input_spikes_array(channel_id)=struct(strcat('Ch_',num2str(channel_id)),[])
%         struct('times',rand(1,100));
    end
end
spikes.STA.spike_times=temp_spike_times2;
% data(1)=struct('times',rand(1,100));
clear temp_spike_times temp_spike_times2 temp trial_id channel_id
%% cut out spike waveforms (i.e. for unit sorting) and spike-triggered average LFP
flags.plot_aligned_waveforms=0;
spikes.STA.padding_time=[50 50];
spikes.STA.padding_time_midpoint=ceil(numel(-1*spikes.STA.padding_time(1):spikes.STA.padding_time(2))/2);

spikes.STA.waveform_array=cell(size(data.filtered_spikes,2),numel(params.last_sweep));
idx=1;
for channel_id=1:64
    for trial_id=1:params.last_sweep
        a=spikes.STA.spike_times{channel_id}(:,trial_id); a(isnan(a))=[];
        if numel(a)>0
            no_spikes=numel(a);
        for spike_id=1:no_spikes % can choose 1st spike only here...
            this_spike=round(a(spike_id)*1000);            
            temp_data=downsample(squeeze(data.filtered_lfp(:,channel_id,trial_id)),20);
            
            temp_LFP=temp_data(this_spike-spikes.STA.padding_time(1):this_spike+spikes.STA.padding_time(2));
            temp_LFP=temp_LFP-temp_LFP(spikes.STA.padding_time_midpoint);
            spikes.STA.waveform_array{channel_id,trial_id}(spike_id,:)=temp_LFP;

%             spikes.waveforms_array{channel_id,trial_id}(spike_id,:)=temp;
%             spikes.clusters.waveforms(idx,:)=temp;
%             spikes.clusters.spiketimes(idx,1)=this_spike./10000; %nb chronux uses seconds not ms
            if flags.plot_aligned_waveforms==1
                plot(temp_LFP)
                hold on
            else
            end
            idx=idx+1;
        end
        end
    end
end

clear idx temp this_spike spike_id trial_id channel_id temp_LFP

%% Average STA per channel
spikes.STA.waveforms_mean=zeros(64,(spikes.STA.padding_time(1)+spikes.STA.padding_time(2)+1));spikes.STA.waveforms_mean(spikes.STA.waveforms_mean==0)=NaN;
spikes.STA.waveforms_sem=zeros(64,(spikes.STA.padding_time(1)+spikes.STA.padding_time(2)+1));spikes.STA.waveforms_sem(spikes.STA.waveforms_sem==0)=NaN;
spikes.STA.waveforms=cell(64,1);
for channel_id=1:64
    temp=spikes.STA.waveform_array{channel_id,1};  

    for trial_id=2:numel(data.sweep_sort.successful_sweeps)
        temp2=spikes.STA.waveform_array{channel_id,trial_id};  
        temp=vertcat(temp,temp2);
    end
    spikes.STA.waveforms{channel_id}=temp;
    if ~isempty(temp)
    spikes.STA.waveforms_mean(channel_id,:)=nanmean(spikes.STA.waveforms{channel_id});
%     spikes.STA.waveforms_mean(channel_id,:)=spikes.STA.waveforms_mean(channel_id,:)-spikes.STA.waveforms_mean(channel_id,spikes.STA.padding_time_midpoint);
%     spikes.STA.waveforms_mean(channel_id,:)=spikes.STA.waveforms_mean(channel_id,:)-mean(spikes.STA.waveforms_mean(channel_id,:));
    spikes.STA.waveforms_sem(channel_id,:) =nansem(spikes.STA.waveforms{channel_id});
    end
end
clear temp temp2 trial_id channel_id
%% plot as subplot
tb=spikes.STA.padding_time(1)*-1:spikes.STA.padding_time(2);
figure; 
for plot_id=1:64
subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
hold on
% cind=[rand rand rand];
cind='b';
plot(tb,spikes.STA.waveforms_mean(plot_id,:)*1000,'color',cind','LineWidth',2)
ciplot(spikes.STA.waveforms_mean(plot_id,:)*1000-spikes.STA.waveforms_sem(plot_id,:)*1000,...
       spikes.STA.waveforms_mean(plot_id,:)*1000+spikes.STA.waveforms_sem(plot_id,:)*1000,...
       tb,cind);
   axis([-1*spikes.STA.padding_time(1) spikes.STA.padding_time(2) -Inf Inf])
   plot([0 0],[-50 50],':k')

end
%% plot overlaid
tb=spikes.STA.padding_time(1)*-1:spikes.STA.padding_time(2);
figure; hold on
for plot_id=1:64
cind=[rand rand rand];
plot(tb,spikes.STA.waveforms_mean(plot_id,:)*1000,'color',cind','LineWidth',2)
% ciplot(spikes.STA.waveforms_mean(plot_id,:)*1000-spikes.STA.waveforms_sem(plot_id,:)*1000,...
%        spikes.STA.waveforms_mean(plot_id,:)*1000+spikes.STA.waveforms_sem(plot_id,:)*1000,...
%        tb,cind)
   plot([0 0],[-50 50],':k')
   axis([-1*spikes.STA.padding_time(1) spikes.STA.padding_time(2) -40 40])

end



