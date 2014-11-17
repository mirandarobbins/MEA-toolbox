function [spikes] = AD_spike_extract(data,params)
spikes=[];
% spikes = MEA_spikes;
% spikes.coeffthreshold=10;
spikes.time_conversion=200000; % Sampling freq
spikes.f=spikes.time_conversion/2; %FFT up to nyquist frequency
spikes.detection_thresh=8; %nb or if using rms, (line 33) change to 7.5
spikes.detection_rearm_time=1; % re-arm time for peakfinsing (20 points= 1ms @20kHz, good for clean unit sorting)

%%%% BLANKING TIMES %%%%
% reject spikes in baseline period  (first 100ms + 5ms to allow stim to return to baseline)
% spikes.blanking_times=1:20; %to ignore
% spikes.blanking_times=1:2099; % for single stim @ 100ms
%%%% blanking times for 5x20Hz stim

% spikes.blanking_times=1:6100; %for 5x20Hz
% spikes.blanking_times=[1:2000,1980:2100,... 
%                               2980:3100,... 
%                               3980:4100,...
%                               4980:5100,...
%                               5980:6100]; %for 5x20Hz to hide stim artifacts

%%%% blanking times for 5x50Hz stim
spikes.blanking_times=1:4100; %for 5x50Hz
% spikes.blanking_times=[1:2000,1980:2060,... 
%                               2380:2460,... 
%                               2780:2860,...
%                               3180:3260,...
%                               3580:3660]; %for 5x20Hz to hide stim artifacts

spikes.padding_time=[10, 20];
% totallength=0;
% spikes_sparse=[];



%flags
if nargin<3
flags.neg_trigger=1;
flags.plot_unfilteredfiltered_comparison=0;
flags.plot_filteredspike_comparison=0;
flags.flag_plotconvolutions=0;
flags.pass2Chronux=1;
flags.plot_aligned_waveforms=0;
flags.plot_cluster_figures=0;
end
% normalise spikes by std, find peak that cross detection thresh using high-BP filtered data
for channel_id=1:size(data.filtered_spikes,2)
%     for trial_id=1:size(data.filtered_spikes,3)
    for trial_id=data.sweep_sort.successful_sweeps
        clear temp
        temp=data.filtered_spikes(:,channel_id,trial_id)./...
            rms(data.filtered_spikes(params.baseline_win_samples(1):params.baseline_win_samples(2),channel_id,trial_id));
        data.filtered_spikes(:,channel_id,trial_id)=temp;
        temp_sweep=data.filtered_spikes(:,channel_id,trial_id);
        clear peak
      
        if flags.neg_trigger==1;
            % For negative-going triggers 
                [peak(:,2) peak(:,1)]= findpeaks(-1*temp_sweep,...
                                         'minpeakheight',spikes.detection_thresh,...
                                         'minpeakdistance',spikes.detection_rearm_time);
                peak(:,2)=peak(:,2)*-1;
        else
            % For positive-going triggers
                [peak(:,2) peak(:,1)]= findpeaks(temp_sweep,...
                                         'minpeakheight',spikes.detection_thresh,...
                                         'minpeakdistance',spikes.detection_rearm_time);
                peak(:,2)=peak(:,2);
        end
        peak(:,1)=peak(:,1);%/spikes.time_conversion*1000 % keep in sample time...
%         peak(lt(peak(:,1),max(spikes.blanking_times)),:)=[];
%         peak(lt(peak(:,1),max(spikes.blanking_times)),:)=[];

    for id=1:numel(spikes.blanking_times)
        
        if(find(spikes.blanking_times(id)==peak(:,1)));
        peak(peak(:,1)==spikes.blanking_times(id),:)=[];
        end
        
    end


spikes.spike_locs{channel_id,trial_id}=peak;
        clear temp_sweep
    end
end
% clear channel_id trial_id peak temp

%%

if flags.plot_filteredspike_comparison==1
% plor all reps overlaid, each channel separate
figure; 
for channel_id=1:size(data.filtered_spikes,2)
    subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
    box off; axis off
%     for trial_id=1:size(data.filtered_spikes,3)
    for trial_id=data.sweep_sort.successful_sweeps
        hold on
%     plot(data.tb,sgolayfilt(data.filtered_spikes(:,channel_id,trial_id),0,15),'k')
        plot(data.tb, data.filtered_spikes(:,channel_id,trial_id),'k')
        plot([min(data.tb) max(data.tb)],[spikes.detection_thresh spikes.detection_thresh],':r')
        plot([min(data.tb) max(data.tb)],[-spikes.detection_thresh -spikes.detection_thresh],':r')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        scatter(spikes.spike_locs{channel_id,trial_id}(:,1)/spikes.time_conversion*10000,spikes.spike_locs{channel_id,trial_id}(:,2),'og')
        axis([0 max(data.tb) -10 10]) % nb yaxis is now SD 
        
    end
end
% plor all reps overlaid
figure; 
for channel_id=1:size(data.filtered_spikes,2)
    for trial_id=data.sweep_sort.successful_sweeps
        hold on
        plot(data.tb, data.filtered_spikes(:,channel_id,trial_id),'k')
        plot([min(data.tb) max(data.tb)],[spikes.detection_thresh spikes.detection_thresh],':r')
        plot([min(data.tb) max(data.tb)],[-spikes.detection_thresh -spikes.detection_thresh],':r')
        scatter(spikes.spike_locs{channel_id,trial_id}(:,1)/spikes.time_conversion*10000,spikes.spike_locs{channel_id,trial_id}(:,2)/500,'og')
    end
end
end; clear channel_id trial_id
%% cut out spike waveforms (i.e. for unit sorting) and spike-triggered averate LFP
spikes.clusters=[];
spikes.waveforms_array=cell(size(data.filtered_spikes,2),size(data.filtered_spikes,3));
spikes.spiketimes=cell(size(data.filtered_spikes,2),size(data.filtered_spikes,3));
spikes.clusters.waveforms=[];
idx=1;
            if flags.plot_aligned_waveforms==1
                figure;
            end
for channel_id=1:size(data.filtered_spikes,2)
%     for trial_id=1:size(data.filtered_spikes,3)
    for trial_id=data.sweep_sort.successful_sweeps

        for spike_id=1:size(spikes.spike_locs{channel_id,trial_id,1})
            this_spike=spikes.spike_locs{channel_id,trial_id}(spike_id,1);
            spikes.spiketimes{channel_id,trial_id}(spike_id)=(spikes.spike_locs{channel_id,trial_id}(spike_id,1)/spikes.time_conversion*1000);
            temp=data.filtered_spikes(this_spike-spikes.padding_time(1):this_spike+spikes.padding_time(2),channel_id,trial_id);
            

            spikes.waveforms_array{channel_id,trial_id}(spike_id,:)=temp;
            spikes.clusters.waveforms(idx,:)=temp;
            spikes.clusters.spiketimes(idx,1)=this_spike./10000; %nb chronux uses seconds not ms
            if flags.plot_aligned_waveforms==1
                plot(temp)
                hold on
            else
            end
            idx=idx+1;
        end
    end

end; clear idx temp this_spike spike_id trial_id channel_id
%% pass to Chronux
if flags.pass2Chronux==1
% Prepare for Chronux
    spikes.clusters.Fs= spikes.time_conversion;
    spikes.clusters.threshV=[-1*spikes.detection_thresh,0];
    spikes.clusters.threshT=spikes.padding_time(1)+1;
%process with Chronux
    spikes.clusters = ss_dejitter(spikes.clusters,'com');
    spikes.clusters = ss_outliers(spikes.clusters);
    spikes.clusters = ss_kmeans(spikes.clusters);
    spikes.clusters = ss_energy(spikes.clusters);
    spikes.clusters = ss_aggregate(spikes.clusters);
    spikes.clusters.assignments = spikes.clusters.hierarchy.assigns;
% plot results
if flags.plot_cluster_figures==1
%     ssg_databrowse3d(spikes.clusters); grid on
%     ssg_databrowse2d(spikes.clusters); grid on
    figure; colormap jet;
        subplot(2,1,1); plot(spikes.clusters.waveforms'); axis tight; title('Centered Data w/ Outliers Removed');
        subplot(2,1,2); histxt(spikes.clusters.waveforms);
    figure;  set(gcf, 'Renderer', 'OpenGL');
        clusterXT(spikes.clusters, spikes.clusters.overcluster.assigns);  title('Local Clusters');
    figure; set(gcf, 'Renderer', 'OpenGL');
        clusterXT(spikes.clusters, spikes.clusters.hierarchy.assigns); title('Final Clusters');
%     figure; colormap jet;
%         showclust(spikes.clusters, spikes.clusters.hierarchy.assigns);
%     figure; 
%         aggtree(spikes.clusters); title('Aggregation Tree');
%     figure;
%         correlations(spikes.clusters);  title('Auto- and Cross- Correlations');
    figure; 
        scatter(spikes.clusters.spiketimes,spikes.clusters.assignments,'ob')
        axis([0 1 0 Inf])
% plot each cluster seperated
else
end
spikes.clusters.cluster_ids=unique(spikes.clusters.overcluster.assigns);
spikes.clusters.no_clusters=numel(spikes.clusters.cluster_ids);
if flags.plot_cluster_figures==1
    %% plot each cluster seperately as a different colour line
    figure;
    for plot_id=1:spikes.clusters.no_clusters-1
%                 subaxis(ceil(sqrt(spikes.clusters.no_clusters)),ceil(sqrt(spikes.clusters.no_clusters)),plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.02); 
%         subplot(ceil(sqrt(spikes.clusters.no_clusters)),ceil(sqrt(spikes.clusters.no_clusters)),plot_id)
    subaxis(5,5,plot_id, 'Spacing', 0)%, 'Padding', 0, 'Margin', 0); 
        hold on
        temp=spikes.clusters.waveforms(spikes.clusters.assignments==spikes.clusters.cluster_ids(plot_id),:)';
%             plot(temp,'color',spikes.clusters.overcluster.colors(plot_id,:))
       for spike_id=1:size(temp,2)
            temp_x=1:size(temp,1);
            temp_y=temp(:,spike_id);
%             temp_y=sgolayfilt(temp(:,spike_id),5,9);
            patchline(temp_x,temp_y,'edgecolor',spikes.clusters.overcluster.colors(plot_id,:),'linewidth',1,'edgealpha',0.1);
        end
    plot(mean(temp,2),'color',spikes.clusters.overcluster.colors(plot_id,:),'LineWidth',2)
    
    axis([-1 spikes.padding_time(2)+1 -Inf Inf]); box off; axis off

    end;     
    set(gcf,'color',[0 0 0])
    clear plot_id temp temp_x temp_y
end
end
assignin('base', 'spikes', spikes);

 %if Chronux...
%%
%  figure
%  for plot_id=1:spikes.clusters.no_clusters
%      hold on
%     spikes.clusters.histo(:,plot_id)=histc(spikes.clusters.spiketimes(spikes.clusters.assignments==spikes.clusters.cluster_ids(plot_id-1),:),1:200:20000);
%     stairs(1:10:1000,spikes.clusters.histo(:,plot_id),'color',spikes.clusters.overcluster.colors(plot_id,:))
%     stairs(1:10:1000,sum(spikes.clusters.histo,2),'LineWidth',1.5) 
%  end; clear plot_id
