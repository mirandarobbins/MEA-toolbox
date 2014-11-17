spikes2=spikes;
%%
trim_id=10;%trim_id=trim_id-1
spikes2.clusters.waveforms(spikes2.clusters.assignments==trim_id,:)=[];
spikes2.clusters.spiketimes(spikes2.clusters.assignments==trim_id)=[];
spikes2.clusters.assignments(spikes2.clusters.assignments==trim_id)=[];
spikes2.clusters.cluster_ids(trim_id)=[];
spikes2.clusters.assignments(spikes2.clusters.assignments>trim_id)=spikes2.clusters.assignments(spikes2.clusters.assignments>trim_id)-1
%% plot each cluster separately as a different colour line

    figure;
    for plot_id=1:12%spikes2.clusters.no_clusters-1
%                 subaxis(ceil(sqrt(spikes2.clusters.no_clusters)),ceil(sqrt(spikes2.clusters.no_clusters)),plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.02); 
%         subplot(ceil(sqrt(spikes2.clusters.no_clusters)),ceil(sqrt(spikes2.clusters.no_clusters)),plot_id)
    
        subaxis(4,3,plot_id, 'Spacing', 0)%, 'Padding', 0, 'Margin', 0); 
        hold on
        temp=spikes2.clusters.waveforms(spikes2.clusters.assignments==spikes2.clusters.cluster_ids(plot_id),:)';
%             plot(temp,'color',spikes2.clusters.overcluster.colors(plot_id,:))
       for spike_id=1:size(temp,2)
            temp_x=1:size(temp,1);
            temp_y=temp(:,spike_id);
%             temp_y=sgolayfilt(temp(:,spike_id),5,9);
            patchline(temp_x,temp_y,'edgecolor',spikes2.clusters.overcluster.colors(plot_id,:),'linewidth',1,'edgealpha',0.1);
        end
    plot(mean(temp,2),'color',spikes2.clusters.overcluster.colors(plot_id,:),'LineWidth',2)
    
    axis([-1 spikes2.padding_time(2)+1 -Inf Inf]); box off; axis off
%     text(0,1,num2str(spikes2.clusters.cluster_ids(plot_id)),'Color',[1 1 1])
    end;     
    set(gcf,'color',[0 0 0])
    clear plot_id temp temp_x temp_y