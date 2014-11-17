% remove shitty clusters
all_clusters(isnan(cluster_ids),:)=[];
channel_ids(isnan(cluster_ids),:)=[];
cluster_ids(isnan(cluster_ids))=[];
%%
no_clusters=numel(unique(cluster_ids));
no_channels=numel(unique(channel_ids));
% sorted_clusters=cell(8,8)
sorted_clusters=cell(no_clusters,1);
sorted_channels=cell(no_channels,1);
% sort by clusters
for clust_id=1:no_clusters
    sorted_clusters{clust_id}=...
        all_clusters(cluster_ids==clust_id,:);
end
% sort by channel
for chan_id=1:64
    sorted_channels{chan_id}=...
        all_clusters(channel_ids==chan_id,:);
end
clear chan_id clust_id
%     sorted_clusters=sorted_clusters';
%% plot by cluster
figure; hold on
for clust_id=1:no_clusters
temp=mean(sorted_clusters{clust_id});
    temp=temp-mean(temp(1:10));
    plot (temp)
end 
clear clust_id
%% plot by channel
figure; hold on
chans=unique(channel_ids);
for chan_id=1:numel(chans)
temp=mean(sorted_channels{chans(chan_id)});
    temp=temp-mean(temp(1:10));
plot (temp)
end
clear chan_id chans
%% plot all clusters by channel grid
figure; 
for plot_id =1:64
    subplot(8,8,plot_id); hold on
    plot(sorted_channels{plot_id})
    
end 
clear plot_id
%% plot clusters by channel

figure; 
for plot_id =1:64
    subplot(8,8,plot_id); hold on
    for clust_id=1:no_clusters


        temp(:,clust_id)=mean((all_clusters(channel_ids==plot_id & cluster_ids==clust_id,:)),1);
       
    end
    
%     temp2= mean(temp(1:18,:),1); % baseline correct 
%     for id=1:size(temp,1)
%         temp(id,:)=temp(id,:)-temp2;
%     end
    temp(:,ge(temp(22,:),0))=[]; % remove clusters showing peaks as current sinks
    temp_mean=nanmean(temp,2);
    temp_sem=nansem(temp,2);

        plot(temp_mean,'LineWidth',1.25)
        ciplot((temp_mean-temp_sem),(temp_mean+temp_sem),1:numel(temp_mean),'b');
 
        
%      plot(temp,'LineWidth',1.25)
        axis([0 71 -0.05 0.05])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(3,0.04,num2str(plot_id),'FontWeight','bold')
end





    
    
