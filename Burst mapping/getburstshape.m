function [burst] = getburstshape(data,params,CSD)
% Uses CSD to find centre of burst
% Thresholds shoulder of burst region by STDEV to estimate number of coactive channels (active radius)
% Reports number and variability of channels participating in burst.


%%
burst=[];
tic
burst.params.no_points=params.search_win_samples(2)-params.search_win_samples(1)+1;
burst.params.detectionthreshold = 0.6; % fraction of array-wide peak amplitude to call channels active at
burst.StrongestChannel.map      =   zeros(8,8,burst.params.no_points,numel(data.sweep_sort.successful_sweeps));
% burst.StrongestChannel.coords   =   zeros(params.no_points,2,numel(data.sweep_sort.successful_sweeps));
burst.StrongestChannel.coords   =   zeros(burst.params.no_points,2,numel(data.sweep_sort.successful_sweeps));

burst.ActiveChannels.map        =   zeros(8,8,burst.params.no_points,numel(data.sweep_sort.successful_sweeps));
burst.ActiveChannels.grand_map  =   zeros(8,8,numel(data.sweep_sort.successful_sweeps));
burst.params.useCSD=isvar('base','CSD'); %else, use LFP data
burst.params.plotYN=1; % online plotting 

burst.Location.useModePosition=0;
burst.burst_window_tb=((params.search_win_samples(1):params.search_win_samples(2))-params.search_win_samples(1))*params.Fs^-1*1000; % in ms
% burst.burst_window_tb=params.search_win(1)-1:params.search_win(2);
if ~burst.params.useCSD; disp('CSD data not found... analysing LFP data'); else disp('using CSD data'); end
for trial_id=data.sweep_sort.successful_sweeps
    disp(strcat('Mapping burst coordinates for trial...',num2str(trial_id)))
%     for tb_id=1:params.no_points
    for tb_id=1:burst.params.no_points
        tb_temp=params.search_win_samples(1):params.search_win_samples(2);
        if burst.params.useCSD
            temp=CSD.csd_array(tb_temp(tb_id),:,trial_id);
        else
            temp=data.filtered_lfp(tb_temp(tb_id),:,trial_id);
        end
        % map position of strongest channel
        burst.StrongestChannel.amp(tb_id,trial_id)=   min(temp);
        burst.StrongestChannel.id(tb_id,trial_id)=   find(temp==min(temp));
        map_temp=zeros(8,8);map_temp(burst.StrongestChannel.id(tb_id,trial_id))=1;
        burst.StrongestChannel.map(:,:,tb_id,trial_id)=map_temp;
        [col row]=find(map_temp,1);        
        burst.StrongestChannel.coords(tb_id,:,trial_id)=[col,row];
        % map no of active channels
        if burst.params.useCSD
           burst.ActiveChannels.map_1D(tb_id,:,trial_id)=temp<burst.params.detectionthreshold*min(min(CSD.csd_array(params.search_win_samples(1):params.search_win_samples(2),:,trial_id))); %min for this sample
        else
           %         burst.ActiveChannels.map_1D(tb_id,:,trial_id)=temp<burst.params.detectionthreshold*min(temp); %min for this sample
           burst.ActiveChannels.map_1D(tb_id,:,trial_id)=temp<burst.params.detectionthreshold*min(min(data.burst_timing.amp{trial_id})/1000); %min for this sample
        end
        clear temp map_temp row col
    end
    % histogram of raw strongest position
    burst.StrongestChannel.hist{trial_id}=nansum(squeeze(burst.StrongestChannel.map(:,:,:,trial_id)),3)';
    burst.StrongestChannel.hist{trial_id}=burst.StrongestChannel.hist{trial_id}./sum(sum(burst.StrongestChannel.hist{trial_id}));
end
    toc

burst.ActiveChannels.no_ActiveChannels=squeeze(sum(burst.ActiveChannels.map_1D,2));
burst.ActiveChannels.no_ActiveChannels_mean=mean(burst.ActiveChannels.no_ActiveChannels,2);
burst.ActiveChannels.no_ActiveChannels_sem=nansem(burst.ActiveChannels.no_ActiveChannels,2);

clear tb_id trial_id tb_temp

% sanity check... plot centre of burst by timestep
figure;hold on
for trial_id=data.sweep_sort.successful_sweeps
plot(burst.StrongestChannel.amp(:,trial_id),burst.StrongestChannel.id(:,trial_id));
end
figure;hold on
subplot(2,1,1);plot(burst.StrongestChannel.amp);
subplot(2,1,2);plot(burst.StrongestChannel.id);
%% recenter burst position to most active channel, get position/speed stats
for trial_id=data.sweep_sort.successful_sweeps
    
%%%% extract search win for analysis
%     burst.Location.peak_position{trial_id}=burst.StrongestChannel.coords(params.search_win_samples(1):params.search_win_samples(2),:,trial_id);
    burst.Location.peak_position{trial_id}=burst.StrongestChannel.coords(:,:,trial_id);
if burst.Location.useModePosition;
   %%%% recenter position such that most active channel has location [0,0]
        % x
        burst.Location.peak_displacement{trial_id}(:,1)=burst.Location.peak_position{trial_id}(:,1)-mode(burst.Location.peak_position{trial_id}(:,1));
        
        burst.Location.peak_displacement{trial_id}(:,2)=burst.Location.peak_position{trial_id}(:,2)-mode(burst.Location.peak_position{trial_id}(:,2));
else
   %%%%OR... recenter by first position:
        % x
        burst.Location.peak_displacement{trial_id}(:,1)=burst.Location.peak_position{trial_id}(:,1)-burst.Location.peak_position{trial_id}(1,1);
        % y
        burst.Location.peak_displacement{trial_id}(:,2)=burst.Location.peak_position{trial_id}(:,2)-burst.Location.peak_position{trial_id}(1,2);
end
%%%%position vector
for tb_id=1:burst.params.no_points%(params.search_win_samples(2)-params.search_win_samples(1)+1)
    burst.Location.peak_displacement_vector(trial_id,tb_id)=sqrt((burst.Location.peak_displacement{trial_id}(tb_id,1))^2 +...
                                                                 (burst.Location.peak_displacement{trial_id}(tb_id,2))^2);
end

burst.Location.peak_displacement_vector(isinf(burst.Location.peak_displacement_vector))=0;
    
%%%% histogram of time spent at each position (for x and y coords)
    burst.Location.displacement_hist{trial_id}(:,1)=...
        hist(burst.Location.peak_displacement{trial_id}(:,1),8)./numel(burst.Location.peak_displacement{trial_id}(:,1));
    burst.Location.displacement_hist{trial_id}(:,2)=...
        hist(burst.Location.peak_displacement{trial_id}(:,2),8)./numel(burst.Location.peak_displacement{trial_id}(:,2));
    
%%%% histogram of time spent at each position (for 8x8 coords)
    burst.Location.map_hist(:,:,trial_id)=...
        sum(burst.StrongestChannel.map(:,:,:,trial_id),3)./size(burst.StrongestChannel.map(:,:,:,trial_id),3);
    
%%%% Convert from Cartesian -> polar coords     
%     burst.peak_displacement_hist_theta{trial_id}=atan2(burst.peak_displacement_hist{trial_id}(:,2),burst.peak_displacement_hist{trial_id}(:,1));
%     burst.peak_displacement_hist_rho{trial_id}=sqrt(burst.peak_displacement_hist{trial_id}(:,1).^2+burst.peak_displacement_hist{trial_id}(:,2).^2);

end
%%%% Grand map position and coordinates
burst.Location.grand_map_hist=sum(burst.Location.map_hist,3)./size(burst.Location.map_hist,3);

[burst.Location.grand_map_maximum(1),burst.Location.grand_map_maximum(2)]=...
    find(burst.Location.grand_map_hist==max(max(burst.Location.grand_map_hist)));

burst.Location.grand_map_xaxis=repmat((1:8),8,1)';
burst.Location.grand_map_yaxis=repmat((1:8),8,1);
burst.Location.grand_map_xyzaxis=horzcat(reshape(burst.Location.grand_map_xaxis,64,1),...
                                        reshape(burst.Location.grand_map_yaxis,64,1),...
                                        reshape((burst.Location.grand_map_hist)',64,1));

burst.Location.grand_map_xyzaxis(:,1)=burst.Location.grand_map_xyzaxis(:,1)-burst.Location.grand_map_maximum(1);
burst.Location.grand_map_xyzaxis(:,2)=burst.Location.grand_map_xyzaxis(:,2)-burst.Location.grand_map_maximum(2);

%%%% mean displacement vector
burst.Location.peak_displacement_vector_mean=nanmean(burst.Location.peak_displacement_vector);
burst.Location.peak_displacement_vector_sem=nansem(burst.Location.peak_displacement_vector);
        
%%%% speed of burst (1st diff of displacement vector)
burst.Location.peak_velocity=diff(burst.Location.peak_displacement_vector,1,2);
for idx=1:numel(burst.Location.peak_velocity)
    if burst.Location.peak_velocity(idx)<0
        burst.Location.peak_velocity(idx)=burst.Location.peak_velocity(idx)*-1;
    else
    end
end % shameful....

burst.Location.peak_velocity=burst.Location.peak_velocity*3;% covert to m/s... this is the same as *150e-6(/(params.Fs^-1);
burst.Location.peak_velocity_mean=nanmean(burst.Location.peak_velocity);
burst.Location.peak_velocity_sem=nansem(burst.Location.peak_velocity);

%%%% histogram of burst speed
burst.Location.peak_velocity_bins=1:10;
temp=reshape(burst.Location.peak_velocity,numel(burst.Location.peak_velocity),1);
burst.Location.peak_velocity_histo=(histc(temp,burst.Location.peak_velocity_bins))./numel(temp);
%% plot 
if burst.params.plotYN==1
figure; 
subplot(4,1,1); hold on 
% plot(burst.burst_window_tb,burst.Location.peak_displacement_vector,'-r','LineWidth',1);
plot(burst.burst_window_tb,burst.Location.peak_displacement_vector_mean*150,'-k','LineWidth',2);
        ciplot(burst.Location.peak_displacement_vector_mean*150+burst.Location.peak_displacement_vector_sem*150,...
            burst.Location.peak_displacement_vector_mean*150-burst.Location.peak_displacement_vector_sem*150,...
            burst.burst_window_tb,'k');
xlabel('Time (ms)'); ylabel('(um)');title('Mean distance from initiation site')
axis([min(burst.Location.peak_velocity_bins) max(burst.burst_window_tb) 0 Inf])


subplot(4,1,2); hold on 
% plot(burst.burst_window_tb(1:numel(burst.burst_window_tb)-1),burst.Location.peak_velocity,'-r','LineWidth',1);  
plot(burst.burst_window_tb(1:numel(burst.burst_window_tb)-1),burst.Location.peak_velocity_mean,'-k','LineWidth',1);
        ciplot(burst.Location.peak_velocity_mean+burst.Location.peak_velocity_sem,...
            burst.Location.peak_velocity_mean-burst.Location.peak_velocity_sem,...
            burst.burst_window_tb(1:numel(burst.burst_window_tb)-1),'k');
xlabel('Time (ms)'); title('Mean burst velocity');ylabel('(m/s)')
axis([min(burst.burst_window_tb) max(burst.burst_window_tb) 0 Inf])

subplot(4,1,3); hold on 
bar(burst.Location.peak_velocity_bins,burst.Location.peak_velocity_histo,'k')
axis([min(burst.Location.peak_velocity_bins) max(burst.Location.peak_velocity_bins) 0 Inf])
xlabel('Burst speed (m/s)'); ylabel('Prob.');title('Burst velocity - probability distribution')
%        set(gca,'YScale','log')

subplot(4,1,4); hold on 
% plot(burst.burst_window_tb(1:numel(burst.burst_window_tb)),burst.ActiveChannels.no_ActiveChannels,'-r','LineWidth',1);  
plot(burst.burst_window_tb(1:numel(burst.burst_window_tb)),burst.ActiveChannels.no_ActiveChannels_mean,'-k','LineWidth',1);
        ciplot(burst.ActiveChannels.no_ActiveChannels_mean+burst.ActiveChannels.no_ActiveChannels_sem,...
            burst.ActiveChannels.no_ActiveChannels_mean-burst.ActiveChannels.no_ActiveChannels_sem,...
            burst.burst_window_tb(1:numel(burst.burst_window_tb)),'k');
        xlabel('Time (ms)'); ylabel('no.'); title('no. coactive channels (m/s)')
        axis([min(burst.burst_window_tb) max(burst.burst_window_tb) 0 Inf])
end
%% number and position of active channels
% for trial_id=data.sweep_sort.successful_sweeps
%         temp=CSD.csd_array(params.search_win_samples(1):params.search_win_samples(2),:,trial_id);
%         temp2(gt(temp,burst.params.detectionthreshold*min(min(temp)))=1
%         burst.ActiveChannels.grand_map(gt(temp,
%     
%% Plot average activity shifted to 0,0 by max activity (modal time)
x=reshape(burst.Location.grand_map_xyzaxis(:,1),8,8);
y=reshape(burst.Location.grand_map_xyzaxis(:,2),8,8);
z=reshape(burst.Location.grand_map_xyzaxis(:,3),8,8);
% z=zeros(8,8);
if burst.params.plotYN==1
figure; hold on
title('Burst position during window: probability of focal location (relative channel location)')
    sh=pcolor(y,x,z);
    set(sh,'edgecolor','none')
    plot([-8 8],[0 0],':r'); plot([0 0],[-8 8],':r'); grid off
    colormap(flipud(hot))
    colorbar;
    caxis([0 0.2])
    axis([-8 +8 -8 +8 ]); view(90,90)
    clear x y z
end
%% Convert Cartesian to Polar and consruct displacement histogram

[burst.Location.grand_map_polaraxis(:,1) ...
 burst.Location.grand_map_polaraxis(:,2) ...
 burst.Location.grand_map_polaraxis(:,3)]   =...
                       cart2pol(burst.Location.grand_map_xyzaxis(:,1),...
                                burst.Location.grand_map_xyzaxis(:,2),...
                                burst.Location.grand_map_xyzaxis(:,3));
if burst.params.plotYN==1
    figure; title('Stability of network activity: displacement of burst position from most frequent channel')
    wind_rose2((-180/pi)*burst.Location.grand_map_polaraxis(:,1),burst.Location.grand_map_polaraxis(:,3),...
               'n',8, 'cmap',hot,'di',[0 1],'lcolor','b');
end
%% Plot burst maxima stability, each trial overlaid
if burst.params.plotYN==1
    figure; title('Stability of network activity: displacement of burst position from most frequent channel')
for trial_id=1:numel(data.sweep_sort.successful_sweeps)

        subaxis(floor(sqrt(numel(data.sweep_sort.successful_sweeps))),ceil(sqrt(numel(data.sweep_sort.successful_sweeps))),trial_id, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0.01);  
        posdata=burst.Location.peak_displacement{data.sweep_sort.successful_sweeps(trial_id)};
hold on
        text(2,5,strcat('Trial no.',num2str(trial_id)))
        plot(posdata(:,2),posdata(:,1));%PlotAxisAtOrigin
        axis([-8 8 -8 8])
%             view(90,-90)
         grid on
end
end
%% Plot burst's fractional dwell time on each channel, one subplot per trial
if burst.params.plotYN==1
    figure;
    for trial_id=1:numel(data.sweep_sort.successful_sweeps)
            subaxis(floor(sqrt(numel(data.sweep_sort.successful_sweeps))),ceil(sqrt(numel(data.sweep_sort.successful_sweeps))),trial_id, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0.01);              title(strcat('Trial no.',num2str(trial_id)))
            plot_data=(burst.Location.map_hist(:,:,data.sweep_sort.successful_sweeps(trial_id)))';
        imagesc(plot_data)
        colormap(flipud(gray)); caxis([0 0.2]); axis square; grid off; axis on; box on
        set(gca,'XTick',[],'YTick',[])
    end

%     % Plot burst moving around in RT
%     figure; hold on
%     for trial_id=2%data.sweep_sort.successful_sweeps
%        figure; hold on
%        set(gcf,'DoubleBuffer','on')
%        posdata=burst.StrongestChannel.coords(params.search_win_samples(1):params.search_win_samples(2),:,trial_id);
%         imgdata=burst.StrongestChannel.map(:,:,params.search_win_samples(1):params.search_win_samples(2),trial_id);
% 
%         for tb_id=1:10:size(posdata,1)
%             cla 
%     %         imgdata=interp2(burst.strongest_channel_map(:,:,tb_id,trial_id), 2, 'linear');
%             imagesc((imgdata(:,:,tb_id)))
%             alpha 0.3
%             colormap(flipud(gray))
%             plot(posdata(1:tb_id,1),posdata(1:tb_id,2)*-1,'LineWidth',1.5)
%     %         axis([0 8 0 8]);
%             drawnow update 
%              grid on
%     % % view(-90,90) % 
%             th=text(7,1,strcat(num2str(tb_id./20),'ms'));
%             set(th, 'fontsize', 14, 'color', 'k')
%         pause(0.001)
%         end
%     end
end
clear trial_id
assignin('base', 'burst', burst);