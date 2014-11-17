figure('Name',strcat('Activity on strongest channel, each repeat staggered.'),'NumberTitle','off')
    StrongestChannel=find(data.max_amp==min(min(data.max_amp)));
%     StrongestChannel=inputdlg('Which channel to plot?','Choose a channel...');
%     StrongestChannel =str2double(StrongestChannel)
    x_shift=000; 
    y_shift=.2;
    temp1  =  staggerplot(squeeze(data.filtered_lfp       (:,StrongestChannel,data.sweep_sort.successful_sweeps)),x_shift,y_shift);
    temp2  =  staggerplot(squeeze(data.filtered_spikes    (:,StrongestChannel,data.sweep_sort.successful_sweeps)),x_shift,y_shift);
    temp3  =  staggerplot(data.WC.chan1_aligned/300       (:,data.sweep_sort.successful_sweeps)                  ,x_shift,y_shift);
    % plot filtered LFP
    subplot(1,3,1); hold on
       text(100,2.2,'Filtered LFP','Fontweight','bold')
        plot((1:size(temp1,1))./20,sgolayfilt(temp1,0,15),'b');
        axis([-20 850 -0.1 2])
        xlabel('time (ms)')%,'Fontweight','bold')
        ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    % plot filtered spikes
    subplot(1,3,2); hold on
        text(100,2.2,'Filtered spikes','Fontweight','bold')
        plot((1:size(temp2,1))./20,sgolayfilt(temp2,0,15),'b');
        axis([-20 850 -0.1 2])
        xlabel('time (ms)')%,'Fontweight','bold')
    ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    subplot(1,3,3), hold on
        text(100,2.2,'Whole cell','Fontweight','bold')
        plot((1:10000)/10,sgolayfilt(temp3,0,15),'b');
        axis([-20 850 -0.1 2])
        xlabel('time (ms)')%,'Fontweight','bold')
    clear x_shift y_shift temp1 temp2 temp3     