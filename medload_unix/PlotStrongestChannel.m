% function PlotStrongestChannel(data)
StrongestChannel=find(data.max_amp==min(min(data.max_amp)));
% StrongestChannel=12;
x_shift=000; 
y_shift=.2;
temp1  =  staggerplot(squeeze(data.raw_data       (:,StrongestChannel,:)),x_shift,y_shift);
temp2  =  staggerplot(squeeze(data.filtered_lfp   (:,StrongestChannel,:)),x_shift,y_shift);
temp3 =  staggerplot(squeeze(data.filtered_spikes(:,StrongestChannel,:)),x_shift,y_shift);
figure; hold on
% plot raw
subplot(1,3,1); hold on
text(100,2.2,'Unfiltered Data','Fontweight','bold')
plot((1:size(temp1,1))./20,sgolayfilt(temp1,0,15),'k')                    
axis([-20 1020 -0.1 2.3])
xlabel('time (ms)','Fontweight','bold')
ylabel('Amplitude (\muV)','Fontweight','bold')
% plot filtered LFP
subplot(1,3,2); hold on
text(100,2.2,'Filtered LFP','Fontweight','bold')
plot((1:size(temp2,1))./20,sgolayfilt(temp2,0,15),'k');
axis([-20 1020 -0.1 2.3])
xlabel('time (ms)','Fontweight','bold')
ylabel('Amplitude (\muV)','Fontweight','bold')
% plot filtered spikes
subplot(1,3,3), hold on
text(100,2.2,'Filtered spikes','Fontweight','bold')
plot((1:size(temp3,1))./20,sgolayfilt(temp3,0,15),'k');
axis([-20 1020 -0.1 2.3])
xlabel('time (ms)','Fontweight','bold')
ylabel('Amplitude (\muV)','Fontweight','bold')
box off

clear x_shift y_shift temp1 temp2 temp3 
