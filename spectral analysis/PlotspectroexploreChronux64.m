%% plot spectrograms    
for trial_id=data.sweep_sort.successful_sweeps%1:10
    spectro.channeltoanalyse=mode(spectro.StrongestChannel);%spectro.StrongestChannel(trial_id);
    
figure;set(gcf,'position',[20 185 400 800])
subplot(311); hold on
    plot(data.tb,sgolayfilt(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,trial_id)),0,15),'-k')%'Color',[0.8 0.8 0.8]); 
    plot(data.tb,mean(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,data.sweep_sort.successful_sweeps)),2),':b')
%     plot(data.tb,mean(squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,data.sweep_sort.failures)),2),':r')
    text(690,0.18,strcat('Activity on strongest channel (Ch',num2str(spectro.channeltoanalyse),')'))
    axis([min(min(spectro.T))*1000 max(max(spectro.T))*1000 -0.3 0.3])
    xlabel('time (ms)')    
    ylabel('LFP (mV)')
    box off
ax1=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);

subplot(312); hold on
%     imagesc(repmat(spectro.T,numel(spectro.params.freqs),1)',...
%             repmat(spectro.params.freqs,numel(spectro.T),1),...
%             spectro.mean_norm');


x= repmat(spectro.T,size(spectro.mean{spectro.channeltoanalyse},2),1)*1000;
y=repmat(spectro.F,size(spectro.mean{spectro.channeltoanalyse},1),1)';
    z=spectro.P{spectro.channeltoanalyse}(:,:,trial_id)';
% z=spectro.mean';
z=smooth2a(z,1,1);
ih=surf(x,y,z); view(2); axis tight;
set(ih, 'edgecolor', 'none');
set(ih, 'facecolor', 'interp'); 

        

    colormap((jet)); 
    caxis([-70 -50]) ; 
%  set(gca,'YScale','log')
    
ax2a=gca;
pos=get(ax2a,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);
pos=get(gca,'pos');
    hc=colorbar('location','eastoutside',...
                'position',[pos(1)+pos(3)+0.01 pos(2) 0.03 pos(4)*0.5],...
                'YAxisLocation','right');
ylabel(hc,'dBuV/Hz')
title('Power Spectral Density Estimate')
ylabel('Frequency (Hz)');
xlabel('Time (ms)');  

ax2a = gca;
ax2b = axes('Position',get(ax2a,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k',...
           'Xtick',[],...
           'YTick',[],...
           'xlim',[0, 1],...
           'ylim',[0, 100]);
       
       %search win
rectangle('Position',[params.search_win(1)/1000,...
                      min(spectro.F),...
                      (params.search_win(2)-params.search_win(1))/1000,...
                      max(spectro.F)-min(spectro.F)],...
                      'LineWidth',2,'LineStyle','-','EdgeColor',[0 0 0],'Parent',ax2b)
                  
       % minima
rectangle('Position',[spectro.burst_centre{spectro.channeltoanalyse}(trial_id),...
                      min(spectro.F),...
                      0.001,...
                      max(spectro.F)-min(spectro.F)],...
                      'LineWidth',2,'LineStyle','-','EdgeColor',[1 0 0],'Parent',ax2b)
       %found win
% rectangle('Position',[spectro.burst_centre{spectro.channeltoanalyse}(trial_id)-spectro.power_ROI(1)/1000,...
%                       min(spectro.F),...
%                       (spectro.burst_centre{spectro.channeltoanalyse}(trial_id)+spectro.power_ROI(2)/1000)-(spectro.burst_centre{spectro.channeltoanalyse}(trial_id)-spectro.power_ROI(1)/1000),...
%                       max(spectro.F)-min(spectro.F)],...
%                       'LineWidth',2,'LineStyle','-','EdgeColor',[0 1 0],'Parent',ax2b)                  
                  
subplot(313)
pos=get(gca,'pos');
   set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);   
channel_to_plot=spectro.channeltoanalyse;
hold on
plot(spectro.F,spectro.average_window_power{spectro.channeltoanalyse}(:,trial_id),'-r','LineWidth',2)
% plot(spectro.F,spectro.average_window_power{spectro.channeltoanalyse},'r')
% plot averages
plot(spectro.F,spectro.average_window_power_mean{spectro.channeltoanalyse},'-k','LineWidth',2);
    ciplot(spectro.average_window_power_mean{spectro.channeltoanalyse}+spectro.average_window_power_sem{spectro.channeltoanalyse},...
       spectro.average_window_power_mean{spectro.channeltoanalyse}-spectro.average_window_power_sem{spectro.channeltoanalyse},...
       spectro.F,'k')

% set(gca,'XScale','log')
% set(gca,'YScale','log')
xlabel('Frequency (Hz)')
axis([1 100 -80 -30])
title('Power spectral density')
ylabel('dBuV/Hz')

end
clear ax1 ax2a ax2b hc ih pos trial_id x y z 
%% plot high power bands
climits=[-50 0];
for trial_id=data.sweep_sort.successful_sweeps
figure('Name',strcat('burst-locked spectral power: trial...', num2str(trial_id)),'NumberTitle','off')
    set(gcf,'pos',[143 707 748 196])
    subaxis(1,4,1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
    temp=spectro.window_power.twenty{trial_id}';
    temp=smooth2a(temp,1,1);
     ih=surf(temp); view(2); axis tight;
    set(ih, 'edgecolor', 'none');
    set(ih, 'facecolor', 'interp'); 
    axis square; 
%     colormap(hot); 
%     caxis(climits);
%     colorbar 
    set(gca,'xtick',[]);set(gca,'ytick',[]); grid off
    text(4.5,7.5,'Power @ 20Hz','Color','g','FontSize',12)
    
    subaxis(1,4,2, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
    temp=spectro.window_power.fifty{trial_id}';
        temp=smooth2a(temp,1,1);
     ih=surf(temp); view(2); axis tight;
    set(ih, 'edgecolor', 'none');
    set(ih, 'facecolor', 'interp'); 
    axis square; 
%     colormap(hot); 
%     caxis(climits);
%     colorbar 
    set(gca,'xtick',[]);set(gca,'ytick',[]); grid off
    text(4.5,7.5,'Power @ 50Hz','Color','g','FontSize',12)

    subaxis(1,4,3, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
    temp=spectro.window_power.onehunderd{trial_id}';
        temp=smooth2a(temp,1,1);
     ih=surf(temp); view(2); axis tight;
    set(ih, 'edgecolor', 'none');
    set(ih, 'facecolor', 'interp'); 
    axis square; 
%     colormap(hot); 
%     caxis(climits);
%     colorbar 
    set(gca,'xtick',[]);set(gca,'ytick',[]); grid off
    text(4.25,7.5,'Power @ 100Hz','Color','g','FontSize',12)
    
    subaxis(1,4,4, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
    temp=spectro.window_power.twohunderd{trial_id}';
        temp=smooth2a(temp,1,1);
     ih=surf(temp); view(2); axis tight;
    set(ih, 'edgecolor', 'none');
    set(ih, 'facecolor', 'interp'); 
    axis square; 
%     colormap(hot); 
%     caxis(climits);
%     colorbar 
    set(gca,'xtick',[]);set(gca,'ytick',[]); grid off
    text(4.25,7.5,'Power @ 200Hz','Color','g','FontSize',12)

end
