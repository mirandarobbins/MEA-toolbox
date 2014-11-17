function AD_replot_multifile(data,params)
%% plot data for all sweeps - one subplot each channel, all trials overlaid

    % new plot for each sweep
%     for sweep_id=data.sweep_sort.successful_sweeps
%         figure; hold on
%         for plot_id=1:params.no_channels
%             subplot(8,8,plot_id); hold on
%             plot(data.tb,sgolayfilt(data.filtered_lfp(:,plot_id,sweep_id),0,15))
%             axis([0 1000 -0.05 0.05])
%             set(gca,'xtick',[])
%             set(gca,'ytick',[])
%             text(3,0.04,num2str(plot_id),'FontWeight','bold')
%         end
%     end; clear sweep_id

%%%%%% LFP all sweeps overlaid
    figure;
    for channel_id=1:params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
%     plot(data.tb,sgolayfilt(squeeze(data.filtered_lfp(:,channel_id,:)),0,15))
        plot(data.tb,squeeze(data.filtered_lfp(:,channel_id,:)))
        axis([0 1000 -0.1 0.1])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(800,0.02,num2str(channel_id),'FontWeight','bold')
    end
%%%%%% MUA all overlaid
    figure;
    for channel_id=1:params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
        plot(data.tb,sgolayfilt(squeeze(data.filtered_spikes(:,channel_id,:)),0,15))
%         plot(data.tb,squeeze(data.filtered_spikes(:,channel_id,:)))
        axis([0 1000 -0.05 0.05])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(800,0.02,num2str(channel_id),'FontWeight','bold')
    end
%%%%%% plot average LFP response
    figure;
    for plot_id=1:params.no_channels
        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
        plot(data.tb,data.mean_channels(:,plot_id),'b','LineWidth',1)
%     ciplot((data.mean_channels(:,plot_id)-data.std_channels(:,plot_id)),(data.mean_channels(:,plot_id)+data.std_channels(:,plot_id)),data.tb,'b');
        axis([0 1000 -0.05 0.05])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(800,0.02,num2str(plot_id),'FontWeight','bold')
    end
%%%%%% plot fEPSP minima as heat map
    figure; 
    imagesc(data.max_amp{params.selected_rep}'); figure(gcf) 
    colormap(flipud(hot));
    caxis([-200 0]); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Peak amplitude (mV)')
    title('fEPSP amplitude minima (mV)')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    clear plot_id channel_id sweep_id cs
%%%%%% plot latency at minima as heat map
    figure;
    imagesc(data.latency{params.selected_rep}'-params.first_stim); figure(gcf) 
    colormap(flipud(jet)); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'time to peak(ms)')
    title('Post-stimulus latency to minima (ms')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 300])

%%%%%% replot to highlight minima
    figure;
    for plot_id=1:params.no_channels

        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
         hold on; box off; axis off
%         plot(data.mean_channels(:,plot_id),'LineWidth',1.25)
    plot(data.tb,sgolayfilt(data.filtered_lfp(:,plot_id,params.selected_rep),params.degree,params.frame),'k')
%     plot(data.tb,data.filtered_lfp(:,plot_id,params.selected_rep),'k')
        plot([params.search_win(1),params.search_win(1)],[-0.05,0.05],'-r')    
        plot([params.search_win(2),params.search_win(2)],[-0.05,0.05],'-r')    
        scatter(data.latency{params.selected_rep}(plot_id),(data.max_amp{params.selected_rep}(plot_id))/1000,'or') 
        axis([0 1000 -0.25 0.25])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    text(800,0.15,num2str(plot_id),'FontWeight','bold')
    end; %clear plot_id channel_id sweep_id c
  %%  
%%%%%% replot LFP for image overlay
    figure;
    for plot_id=1:params.no_channels

        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
         hold on; box off; axis off
%         plot(data.mean_channels(:,plot_id),'LineWidth',1.25)
    plot(data.tb,sgolayfilt(data.filtered_lfp(:,plot_id,params.selected_rep),params.degree,params.frame),'r')
%     plot(data.tb,data.filtered_lfp(:,plot_id,params.selected_rep),'k')
%         plot([params.search_win(1),params.search_win(1)],[-0.05,0.05],'-r')    
%         plot([params.search_win(2),params.search_win(2)],[-0.05,0.05],'-r')    
%         scatter(data.latency(plot_id),(data.max_amp(plot_id))/1000,'or') 
        axis([0 1000 -0.2 0.1])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
%     text(800,0.15,num2str(plot_id),'FontWeight','bold')
    end; %clear plot_id channel_id sweep_id c

%%%%%% plot fEPSP minima as heat map
    figure; 
    imagesc(data.max_amp{params.selected_rep}'); figure(gcf) 
    colormap(hot);
    caxis([-200 0]); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Peak amplitude (mV)')
%     title('fEPSP amplitude minima (mV)')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    clear plot_id channel_id sweep_id cs    
    
% %% side-by side MEA and WC data - VC data
% MEA_closest_channel=input('which MEA channel to plot against?');
% temp_MEA=squeeze(data.raw_data(:, MEA_closest_channel,:));
% temp_MEA=sgolayfilt(temp_MEA,0,15);
% temp_MEA(gt(temp_MEA,0.1)|lt(temp_MEA,-0.2))=NaN;
% temp_MEA=staggerplot(temp_MEA,0, 0.3);
% 
% temp_WC=0.1*sgolayfilt(data.WC.chan1_aligned,0,15);
% temp_WC=staggerplot(temp_WC,0,30);
% figure; hold on
% plot(.001*repmat(data.tb(1:17000),1,size(temp_MEA,2)),temp_MEA(1:17000,:),'r')
% plot(1+repmat(params.tb_WC(1:8500)',1,size(temp_WC,2)),.01*temp_WC,'k')
% axis([ 0 2 0 3.5])
% 
% %% side-by side MEA and WC data - IC data
% MEA_closest_channel=input('which MEA channel to plot against?');
% temp_MEA=squeeze(data.raw_data(:, MEA_closest_channel,:));
% temp_MEA=sgolayfilt(temp_MEA,0,15);
% temp_MEA(gt(temp_MEA,0.1)|lt(temp_MEA,-0.2))=NaN;
% temp_MEA=staggerplot(temp_MEA,0, 0.3);
% 
% temp_WC=sgolayfilt(data.WC.chan1_aligned,0,15);
% temp_WC=staggerplot(temp_WC,0,30);
% figure; hold on
% plot(.001*repmat(data.tb(1:17000),1,size(temp_MEA,2)),temp_MEA(1:17000,:),'r');
% 
% plot(1+repmat(params.tb_WC(1:8500)',1,size(temp_WC,2)),.01*temp_WC,'k')
% axis([ 0 2 0 3.5])
end