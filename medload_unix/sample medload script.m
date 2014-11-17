%% search params
params.Fs=20000;
params.Nyquist=params.Fs/2;
params.selected_rep=1;
params.last_sweep=10;
params.dead_channels=[];
params.frame = 15; %for sgolay
params.degree = 1; 
params.window_retained=1:20000;% keep 1st 1000ms only
params.baseline_win=([10 99]);
params.search_win=([200, 800]);
params.baseline_win_samples=params.baseline_win*params.Fs/1000+1;
params.search_win_samples=params.search_win*params.Fs/1000+1;
params.detection_threshold=5;
params.first_stim=100 ;%(ms)

params.flags.denoise=1;
params.flags.prune_failures=1;

% data.this_file='20120511016';
% data.this_file='20120426012';
data.this_file='20120707003';
%% open 1 channel,filter, plot
% 
% figure % load all + plot
% data.raw=mload(data.this_file,params.selected_rep,20);
% data.tb=data.raw(:,1); data.raw(:,1)=[];
% data.raw(:,params.dead_channels)=NaN;
% params.no_channels = size(data.raw,2);
% 
% for plot_id=1:params.no_channels
%     subplot(8,8,plot_id); hold on
%     plot(data.tb,data.raw(:,plot_id)) % unfiltered
%     plot(data.tb,sgolayfilt(data.raw(:,plot_id),params.degree,params.frame)) %filtered
%     plot([params.search_win(1),params.search_win(1)],[-0.05,0.05],'-r')    
%     plot([params.search_win(2),params.search_win(2)],[-0.05,0.05],'-r')    
%      axis([0 1000 -0.05 0.05])
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         text(800,0.03,num2str(plot_id),'FontWeight','bold')
% end ;clear plot_id
%% load and compute mean response
data.raw=[];

for sweep_id=1:params.last_sweep
    temp=mload(data.this_file,sweep_id,100); % open,load to {sweep}(timebase,channel)
    temp=temp(params.window_retained,:);
    data.raw{sweep_id}=temp;
end
clear temp

data.tb=data.raw{1}(:,1); 
for sweep_id=1:params.last_sweep
    data.raw{sweep_id}(:,1)=[];%remove timebase column
end

% some parameters
params.no_points=size(data.raw{1},1);
params.no_channels=size(data.raw{1},2);
for sweep_id = 1:params.last_sweep
    data.raw_data(:,:,sweep_id)=data.raw{sweep_id}; % wrap 'data.raw' to 3D array(timebase,channel,sweep)
end



% denoise 50Hz with Chronux
params.denoise.tapers=[3 5];
params.denoise.fpass=[0 100];
if params.flags.denoise==1
for sweep_id=1:params.last_sweep
    data.raw_data(:,:,sweep_id)=rmlinesc(data.raw_data(:,:,sweep_id),params.denoise,0.001/params.Fs,'n',50);
end;
end
data.raw_data(:,params.dead_channels,:)=NaN;
data.raw=[];  clear sweep_id

%% Automatic removal of failure traces and dead channels
if params.flags.prune_failures==1;
    
%   temp_win      =  sgolayfilt(data.denoised(params.search_win_samples(1):params.search_win_samples(2),:,:),0,15);
%   temp_baseline =  sgolayfilt(data.raw_data(params.baseline_win_samples(1):params.baseline_win_samples(2),:,:),0,15);

    temp_win      =  data.raw_data(params.search_win_samples(1):params.search_win_samples(2),:,:);
    temp_baseline =  data.raw_data(params.baseline_win_samples(1):params.baseline_win_samples(2),:,:);

%   for plot_id=1:64
%       subplot(8,8,plot_id)
%       plot(squeeze(temp_win(:,plot_id,:)))
%   end

    for channel_id=1:params.no_channels
        for sweep_id=1:params.last_sweep
            data.sweep_sort.window_std(channel_id,sweep_id) =nanstd(temp_win(:,channel_id,sweep_id))./...
                                                             nanstd(temp_baseline(:,channel_id,sweep_id));
        end
    end
    data.sweep_sort.window_std_max=nanmax(data.sweep_sort.window_std);

    figure; hold on; plot(data.sweep_sort.window_std_max); 
        plot ([1 params.last_sweep],[params.detection_threshold params.detection_threshold],':k')
        axis([1 params.last_sweep 0 Inf])
        xlabel('Sweep no.')
        ylabel('Normalised window std.dev.')
        temp=1:params.last_sweep;
        
    data.sweep_sort.successful_sweeps=temp(gt(data.sweep_sort.window_std_max,params.detection_threshold));
    params.last_sweep=numel(data.sweep_sort.successful_sweeps);
    data.trimmed=data.raw_data(:,:,data.sweep_sort.successful_sweeps);
    
else %otherwise, include all
    data.sweep_sort.successful_sweeps=1:params.last_sweep;
    data.trimmed=data.raw_data(:,:,data.sweep_sort.successful_sweeps);

end
clear temp_win temp_baseline temp
%% bandpass filter to separate LFP and spikes

params.bp.wn_low=[1/params.Nyquist 200/params.Nyquist];
params.bp.wn_high=[100/params.Nyquist 5000/params.Nyquist];
    
[params.bp.B_low,  params.bp.A_low]   = butter(1,params.bp.wn_low,'bandpass');
[params.bp.B_high, params.bp.A_high]  = butter(5,params.bp.wn_high,'bandpass');

%     fvtool(params.bp.B_low,params.bp.A_low,'Fs',params.Fs)  
%     fvtool(params.bp.B_high,params.bp.A_high,'Fs',params.Fs)    

for trial_id=1:params.last_sweep
    for channel_id=1:64
        clear temp; 
        temp=data.trimmed(:,channel_id,trial_id);
        data.filtered_lfp(:,channel_id,trial_id)      = filter(params.bp.B_low,params.bp.A_low,temp);
        data.filtered_spikes(:,channel_id,trial_id)   = filter(params.bp.B_high,params.bp.A_high,temp);
    end
end
          
          
clear temp trial_id


%% calculate mean,std LFP responses

for channel_id=1:params.no_channels
    data.mean_channels(:,channel_id)=squeeze(nanmean(data.filtered_lfp(:,channel_id,:),3));
    data.std_channels(:,channel_id) =squeeze( nanstd(data.filtered_lfp(:,channel_id,:),3));
%     
%     data.mean_channels(:,channel_id)=sgolayfilt(mean(data.filtered_lfp(:,channel_id,:),3),params.degree,params.frame);
%     data.std_channels(:,channel_id) =sgolayfilt(nanstd(data.filtered_lfp(:,channel_id,:),3),params.degree,params.frame);
end
clear data.raw channel_id sweep_id
%% plot data for all sweeps - one subplot each channel, all trials overlaid
% % clf; 
% for sweep_id=data.sweep_sort.successful_sweeps
%     figure; hold on
%     for plot_id=1:params.no_channels
%         subplot(8,8,plot_id); hold on
%         plot(data.tb,sgolayfilt(data.filtered_lfp(:,plot_id,sweep_id),0,15))
%         axis([0 1000 -0.05 0.05])
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         text(3,0.04,num2str(plot_id),'FontWeight','bold')
%     end
% end
clear sweep_id

figure;% LFP all overlaid
for channel_id=1:params.no_channels
    subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
%  plot(data.tb,sgolayfilt(squeeze(data.filtered_lfp(:,channel_id,:)),0,15))
    plot(data.tb,squeeze(data.filtered_lfp(:,channel_id,:)))
    axis([0 1000 -0.1 0.1])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    text(800,0.02,num2str(channel_id),'FontWeight','bold')
end
clear sweep_id
figure;% spikes all overlaid
for channel_id=1:params.no_channels
    subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
    plot(data.tb,sgolayfilt(squeeze(data.filtered_spikes(:,channel_id,:)),0,15))
%  plot(data.tb,squeeze(data.filtered_spikes(:,channel_id,:)))
    axis([0 1000 -0.05 0.05])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    text(800,0.02,num2str(channel_id),'FontWeight','bold')
end

%% plot average response
figure;
for plot_id=1:params.no_channels
    subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
    plot(data.tb,data.mean_channels(:,plot_id),'b','LineWidth',1)
%     ciplot((data.mean_channels(:,plot_id)-data.std_channels(:,plot_id)),(data.mean_channels(:,plot_id)+data.std_channels(:,plot_id)),data.tb,'b');
    axis([0 1000 -0.05 0.05])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    text(800,0.02,num2str(channel_id),'FontWeight','bold')
end
%% fEPSP minima, latency 

data.max_amp=zeros(8,8);
data.latency=zeros(8,8);
for channel_id=1:params.no_channels
%     [data.max_amp(channel_id)
%     data.latency(channel_id)]=min(data.mean_channels(params.search_win(1):params.search_win(2),channel_id)); % plot mean
    [data.max_amp(channel_id) data.latency(channel_id)]=min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),channel_id,params.selected_rep)); % plot selected rep
end
data.max_amp=data.max_amp*1000; % convert amplitude from mV to uV
data.latency=((data.latency)+params.search_win_samples(1))/params.Fs*1000;
data.max_amp(params.dead_channels)=10; data.latency(params.dead_channels)=NaN;%hide dead channels

%% plot spread of activity 
% plot fEPSP minima as heat map
figure; 

imagesc(data.max_amp'); figure(gcf) 
colormap(flipud(hot));
caxis([-200 0]); c= colorbar%('title','fEPSP minima')
ylabel(c,'Peak amplitude (mV)')
title('fEPSP amplitude minima (mV)')
box off

        set(gca,'xtick',[])
        set(gca,'ytick',[])

clear plot_id channel_id sweep_id cs

% plot latency at minima as heat map
figure;
imagesc(data.latency'-params.first_stim); figure(gcf) 
colormap(flipud(jet)); c= colorbar%('title','fEPSP minima')
ylabel(c,'time to peak(ms)')
title('Post-stimulus latency to minima (ms')
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

clear plot_id channel_id sweep_id c
caxis([0 300])
%% replot to highlight minima
% plot(data.tb,data.trimmed(:,32,params.selected_rep))
figure;
for plot_id=1:params.no_channels
%     subplot(8,8,plot_id); hold on
        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
        box off; axis off
%     plot(data.mean_channels(:,plot_id),'LineWidth',1.25)
% plot(data.tb,sgolayfilt(data.filtered_lfp(:,plot_id,params.selected_rep),params.degree,params.frame),'k')
plot(data.tb,data.filtered_lfp(:,plot_id,params.selected_rep),'k')
    plot([params.search_win(1),params.search_win(1)],[-0.05,0.05],'-r')    
    plot([params.search_win(2),params.search_win(2)],[-0.05,0.05],'-r')    
scatter(data.latency(plot_id),(data.max_amp(plot_id))/1000,'or') 
axis([0 1000 -0.25 0.25])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
text(800,0.15,num2str(plot_id),'FontWeight','bold')
end; clear plot_id 