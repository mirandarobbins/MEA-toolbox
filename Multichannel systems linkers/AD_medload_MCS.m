function [params, data] = AD_medload_Miranda(data,params)



% if nargin <2 % unless you specify the "params" structure as an argument when calling, these will default 

params.frame = 15; %for savitzy-golay filter kernal
params.degree = 1;  %for savitzy-golay filter kernal

% Frequencies for LFP/MUA separation
params.bp.lowpass_window=[1 200];
params.bp.highpass_window=[1000 5000];

% these parameters are registered in units real time
% to find these times, run: plot(data.tb,data.raw_data(:,2,1))

params.window_retained  =1:1000;       % keep 1st 50ms only
params.baseline_win     =([1 9]);       % search time window for baseline calculation
params.monosynaptic_win =([15 40]); % search time window for monosynaptic stim calculation
params.search_win       =([15, 40]);% search time window for hunting for polysynaptic activity

%convert the above into sample 
params.baseline_win_samples=params.baseline_win*params.Fs/1000;
params.monosynaptic_win_samples=params.monosynaptic_win*params.Fs/1000;
params.search_win_samples=params.search_win*params.Fs/1000;
params.detection_threshold=3; % noise threshold for multi/single unit finding, e.g. 3x rms basline
params.spontaneous_threshold=0.15; % positive-going threshold (in mV) to tag as contaminated with spontaneous activity)
params.first_stim=10 ;%(ms)
params.channel_index=reshape(1:64,8,8)';
% params.channel_index([1 8 57 64])=NaN;

%Some internal flags
params.flags.rotate=0;          % needed to rotate the array 90 degrees, if the headstage is bolted sideways on the scope
params.flags.denoise=0;         % runs 50Hz noise removal (experimental)
params.flags.prune_failures=1;  % hunts for and removes response failures (i.e. no fEPSP activity in seaerch_win
params.flags.prune_spontaneous=1;
params.flags.plot_online=1;
params.flags.mask_channels=1;
% end
%% Choose stim channels
if params.flags.plot_online==1;
   figure;
   for channel_id=1:64
       for trial_id=1:10
            subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
            plot(data.tb,data.raw_data(:,channel_id,trial_id))
            
            axis([0 max(data.tb) -0.5 0.5]) %axes of the individual graphs
%             axis([0 1000 -inf inf]) %axes of the individual graphs

            set(gca,'xtick',[])
            set(gca,'ytick',[])
       end
       
       if channel_id==1 | channel_id==8 | channel_id==57 |channel_id==64
            set(gca, 'color', [0.6 0.6 0.6])
       else
            text(max(data.tb)*0.8,0.0002,num2str(channel_id),'FontWeight','bold','FontSize',8)
       end
   end
end   

if params.flags.mask_channels==1;
    temp=inputdlg('Choose dead channels or stimulation channels to blank ( e.g. 1 35:37 64)','channe blanking',1,{'33'});
    if ~isempty(temp)
        params.dead_channels= str2num(temp{1});
        data.raw_data(:,params.dead_channels,:)=0;
    end
end; clear temp

%% denoise 50Hz with Chronux
params.denoise.tapers=[30 5]; %[3 5];
params.denoise.fpass=[0 100];
params.denoise.Fs=20000;
if params.flags.denoise==1
for sweep_id=1:params.last_sweep
    data.raw_data(:,:,sweep_id)=rmlinesc(data.raw_data(:,:,sweep_id),params.denoise,0.001/params.denoise.Fs,'n',50);
end;
end
data.raw_data(:,params.dead_channels,:)=NaN;
clear sweep_id raw
%
%% Automatic detection of failure traces
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
%     if params.flags.plot_online==1
        figure; hold on; plot(data.sweep_sort.window_std_max); 
            plot ([1 params.last_sweep],[params.detection_threshold params.detection_threshold],':k')
            axis([1 params.last_sweep 0 Inf])
            xlabel('Sweep no.')
            ylabel('Normalised window std.dev.')
%     end
    temp=1:params.last_sweep;
    data.sweep_sort.successful_sweeps=temp(gt(data.sweep_sort.window_std_max,params.detection_threshold));
    data.sweep_sort.failures=temp(lt(data.sweep_sort.window_std_max,params.detection_threshold));
%     params.last_sweep=numel(data.sweep_sort.successful_sweeps);
%     data.trimmed=data.raw_data(:,:,data.sweep_sort.successful_sweeps);    
else %otherwise, include all
    data.sweep_sort.successful_sweeps=1:params.last_sweep;
    data.sweep_sort.successful_sweeps=[];
%     data.raw_data=data.raw_data(:,:,data.sweep_sort.successful_sweeps);
end

% Automatic detection of trace swith spontaneous activity
if params.flags.prune_spontaneous==1;
    temp_all=[];
    for trial_id=1:size(data.raw_data,3)
        temp_all(:,trial_id)=max(data.raw_data(:,:,trial_id));
    end
    
    temp=1:params.last_sweep;
    temp(gt(nanmean(temp_all),params.spontaneous_threshold));
    data.sweep_sort.spontaneous_sweeps = temp(gt(nanmean(temp_all),params.spontaneous_threshold));
    data.sweep_sort.successful_sweeps(ismember(data.sweep_sort.successful_sweeps,data.sweep_sort.spontaneous_sweeps))=[]

end

clear temp_win temp_baseline temp    
%% bandpass filter to separate LFP and spikes
params.bp.wn_low=params.bp.lowpass_window/params.Nyquist;     %Freqs for low pass filter (Hz/Nyquist) ->LFP
params.bp.wn_high=params.bp.highpass_window/params.Nyquist; %Freqs for low pass filter (Hz/Nyquist) ->MUA
    
[params.bp.B_low,  params.bp.A_low]   = butter(1,params.bp.wn_low,'bandpass');
[params.bp.B_high, params.bp.A_high]  = butter(5,params.bp.wn_high,'bandpass');

%     fvtool(params.bp.B_low,params.bp.A_low,'Fs',params.Fs)  
%     fvtool(params.bp.B_high,params.bp.A_high,'Fs',params.Fs)    

for trial_id=1:params.last_sweep
    for channel_id=1:64
        clear temp; 
        temp=data.raw_data(:,channel_id,trial_id);
        data.filtered_lfp(:,channel_id,trial_id)      = filter(params.bp.B_low,params.bp.A_low,temp);
        data.filtered_spikes(:,channel_id,trial_id)   = filter(params.bp.B_high,params.bp.A_high,temp);
    end
end
          
          
clear temp trial_id


%% calculate mean,std LFP responses

for channel_id=1:params.no_channels
    data.mean_channels(:,channel_id)=squeeze(nanmean(data.filtered_lfp(:,channel_id,data.sweep_sort.successful_sweeps),3));
    data.std_channels(:,channel_id) =squeeze( nanstd(data.filtered_lfp(:,channel_id,data.sweep_sort.successful_sweeps),3));
%     
%     data.mean_channels(:,channel_id)=sgolayfilt(mean(data.filtered_lfp(:,channel_id,:),3),params.degree,params.frame);
%     data.std_channels(:,channel_id) =sgolayfilt(nanstd(data.filtered_lfp(:,channel_id,:),3),params.degree,params.frame);
end

% fEPSP minima, latency
data.max_amp_mono=zeros(8,8);
data.latency_mono=zeros(8,8);
data.max_amp=zeros(8,8);
data.latency=zeros(8,8);
for channel_id=1:params.no_channels
%     [data.max_amp(channel_id) data.latency(channel_id)]=min(data.mean_channels(params.search_win(1):params.search_win(2),channel_id)); % plot mean
       [data.max_amp(channel_id) data.latency(channel_id)]=min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),channel_id,params.selected_rep)); % plot selected rep
    if gt(mean(data.mean_channels(params.monosynaptic_win_samples(1):params.monosynaptic_win_samples(2),channel_id)),...
          mean(data.mean_channels(params.baseline_win_samples(1):params.baseline_win_samples(2),channel_id))); % plot selected rep
       [data.max_amp_mono(channel_id) data.latency_mono(channel_id)]=max(data.mean_channels(params.monosynaptic_win_samples(1):params.monosynaptic_win_samples(2),channel_id)); % plot selected rep
    else
       [data.max_amp_mono(channel_id) data.latency_mono(channel_id)]=min(data.mean_channels(params.monosynaptic_win_samples(1):params.monosynaptic_win_samples(2),channel_id)); % plot selected rep
    end
end
data.max_amp=data.max_amp*1000; % convert amplitude from mV to uV
data.latency=((data.latency)+params.search_win_samples(1))/params.Fs*1000;
data.max_amp(params.dead_channels)=10; data.latency(params.dead_channels)=NaN;%hide dead channels

data.max_amp_mono=data.max_amp_mono*1000; % convert amplitude from mV to uV
data.latency_mono=((data.latency_mono))/params.Fs*1000;
data.max_amp_mono(params.dead_channels)=10; data.latency_mono(params.dead_channels)=NaN;%hide dead channels

data.burst_timing.latency=cell(params.last_sweep,1);
data.burst_timing.amp=cell(params.last_sweep,1);
for trial_id=1:params.last_sweep
% for trial_id=data.sweep_sort.successful_sweeps % <-- careful here!
data.burst_timing.latency{trial_id}=zeros(8,8); 
data.burst_timing.amp{trial_id}=zeros(8,8); 
for channel_id=1:params.no_channels
    [data.burst_timing.amp{trial_id}(channel_id) data.burst_timing.latency{trial_id}(channel_id)]=...
        min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),channel_id,trial_id)); % plot selected rep
end
data.burst_timing.amp{trial_id}=data.burst_timing.amp{trial_id}*1000; % convert amplitude from mV to uV
data.burst_timing.latency{trial_id}=(data.burst_timing.latency{trial_id}+params.search_win_samples(1))/params.Fs*1000;

end
clear data.raw_data channel_id sweep_id trial_id
%% plot data for all sweeps - one subplot each channel, all trials overlaid
if params.flags.plot_online==1
%%%%%% LFP all sweeps overlaid
    figure;
    for channel_id=1:params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
%     plot(data.tb,sgolayfilt(squeeze(data.filtered_lfp(:,channel_id,:)),0,15))
        plot(data.tb,squeeze(data.filtered_lfp(:,channel_id,:)))
% plot(...
%      sgolayfilt(staggerplot(squeeze(data.filtered_lfp(:,channel_id,:)),1000,0.05),...
%      0,15));  
        axis([0 max(data.tb) -0.5 0.5])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(max(data.tb)*0.8,0.02,num2str(channel_id),'FontWeight','bold','FontSize',8)
    end
%%%%%% MUA all overlaid
    figure;
    for channel_id=1:params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
        plot(data.tb,sgolayfilt(squeeze(data.filtered_spikes(:,channel_id,:)),0,15))
%         plot(data.tb,squeeze(data.filtered_spikes(:,channel_id,:)))
        axis([0 max(data.tb) -0.01 0.01])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(max(data.tb)*0.8,0.02,num2str(channel_id),'FontWeight','bold','FontSize',8)
    end
%%%%%% plot average LFP response
    figure;
    for plot_id=1:params.no_channels
        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
        plot(data.tb,data.mean_channels(:,plot_id),'b','LineWidth',1)
%     ciplot((data.mean_channels(:,plot_id)-data.std_channels(:,plot_id)),(data.mean_channels(:,plot_id)+data.std_channels(:,plot_id)),data.tb,'b');
        axis([0 max(data.tb) -0.5 0.5])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(max(data.tb)*0.8,0.02,num2str(channel_id),'FontWeight','bold','FontSize',8)
    end
%%%%%% plot fEPSP minima as heat map
    figure; 
    imagesc(data.max_amp'); figure(gcf) 
    colormap(flipud(hot));
%     caxis([-100 0]);
    c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Peak amplitude (mV)')
    title('fEPSP amplitude minima (uV)')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    clear plot_id channel_id sweep_id cs
%%%%%% plot latency at fEPSP minima as heat map
    figure;
    imagesc(data.latency'-params.first_stim); figure(gcf) 
    colormap(flipud(summer)); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'time to peak(ms)')
    title('Post-stimulus latency to minima (ms')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 
        max(max(data.latency'-params.first_stim))])

%%%%%% replot to highlight minima
    figure;
    for plot_id=1:params.no_channels

        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
         hold on; box off; axis off
%         plot(data.mean_channels(:,plot_id),'LineWidth',1.25)
    plot(data.tb,sgolayfilt(data.filtered_lfp(:,plot_id,params.selected_rep),params.degree,params.frame),'b')
%     plot(data.tb,data.filtered_lfp(:,plot_id,params.selected_rep),'k')
        plot([params.search_win(1),params.search_win(1)],[-0.05,0.05],'-r')    
        plot([params.search_win(2),params.search_win(2)],[-0.05,0.05],'-r')    
        scatter(data.latency(plot_id),(data.max_amp(plot_id))/1000,'or') 
        axis([0 max(data.tb) -0.5 0.5])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(max(data.tb)*0.8,0.02,num2str(plot_id),'FontWeight','bold','FontSize',8)
    end; %clear plot_id channel_id sweep_id c
end
% vargout write
assignin('base', 'data', data) ;
assignin('base', 'params', params);
end
