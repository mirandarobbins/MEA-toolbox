function [params, data] = AD_medload_multifile_imaging(parent,first,last,params)
%% Prepare filenames for opening
% data.this_file=fname; clear fname;
file.currentdir= pwd;

if nargin<3
% construct filename
answer=inputdlg({'File parent sequence:','First file in sequence:','last file in sequence:','OR, file numbers:'},...
                'File sequence constructor',1,...
                {'e.g. 20120628','first file number','last file number','eg 1 4 201'});
            
file.parent=answer{1};
file.first =answer{2};
file.last = answer{3};
file.sequence = str2num(answer{4});file.sequence=file.sequence';
else
file.parent=num2str(parent);
file.first=num2str(first);
file.last=num2str(last);

end

% make sequence of filenames... if not already supllied
if isempty(file.sequence)
    file.sequence = str2double(file.first):str2double(file.last);
else
    file.first=file.sequence(1);file.last=file.sequence(numel(file.sequence));
end
for file_id=1:numel(file.sequence)
    this_file=file.sequence(file_id);
    switch numel(num2str(this_file));
        case{1}
            file.this=strcat('00',num2str(this_file));
        case{2}
           file.this=strcat('0',num2str(this_file)); 
        case{3}
           file.this=num2str(this_file); 
    end
    
    file.FileListFull{file_id}=strcat(file.currentdir,'/',file.parent,file.this);
    file.FileListShort{file_id}=strcat(file.parent,file.this);
end
params.files.MEA.dir=file.currentdir;
params.files.MEA.parent=file.parent;
params.files.MEA.FileListShort=file.FileListShort';
params.files.MEA.FileListFull=file.FileListFull';
data.this_file=strcat(params.files.MEA.FileListShort(1),'-',params.files.MEA.FileListShort(numel(params.files.MEA.FileListShort)));
data.this_file=data.this_file{1};
clear file this_file file_id 

if nargin <5
params.Fs=20000;
params.Nyquist=params.Fs/2;
params.selected_rep=1;
params.last_sweep=numel(params.files.MEA.FileListShort);
params.dead_channels=[];
params.frame = 15; %for sgolay
params.degree = 1; 
% params.window_retained=1:40000;% keep 1st 1000ms only
params.start_shift=835;% time of 1st stim in ms (note, will keep and extra 100ms of preceding baseline)

params.window_retained=(params.start_shift-100)*(params.Fs/1000)+1:...
                       (params.start_shift-100)*(params.Fs/1000+params.Fs);
params.baseline_win=([10 99]);
params.monosynaptic_win=([103 115]);

params.search_win=([200, 400]);%for 5x50Hz
% params.search_win=([310, 800]);%for 5x20Hz

params.baseline_win_samples=params.baseline_win*params.Fs/1000+1;
params.monosynaptic_win_samples=params.monosynaptic_win*params.Fs/1000+1;
params.search_win_samples=params.search_win*params.Fs/1000+1;
params.detection_threshold=2;
params.first_stim=100 ;%(ms)
params.channel_index=reshape(1:64,8,8)';

params.flags.rotate=1;
params.flags.denoise=1;
params.flags.prune_failures=1;
params.flags.plot_online=0;
end
%% load and compute mean response
raw=[]; 

for file_id=1:params.last_sweep
    this_file=params.files.MEA.FileListShort{file_id};
    for sweep_id=1:100
        try
            temp{sweep_id}=mload(this_file,sweep_id,500); % open,load to {sweep}(timebase,channel)
%             temp{sweep_id}=temp{sweep_id}(params.window_retained,:);
            catch 
            break
        end
    end
    sweep_id=sweep_id-1;
    data.tb=temp{1}(:,1)-params.start_shift+100; %remove timebase column
    for sweep_id=1:sweep_id
        temp{sweep_id}(:,1)=[];   
    end
    raw=temp;last_sweep=sweep_id; clear temp sweep_id
    params.no_points=size(raw{1},1);
    params.no_channels=size(raw{1},2);
    for sweep_id = 1:last_sweep
        file_temp{file_id}(:,:,sweep_id)=raw{sweep_id}; % wrap 'data.raw' to 3D array(timebase,channel,sweep)
    end
    
end
data.raw_data=file_temp{1};
for file_id=2:params.last_sweep
    data.raw_data=cat(3,data.raw_data,file_temp{file_id}); 
end
params.last_sweep=size(data.raw_data,3);
clear raw file_id file_temp last_sweep sweep_id this_file

% rotate data.raw_data by 90 degrees CW - headstage is now bolted on facing out....
if params.flags.rotate==1
   temp2=cell(8,8);
   for chan_id=1:params.no_channels
       temp(:,chan_id)=reshape(data.raw_data(:,chan_id,:),1,params.no_points*params.last_sweep);
       temp2{chan_id}=temp(:,chan_id);
   end
   clear temp
    temp2=rot90(temp2);
   for chan_id=1:params.no_channels
       temp(:,chan_id)=temp2{chan_id};
       temp3(:,chan_id,:)=reshape(temp(:,chan_id),params.no_points,1,params.last_sweep);
   end
  data.raw_data=temp3;
  clear temp temp2 temp3
end

% denoise 50Hz with Chronux
params.denoise.tapers=[30 5];
params.denoise.fpass=[0 100];
params.denoise.Fs=20000;
if params.flags.denoise==1
for sweep_id=1:params.last_sweep
    data.raw_data(:,:,sweep_id)=rmlinesc(data.raw_data(:,:,sweep_id),params.denoise,0.001/params.denoise.Fs,'n',50);
end;
end
data.raw_data(:,params.dead_channels,:)=NaN;
clear sweep_id


plot(data.filtered_lfp(:,1,1))
%% Automatic detection of failure traces and dead channels
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
%     data.raw=data.raw_data(:,:,data.sweep_sort.successful_sweeps);

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

clear data.raw channel_id sweep_id
%% plot data for all sweeps - one subplot each channel, all trials overlaid
if params.flags.plot_online==1
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
    imagesc(data.max_amp'); figure(gcf) 
    colormap(flipud(hot));
    caxis([-200 0]); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Peak amplitude (mV)')
    title('fEPSP amplitude minima (mV)')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    clear plot_id channel_id sweep_id cs
    % plot latency at minima as heat map
    figure;
    imagesc(data.latency'-params.first_stim); figure(gcf) 
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
        scatter(data.latency(plot_id),(data.max_amp(plot_id))/1000,'or') 
        axis([0 1000 -0.25 0.25])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    text(800,0.15,num2str(plot_id),'FontWeight','bold')
    end; %clear plot_id channel_id sweep_id c
end
% vargout write
assignin('base', 'data', data) 
assignin('base', 'params', params)
end
