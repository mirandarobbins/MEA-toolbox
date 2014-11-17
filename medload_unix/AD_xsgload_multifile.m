function [params, data]=AD_xsgload_multifile(data,params)
 
% rmfield(data,'WC')
%% Prepare filename for opening
% construct filename
file.currentdir= pwd;
file.parent=('AA0001AAAA');
if nargin<3
params.flags.plot_online=0;
end
if nargin<4
% file.directory=dir(pwd)
% filename='AA0001AAAA0001.xsg';
answer=inputdlg({'First file in sequence:','Last file in sequence:','Channel to import:','Nearest MEA channel:'},...
                'File sequence constructor',1,...
                {'first file number','last file number','1 or 2','(e.g. 16)'});
     
file.first          = num2str(answer{1});
file.last           = num2str(answer{2});
file.channel        = num2str(answer{3});
data.WC.closestMEA  = str2num(answer{4});
else
file.first          = num2str(first);
file.last           = num2str(last);
file.channel        = num2str(channel);
end
switch lower (file.channel)
    case{'1'}
%         file.channel='data.ephys.trace_1';
        disp('analysing channel 1...')
    case{'2'}
%         file.channel='data.ephys.trace_2';
        disp('analysing channel 2...')
    otherwise
        disp('Unknown channel requested!!')
end

% make sequence of filenames
file.sequence = str2double(file.first):str2double(file.last);
clear answer
%% load loop
params.files.WC.data_time_range=1:10000;
for file_id=1:numel(file.sequence)
    this_file=file.sequence(file_id);
    switch numel(num2str(this_file));
        case{1}
            file.this=strcat('000',num2str(this_file));
        case{2}
            file.this=strcat('00',num2str(this_file));
        case{3}
            file.this=strcat('0',num2str(this_file));  
    end
    
    file.currentfile=strcat(file.currentdir,'/',file.parent,file.this,'.xsg');
    temp=load(file.currentfile, '-mat');
    
    params.files.WC.dir=file.currentdir;
    params.files.WC.parent=file.parent;
    params.files.WC.FileListFull{file_id,1}=strcat(file.currentfile);
    params.files.WC.FileListShort{file_id,1}=strcat(file.parent,file.this);
    
    switch (file.channel)
        case{'1'}
            output(:,file_id)=temp.data.ephys.trace_1(1:10000);
            disp(strcat('importing channel 1: ',file.parent,file.this,'.xsg'))
            data.WC.chan_raw=output;

        case{'2'}
            output(:,file_id)=temp.data.ephys.trace_2(1:10000);
            disp(strcat('importing channel 2: ',file.parent,file.this,'.xsg'))
            data.WC.chan_raw=output;
        otherwise
            disp('Unknown channel requested!!')
    end

end
params.fs_WC=10000;
params.tb_WC=(0:size(output,1)-1)/params.fs_WC;
% example plot
if params.flags.plot_online==1
   figure; hold on
   for plot_id=1:size(output,2)
       plot(params.tb_WC,output(:,plot_id),'b')
   end
   title('Imported data')
   xlabel('time (s)')
   ylabel('Vm (mV)')
   grid on
end
clear as output  file_id  file.this this_file

%% clean up whole cell data
% is data VC or IC?
if max(max(diff(data.WC.chan_raw)))>1000 % look for v fast resolution at stim artifacts
    params.flags.WCisVoltageClamp=1;
    thresh=5;
    disp('data is from a Voltage clamp recording')
else
    params.flags.WCisVoltageClamp=0;
    disp('data is from a Current clamp recording')
    thresh=3;
end

data.WC.stim_time=[];
data.WC.chan_aligned=[];
for sweep_id=1:size(data.WC.chan_raw,2)
data.WC.chan_aligned(:,sweep_id)=data.WC.chan_raw(:,sweep_id)-(mean(data.WC.chan_raw(10:990,sweep_id)));

data.WC.chan_dVdt(:,sweep_id)=diff(data.WC.chan_aligned(:,sweep_id)); data.WC.chan_dVdt(:,sweep_id)=data.WC.chan_dVdt(:,sweep_id)./std(data.WC.chan_dVdt(:,sweep_id));
[pks(sweep_id),data.WC.stim_time(sweep_id)]=findpeaks(data.WC.chan_dVdt(:,sweep_id),'MINPEAKHEIGHT',thresh,'NPEAKS',1);%'THRESHOLD',1,
% data.WC.stim_time(sweep_id)=data.WC.stim_time(sweep_id)+2000;
end

data.WC.TTL_jitter=std(data.WC.stim_time)/10; %timing jitter in ms
% % % % plot to check stim artifacts have been correctly tagged
% figure; hold on
%    title('Imported data')
%    for plot_id=1:size(data.WC.chan_dVdt,2)
%        subplot(ceil(sqrt(size(data.WC.chan_aligned,2))),...
%        ceil(sqrt(size(data.WC.chan_aligned,2))),plot_id); hold on
%        plot(data.WC.chan_dVdt(:,plot_id),'b')
%        scatter(data.WC.stim_time(plot_id),pks(sweep_id),'xr')
% %        plot(data.WC.chan_aligned(:,plot_id),'b')
% %        scatter(data.WC.stim_time(plot_id),0,'xr')
%        xlabel('time (s)')
%        ylabel('Vm (mV)')
%        grid on
%    end
   
for sweep_id=1:size(data.WC.chan_raw,2)
    data.WC.time_shift(sweep_id)=(data.WC.stim_time(sweep_id)-1000); % trim so 1st stim is at 100ms
clear temp
temp=data.WC.chan_aligned(:,sweep_id); temp(1:data.WC.time_shift(sweep_id))=[];
temp=vertcat(temp,(zeros(data.WC.time_shift(sweep_id),1)));
data.WC.chan_aligned(:,sweep_id)=temp;
end
data.WC.chan_aligned(8501:10000,:)=NaN;
% if params.flags.plot_online==1
if params.flags.WCisVoltageClamp==1
   figure; hold on
   plot(staggerplot(data.WC.chan_aligned,0,200),'k')
   axis([0 8500 -400 (params.last_sweep)*240])
else 
     figure; hold on
   plot(staggerplot(data.WC.chan_aligned,0,10),'k')
   axis([0 8500 -20 (params.last_sweep/10)*120])
end
% end
clear sweep_id temp pks plot_id thresh
%% Measure activity in analysis window
% use dVdt to find stim artifacts
for sweep_id=1:size(data.WC.chan_raw,2)
    data.WC.chan_dVdt(1:8000,sweep_id)=diff(data.WC.chan_aligned(1:8001,sweep_id)); 
end
% params.WC_stim_blank=[997:1005, 1197:1205, 1397:1405, 1597:1605, 1797:1805];%... or use defaults for 50Hz
% temp_WC=data.WC.chan_aligned; temp_WC(params.WC_stim_blank,:)=NaN; temp_WC=temp_WC(1000:2000,:);

params.WC_stim_blank=[997:1005, 1497:1505, 1997:2005, 2497:2505, 2997:3005];%... or use defaults for 20Hz
temp_WC=data.WC.chan_aligned; temp_WC(params.WC_stim_blank,:)=NaN; temp_WC=temp_WC(1000:3200,:);

if params.flags.WCisVoltageClamp==0
   data.WC.window_max       = max(data.WC.chan_aligned(params.search_win(1)*10:params.search_win(2)*10,:))-...
                             mean(data.WC.chan_aligned(params.baseline_win(1)*10:params.baseline_win(2)*10,:));
    
   data.WC.monosynaptic_max = max(data.WC.chan_aligned(params.monosynaptic_win(1)*10:params.monosynaptic_win(2)*10,:))-...
                             mean(data.WC.chan_aligned(params.baseline_win(1)*10:params.baseline_win(2)*10,:));
    
   data.WC.win_summation_max =max(temp_WC)-mean(temp_WC(10:20,:));
   for sweep_id=1:params.last_sweep
       data.WC.closestMEA_window_max_amp(sweep_id)=-1*data.burst_timing.amp{sweep_id}(data.WC.closestMEA);
   end
   
%    data.WC.win_summation_max(gt(data.WC.win_summation_max,30))=NaN; % to remove any spikes...
   
   figure; 
   subplot(1,3,1); scatter(data.WC.monosynaptic_max,         data.WC.window_max); xlabel('Monosynaptic amplitude (mV)'); ylabel('Polysynaptic amplitude (mV)'); axis([0 20 0 50])
    f=ezfit('linear');  showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
    subplot(1,3,2); scatter(data.WC.win_summation_max,        data.WC.window_max); xlabel('Maximum summation (mV)');      ylabel('Polysynaptic amplitude (mV)');axis([0 50 0 50])
    f=ezfit('linear');  showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
    subplot(1,3,3); scatter(data.WC.closestMEA_window_max_amp,data.WC.window_max); xlabel('LFP amplitude (uV)');          ylabel('Polysynaptic amplitude (mV)');axis([0 50 0 50]) 
    f=ezfit('linear');  showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
   
   
   
   
else
   data.WC.window_max       = min(data.WC.chan_aligned(params.search_win(1)*10:params.search_win(2)*10,:))-...
                             mean(data.WC.chan_aligned(params.baseline_win(1)*10:params.baseline_win(2)*10,:));

   data.WC.monosynaptic_max = min(data.WC.chan_aligned(params.monosynaptic_win(1)*10:params.monosynaptic_win(2)*10,:))-...
                            mean(data.WC.chan_aligned(params.baseline_win(1)*10:params.baseline_win(2)*10,:));
   data.WC.win_summation_max =min(temp_WC)-mean(temp_WC(10:20,:));
   for sweep_id=1:params.last_sweep
       data.WC.closestMEA_window_max_amp(sweep_id)=-1*data.burst_timing.amp{sweep_id}(data.WC.closestMEA);
   end
   figure; 
   subplot(1,2,1);
   scatter(-1*data.WC.monosynaptic_max,-1*data.WC.window_max); xlabel('W/C Monosynaptic amplitude (pA)');ylabel('W/C Polysynaptic amplitude (pA)');    axis([0 500 0 500])
   f=ezfit('linear');  showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
   subplot(1,2,2); 
% figure;
   scatter(data.WC.closestMEA_window_max_amp,-1*data.WC.window_max); xlabel('LFP burst amplitude (uV)');ylabel('W/C burst amplitude (pA)'); axis([0 500 0 500]) 
   g=ezfit('linear');  showfit(g, 'fitcolor', 'blue', 'fitlinewidth', 2);
    
end

clear temp_WC
%% find spikes and count them
if params.flags.WCisVoltageClamp==0
data.WC.APs.AP_time=cell(params.last_sweep,1);
data.WC.APs.AP_amp=cell(params.last_sweep,1);
data.WC.APs.no_APs=zeros(params.last_sweep,1);
thresh=60; % only zero-crossing events are tagged
figure;
subaxis(1,2,1, 'Spacing', 0.01, 'Padding', 0.03, 'Margin', 0.02);  
hold on;axis square
for trial_id=1:params.last_sweep
    
    [data.WC.APs.AP_amp{trial_id},data.WC.APs.AP_time{trial_id}]=...
        findpeaks(data.WC.chan_aligned(:,trial_id),'MINPEAKHEIGHT',thresh);%'THRESHOLD',1,
    
    data.WC.APs.no_APs(trial_id)=numel(data.WC.APs.AP_time{trial_id});
    plot(data.WC.chan_aligned(:,trial_id),'k')
    scatter(data.WC.APs.AP_time{trial_id},data.WC.APs.AP_amp{trial_id})
end
axis([0 8000 -20 Inf])

subaxis(1,2,2, 'Spacing', 0.01, 'Padding', 0.03, 'Margin', 0.02);  
hold on;  axis square
scatter(data.WC.closestMEA_window_max_amp,data.WC.APs.no_APs)
axis([0 500 0 10])
xlabel('LFP amplitude (uV)'); ylabel('No.APs') 
 f=ezfit('linear');  showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
clear thresh trial_id
end
%% plot closest MEA and Whole-cell alongside
% if params.flags.plot_online==1
if params.flags.WCisVoltageClamp==1
  figure; hold on
   subaxis(1,3,1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
   plot(staggerplot(squeeze(data.filtered_lfp(:,data.WC.closestMEA,:)),0,0.1),'k');
   axis([0 16000 -0.1 (params.last_sweep/10)*1.2])
   
   subaxis(1,3,2, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
   plot(sgolayfilt(staggerplot(squeeze(data.filtered_spikes(:,data.WC.closestMEA,:)),0,0.1),0,15),'r');
   axis([0 16000 -0.1 (params.last_sweep/10)*1.2])

   subaxis(1,3,3, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
   plot(sgolayfilt(staggerplot(data.WC.chan_aligned/2000,0,0.1),0,15),'k')
   axis([0 8500 -0.1 (params.last_sweep/10)*1.2])
   
else 
   figure; hold on
   subaxis(1,3,1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
   plot(staggerplot(squeeze(data.filtered_lfp(:,data.WC.closestMEA,:)),0,0.1),'k');
   axis([0 16000 -0.1 (params.last_sweep/10)*1.2])
   
   subaxis(1,3,2, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
   plot(sgolayfilt(staggerplot(squeeze(data.filtered_spikes(:,data.WC.closestMEA,:)),0,0.1),0,15),'r');
   axis([0 16000 -0.1 (params.last_sweep/10)*1.2])
   
   subaxis(1,3,3, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
   plot(staggerplot(data.WC.chan_aligned/100,0,0.1),'b')
   axis([0  8000 -0.1 (params.last_sweep/10)*1.2])
end
% end

%% vargout write
assignin('base', 'data', data) 
assignin('base', 'params', params)
end

