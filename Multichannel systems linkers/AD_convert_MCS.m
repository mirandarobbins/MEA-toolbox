function [data params] = AD_convert_MCS(filename)
% APF Domanski (2014) adomansk@exseed.ed.ac.uk
%
% Linker script for Multi-electrode array analysis
%
% Uses the Neuroshare API to import Multichannel System data into MED64
% data analysis pipeline
% 
% Assumes continous sampling @20kHz with regularly spaced stimulation on 1 or
% more MEA channels
% 
% If "params.flags.AutoAlign" feature is used, time of 1st stimulus is automatically
% detected by hunting for the 1st derivative of the stim artefact
% 
% N.B. this assumes 1st stim within first 5 seconds of acquisition
%
% Timestamped MEA data is imported with timebase in "data.tb"
% 
% Imported data is aligned with 10ms baseline (i.e. stim onset at 10ms)
% and 100ms total length per sweep.
% 
% Data output is in a 3D array "data.raw_data" with format [time x channel x trial]

%% Initial set-up and choose data
params.flags.AutoAlign=1;
params.flags.PlotOnline=0;

%%%%%% load neuroshare DLL %%%%%%
[ns_RESULT]                                                = mcs_SetLibrary('nsMCDLibrary.dll');
[imported.nsresult,imported.file.Neuroshare_library_info]  = ns_GetLibraryInfo();

%%%%%% choose & open a datafile %%%%%%
if nargin==0
    [imported.file.FileName,imported.file.FilePath,FilterIndex]  = uigetfile('.mcd'); clear FilterIndex
    [ns_RESULT,imported.file.hfile]                              = ns_OpenFile(strcat(imported.file.FilePath,imported.file.FileName));
else
    imported.file.FileName = filename;
    [ns_RESULT,imported.file.hfile]                              = ns_OpenFile(imported.file.FileName);
end
[ns_RESULT,imported.file.File_overview]                      = ns_GetFileInfo(imported.file.hfile);

%%%%%% get some information %%%%%%
for entity_id=1:imported.file.File_overview.EntityCount
    [temp1,temp2]= ns_GetEntityInfo(imported.file.hfile,entity_id);
    if temp1==0
        imported.file.File_info.EntityList{entity_id,1}=temp2.EntityLabel;
        [ns_RESULT,temp] = ns_GetEntityInfo(imported.file.hfile,3);
        imported.file.File_info.ItemCount(entity_id)=temp.ItemCount;
    end
end; clear temp1 temp2 temp entity_id ns_RESULT
imported.file.File_info.MCS_EntityInfo         = mcs_Info(imported.file.hfile); %n.b. in Neuroshare, type 2 data is analogue, type 3 data is segment
imported.file.File_info.MCS_EntityNumbers      = mcs_GetEntities(imported.file.hfile,'elec0001');

%%%%%% choose 1 channel to import sampling parameters %%%%%%
[ns_RESULT,imported.file.File_info.AnalogInfo] = ns_GetAnalogInfo(imported.file.hfile,3); 

%% sneak preview first 5 seconds
if params.flags.AutoAlign==1
    temp_data=[];
    temp_cutout=[0 100000]; % import 1st 5 seconds for eyeballing
    no_samples=diff(temp_cutout);

    for ch_id=1:64% 
        temp=[]; 
        [ns_RESULT,count(ch_id,1),temp]=ns_GetAnalogData(imported.file.hfile,ch_id,...
                                                         temp_cutout(1),no_samples);
        temp_data(:,ch_id) =temp;

    end
    temp_data(:,1)=[];    % trim off the first two empty channels of trigger data

    tb_temp=(temp_cutout(1):temp_cutout(2)-1)'; tb_temp=tb_temp/20000;
    temp_data(1,:)=[]; tb_temp(1:2)=[];

    temp_data=diff(temp_data); % plot(tb_temp,temp_data)
    first_stim=tb_temp(mean(temp_data,2)==min(min(mean(temp_data,2))));
    first_stim=(first_stim-1.7e-3)/5e-5;
    % axis([0 5 -0.0001 0.0001])

else
    first_stim=inputdlg('Time of first stim in seconds','stim aligning',1,{'1.0412'});
    first_stim=str2num(first_stim{1})/5e-5;
end
disp(['first stim is at ', num2str(first_stim*5e-5) ,'s (',num2str(first_stim),' samples @20kHz).'])
params.first_stim_original=first_stim;
    clear ch_id channel_id count no_samples ns_RESULT temp trial_id temp_data temp_cutout tb_temp  first_stim
%%
%%%%%% import the sodding data... %%%%%%
temp_data=[];
% cutout=[1 imported.file.File_info.ItemCount(ch_id)]; % import all



imported.file.File_info.samples_cutout=[params.first_stim_original-200 params.first_stim_original+1800]; % samples: import 100ms epoch @ 20kHz sample rate
% imported.file.File_info.samples_cutout=[13000 15000]; % samples: import 100ms epoch @ 20kHz sample rate
imported.file.File_info.samples_advance=1200204;  % samples: (60s inter-trial interval x 20KHz)

imported.file.File_info.import_window=zeros(10,2);
for trial_id=1:10
    imported.file.File_info.import_window(trial_id,:)=(imported.file.File_info.samples_cutout+...
                                                       imported.file.File_info.samples_advance*(trial_id-1));
end
    imported.file.File_info.import_window_RealTime=imported.file.File_info.import_window/20000; %time in seconds

for ch_id=1:62% 
    for trial_id=1:10
    no_samples=imported.file.File_info.import_window(trial_id,2)-imported.file.File_info.import_window(trial_id,1);
    temp=[]; 
    [ns_RESULT,count(ch_id,1),temp]=ns_GetAnalogData(imported.file.hfile,ch_id,...
                                                     imported.file.File_info.import_window(trial_id,1),no_samples);

%     temp_data(1:imported.file.File_info.ItemCount(ch_id),ch_id)=temp(1:20000);
%     temp_data(1:imported.file.File_info.import_window(2),ch_id)=temp;
      temp_data((1:no_samples) +(no_samples*(trial_id-1)),ch_id) =temp;
    end
end
imported.file.File_info.no_samples=no_samples;
temp_data(:,1)=[];    % trim off the first two empty channels of trigger data
temp_data(:,61:64)=NaN; % add 4 empty channels (placeholders for corner data 
clear ch_id channel_id count no_samples ns_RESULT temp trial_id

% plot(reshape(temp_data(:,1),2000,10))
% axis([100 450 -0.8e-5 15e-5])
%% rearrange channel list...
% imported.file.File_info.EntityList=imported.file.File_info.EntityList(3:6
% 2);
if size(imported.file.File_info.EntityList,1)>60
imported.file.File_info.EntityList=imported.file.File_info.EntityList(2:61);
else

end
no_chars=numel(imported.file.File_info.EntityList{1}); % how many letters in entitiy string?
for channel_id=1:60
    imported.file.File_info.channel_no(channel_id,1) = ...
        str2num(imported.file.File_info.EntityList{channel_id}(no_chars-1:no_chars));
end
    
imported.file.File_info.channel_no(:,2) =1:60;
temp_order = sortrows(imported.file.File_info.channel_no);
final_data = zeros(size(temp_data)); final_data(final_data==0)=NaN;

% first row
final_data(:,2:7)=temp_data(:,temp_order(gt(temp_order(:,1),10) & lt(temp_order(:,1),19),2));
% centre rows
final_data(:,9:56)=temp_data(:,temp_order(7:54,2));
% last row
final_data(:,58:63)=temp_data(:,temp_order(55:60,2));

temp_cell=cell(8,8);
for channel_id=1:64
    temp_cell{channel_id}=final_data(:,channel_id);
end
temp_cell=flipud(temp_cell);
temp_cell=rot90(temp_cell);temp_cell=rot90(temp_cell);temp_cell=rot90(temp_cell); 

% reshape(final_data,10000,8,8);
% reshape(final_data,imported.file.File_info.ItemCount(1),8,8);
% reshape(final_data,imported.file.File_info.no_samples*10,8,8);
for channel_id=1:64
    final_data(:,channel_id)=temp_cell{channel_id};
end
clear temp temp_cell temp_data temp_order
% reshuffle data into Aleks' format
imported.file.File_info.no_trials=10;
imported.file.File_info.trial_length=imported.file.File_info.no_samples;%(diff(imported.file.File_info.import_window)+1)/imported.file.File_info.no_trials;
imported.final_data=zeros(imported.file.File_info.trial_length,64,...
                          imported.file.File_info.no_trials); 
imported.final_data(imported.final_data==0)=NaN;
for trial_id=1:imported.file.File_info.no_trials  
    imported.final_data(:,:,trial_id) = final_data(1:imported.file.File_info.trial_length,:); %final layout is a 3D array: samples x channels x trials
    final_data(1:imported.file.File_info.trial_length,:)=[];
end


% convert some more parameters -> Aleks format   
data.raw_data=imported.final_data*1e3; %n.b. convert to uV from Volts

data.this_file=imported.file.FileName;

params.Fs=imported.file.File_info.AnalogInfo.SampleRate;
params.Nyquist=params.Fs/2;
params.selected_rep=1;
params.last_sweep=size(data.raw_data,3);
params.dead_channels=[]; %channels to blank (implemented in main script)
params.no_points=size(data.raw_data,1);
params.no_channels=size(data.raw_data,2);
data.tb=((0:size(data.raw_data,1)-1)/params.Fs*1000)'; % create timebase vector in seconds

data.sweep_sort.successful_sweeps=1:10;

clear channel_id final_data no_chars trial_id
%%
% 
if params.flags.PlotOnline==1;
   figure;
   for channel_id=1:64
%        for trial_id=1:10
         for trial_id=1%:10%
            subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
            plot(data.tb,data.raw_data(:,channel_id,trial_id))
            
            axis([0 max(data.tb) -0.5 0.5]) %axes of the individual graphs
%             axis([0 inf -inf inf]) %axes of the individual graphs

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
%% clean up... run main import script
params.file=imported.file;
assignin('base', 'params', params);
assignin('base', 'data', data);
clear ans ch_id channel_id count final_data no_chars ns_RESULT temp_cell temp_data temp_order trial_id imported

[params, data] = AD_medload_MCS(data,params);
   