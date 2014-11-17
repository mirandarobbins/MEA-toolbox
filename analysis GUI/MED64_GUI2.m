function varargout = MED64_GUI2(varargin)

% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MED64_GUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @MED64_GUI2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%%%%%%%%%%%%%%%% Init CALLBACKS
function MED64_GUI2_OpeningFcn(hObject, eventdata, handles, varargin) 

handles.output = hObject; % Choose default command line output for MED64_GUI
guidata(hObject, handles); % Update handles structure
function varargout = MED64_GUI2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);


% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%% menu CALLBACKS
function FileMenu_Callback(hObject, eventdata, handles)
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

%%%%%%%%%%%%%%%% Load CALLBACKS
function LoadButton_Callback(hObject, eventdata, handles)
tic
disp('Loading........')
AD_medload();
handles.data = evalin('base','data');
handles.params = evalin('base','params');
t2=toc;     disp(strcat('Completed loading in........', num2str(t2),'s'))
displaystats(hObject, eventdata, handles)
LFP_all_overlaid(hObject, eventdata, handles) % LFP all sweeps overlaid
function LoadButtonMulti_Callback(hObject, eventdata, handles)
tic
disp('Loading........')
AD_medload_multifile2();
handles.data = evalin('base','data');
handles.params = evalin('base','params');
t2=toc;     disp(strcat('Completed loading in........', num2str(t2),'s'))
displaystats(hObject, eventdata, handles)
LFP_all_overlaid(hObject, eventdata, handles) % LFP all sweeps overlaid

function displaystats(hObject, eventdata, handles)
handles.data = evalin('base','data');
handles.params = evalin('base','params');
msgbox({'Finished analysing data!!!',...
        strcat('filename........',          handles.data.this_file)...
        strcat('no. sweeps.......',        num2str(handles.params.last_sweep)),...
        strcat('successful sweeps.......', num2str(handles.data.sweep_sort.successful_sweeps))...
        strcat('no successful sweeps.......', num2str(numel(handles.data.sweep_sort.successful_sweeps)))...
        strcat('fraction successful.......', num2str(numel(handles.data.sweep_sort.successful_sweeps)/handles.params.last_sweep))...
        strcat('Baseline window.......',   num2str(handles.params.baseline_win),'ms')...
        strcat('Search window.......',     num2str(handles.params.search_win),'ms')})
    
figure; hold on; plot(handles.data.sweep_sort.window_std_max); 
        plot ([1 handles.params.last_sweep],[handles.params.detection_threshold handles.params.detection_threshold],':k')
        axis([1 handles.params.last_sweep 0 Inf])
        xlabel('Sweep no.')
        ylabel('Normalised window std.dev.')

%%%%%%%%%%%%%%%% Plot CALLBACKS
function PlotMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot all LFP overlaid',...
                        'plot all LFP staggered',...
                        'plot all MUA overlaid',...
                        'plot all MUA staggered',...
                        'plot average LFP',...
                        'plot average LFP with minima markers',...
                        'plot LFP ampitude as heat map',...
                        'plot peak LFP latency as heat map',...
                        'plot first monosynaptic response',...
                        'plot activity on strongest channel',...
                        'plot activity on specified channel',...
                        'plot monosynaptic fEPSP amplitude as heat map',...
                        'plot monosynaptic fEPSP latency as heat map',...
                        'plot LFP for a selected trial',...
                        'plot monosynaptic CSD as heat map',...
                        'plot CSD for a chosen trial',...
                        'plot CSD all overlaid'});
function PlotMenu_Callback(hObject, eventdata, handles)
function CloseAllPB_Callback(hObject, eventdata, handles)
figs = get(0,'children');
figs(figs == gcf) = []; % delete your current figure from the list
close(figs)
function PlotButton_Callback(hObject, eventdata, handles)

handles.params = evalin('base','params');
handles.data = evalin('base','data');

popup_sel_index = get(handles.PlotMenu, 'Value');
switch popup_sel_index
    case 1
        LFP_all_overlaid(hObject, eventdata, handles)
    case 2
        LFP_all_staggered(hObject, eventdata, handles)
    case 3
        MUA_all_overlaid(hObject, eventdata, handles)
    case 4
        MUA_all_staggered(hObject, eventdata, handles)
    case 5
        average_LFP(hObject, eventdata, handles)
    case 6
        LFP_with_minima(hObject, eventdata, handles) 
    case 7
        fEPSP_heat(hObject, eventdata, handles)
    case 8
        fEPSP_heat_latency(hObject, eventdata, handles) 
    case 9
        monosynaptic_zoom(hObject, eventdata, handles)
    case 10
        PlotStrongestChannel(hObject, eventdata, handles) 
    case 11
        PlotChosenChannel(hObject, eventdata, handles)
    case 12
         mono_amp_heat(hObject, eventdata, handles)
    case 13
         mono_latency_heat(hObject, eventdata, handles)
    case 14
        LFP_selected_rep(hObject, eventdata, handles)
    case 15
        mono_CSD(hObject, eventdata, handles) 
    case 16
         CSD_selected_rep(hObject, eventdata, handles)
    case 17
        CSD_all_overlaid(hObject, eventdata, handles)
        
end

%%%%%%%%%%%%%%%% Analysis CALLBACKS
function SpikeProcessPB_Callback(hObject, eventdata, handles)
disp('caching data...')
handles.data = evalin('base','data');
handles.params = evalin('base','params');
disp('processing...')
 AD_spike_extract(handles.data,handles.params)
function MUAProcessPB_Callback(hObject, eventdata, handles)
disp('caching data...')
handles.data = evalin('base','data');
handles.params = evalin('base','params');
handles.spikes = evalin('base','spikes');
disp('processing...')
MUA64plot(handles.data,handles.params,handles.spikes)
function LFPspecProcessPB_Callback(hObject, eventdata, handles)
disp('caching data...')
handles.data = evalin('base','data');
handles.params = evalin('base','params');
disp('processing...')
% spectroexploreChronux64(handles.data,handles.params)
spectroexploreChronux64_AllChannels(handles.data,handles.params)
function CSDprocessPB_Callback(hObject, eventdata, handles)
disp('caching data...')
handles.data = evalin('base','data');
handles.params = evalin('base','params');
disp('processing...')
AD_iCSD(handles.data,handles.params)
function burstshapePB_Callback(hObject, eventdata, handles)
disp('caching data...')
handles.data = evalin('base','data');
handles.params = evalin('base','params');
handles.CSD = evalin('base','CSD');
disp('processing...')
getburstshape(handles.data,handles.params,handles.CSD)
function GrangerPB_Callback(hObject, eventdata, handles)
disp('caching data...')
disp('processing...')

%%%%%%%%%%%%%%%% movie rendering CALLBACKS
function LFPmoviePB_Callback(hObject, eventdata, handles)
disp('rendering movie...')
evalin('base','AD_MED64_LFP_movie(data)')
function CSDmoviePB_Callback(hObject, eventdata, handles)
disp('rendering movie...')
evalin('base','AD_MED64_CSD_movie(CSD,data)')

%%%%%%%%%%%%%%%% save CALLBACKS
function SaveAllPB_Callback(hObject, eventdata, handles)
disp('caching data...')
data = evalin('base','data');
params = evalin('base','params');
spikes = evalin('base','spikes');
CSD = evalin('base','CSD');
burst = evalin('base','burst');
disp('saving...')
fname=data.this_file; fname(numel(fname)-3:numel(fname))=[];
uisave({'data','params','spikes','CSD','burst'},strcat(fname,'.MAT'))
disp((strcat('saved all data as... "',fname,'.MAT"')))
function SaveMeasurementsPB_Callback(hObject, eventdata, handles)
disp('caching data...')
data = evalin('base','data');
params = evalin('base','params');
spikes = evalin('base','spikes');
CSD = evalin('base','CSD');
burst = evalin('base','burst');
% Strip big arrays
if isfield(data,'tb')
   data=rmfield(data,{'tb','raw_data','filtered_lfp','filtered_spikes'});
   CSD=rmfield(CSD,'csd_array');
end
disp('saving...')
fname=data.this_file; fname(numel(fname)-3:numel(fname))=[];
uisave({'data','params','spikes','CSD','burst'},strcat(fname,'_measurements.MAT'))
disp((strcat('saved all data as... "',fname,'.MAT"')))

                



% % --- PLOTTING OPTIONS SUBFUNCTIONS
function LFP_all_overlaid(hObject, eventdata, handles) % LFP all sweeps overlaid
figure('Name',strcat('LFP, each trial overlaid'),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
%         cla; 
        hold on,%box off; axis off
%         plot(handles.data.tb,squeeze(handles.data.filtered_lfp(:,channel_id,:)))
        plot(handles.data.tb,squeeze(handles.data.filtered_lfp(:,channel_id,handles.data.sweep_sort.successful_sweeps)),'k')
        if ~isempty(handles.data.sweep_sort.failures)
            plot(handles.data.tb,squeeze(handles.data.filtered_lfp(:,channel_id,handles.data.sweep_sort.failures)),'color',[0.3 0.3 0.3])
        end
        axis([0 1000 -0.2 0.2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,0.1,num2str(channel_id),'FontWeight','bold')
end

function LFP_selected_rep(hObject, eventdata, handles) % LFP all sweeps overlaid
  this_rep=inputdlg('Which rep to plot?','Choose a rep...');
  this_rep =str2double(this_rep);
  figure('Name',strcat('LFP, trial ',num2str(this_rep)),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
%         cla; 
        hold on,%box off; axis off
        if handles.data.sweep_sort.window_std(channel_id,this_rep)>handles.params.detection_threshold
        plot(handles.data.tb,squeeze(handles.data.filtered_lfp(:,channel_id,this_rep)),'r')
        else
        plot(handles.data.tb,squeeze(handles.data.filtered_lfp(:,channel_id,this_rep)),'color',[0.3 0.3 0.3])
        end
        axis([0 1000 -0.2 0.2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,0.1,num2str(channel_id),'FontWeight','bold')
end

function LFP_all_staggered(hObject, eventdata, handles) % LFP all sweeps overlaid
figure('Name',strcat('LFP, each trial staggered.'),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); 
        cla; hold on,%box off; axis off
        plot(sgolayfilt(staggerplot(squeeze(handles.data.filtered_lfp(:,channel_id,:)),2000,0.05),0,15));  
        axis([-1000 40000 -0.01 0.6])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,-0.08,num2str(channel_id),'FontWeight','bold')
end

function MUA_all_overlaid(hObject, eventdata, handles) %  MUA all overlaid
figure('Name',strcat('MUA, each trial overlaid.'),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); 
        cla; hold on,%box off; axis off
%         plot(handles.data.tb,sgolayfilt(squeeze(handles.data.filtered_spikes(:,channel_id,channel_id,:)),0,15))
        plot(handles.data.tb,...
            sgolayfilt(squeeze(handles.data.filtered_spikes(:,channel_id,handles.data.sweep_sort.successful_sweeps)),0,15)...
                            ,'r')
        plot(handles.data.tb,...
            sgolayfilt(squeeze(handles.data.filtered_spikes(:,channel_id,handles.data.sweep_sort.failures)),0,15)...
            ,'color',[0.3 0.3 0.3])
        axis([0 1000 -0.05 0.05])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,0.01,num2str(channel_id),'FontWeight','bold')
end

function MUA_all_staggered(hObject, eventdata, handles) % MUA all sweeps overlaid
figure('Name',strcat('MUA, each trial staggered.'),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); 
        cla; hold on,%box off; axis off
        plot(sgolayfilt(staggerplot(squeeze(handles.data.filtered_spikes(:,channel_id,:)),200,0.02),0,15));  
        axis([-1000 25000 -0.01 0.25])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,-0.08,num2str(channel_id),'FontWeight','bold')
end

function average_LFP(hObject, eventdata, handles) %  plot average LFP response
figure('Name',strcat('Average LFP response (mean +/- St.Dev.).'),'NumberTitle','off')
for plot_id=1:handles.params.no_channels
        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
        cla; hold on,%box off; axis off
        plot(handles.data.tb,handles.data.mean_channels(:,plot_id),'b','LineWidth',1)
%         ciplot((handles.data.mean_channels(:,plot_id)-handles.data.std_channels(:,plot_id)),...
%               (handles.data.mean_channels(:,plot_id)+handles.data.std_channels(:,plot_id)),handles.data.tb,'b');
%          if isvar('base','spikes')
%               evalin('caller','plot(0.5*spikes.PDF.mean_PDF(plot_id,:)-0.2)')
%          end
        axis([0 1000 -0.2 0.2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,0.05,num2str(plot_id),'FontWeight','bold')
        box on
end

function LFP_with_minima(hObject, eventdata, handles) % LFP highlighting minima
figure('Name',strcat('Mean LFP response showing peak of burst. '),'NumberTitle','off')
    for plot_id=1:handles.params.no_channels

        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
        cla; hold on,%box off; axis off
%         plot(data.mean_channels(:,plot_id),'LineWidth',1.25)
    plot(handles.data.tb,sgolayfilt(handles.data.filtered_lfp(:,plot_id,handles.params.selected_rep),handles.params.degree,handles.params.frame),'b')
%     plot(data.tb,data.filtered_lfp(:,plot_id,params.selected_rep),'k')
        plot([handles.params.search_win(1),handles.params.search_win(1)],[-0.05,0.05],'-r')    
        plot([handles.params.search_win(2),handles.params.search_win(2)],[-0.05,0.05],'-r')    
        scatter(handles.data.latency(plot_id),(handles.data.max_amp(plot_id))/1000,'or') 
        axis([0 1000 -0.2 0.2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    text(50,-0.08,num2str(plot_id),'FontWeight','bold')
    end; %clear plot_id channel_id sweep_id c

function fEPSP_heat(hObject, eventdata, handles) % plot fEPSP minima as heat map
figure('Name',strcat('Peak amplitude of LFP burst.'),'NumberTitle','off')
imagesc(handles.data.max_amp'); figure(gcf) 
cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/hot_cutoff.mat');
cmap=cmap.cmap;
colormap(cmap);
    colormap(flipud(hot));
    caxis([-400 -100]);
    c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Peak amplitude (mV)')
    title('fEPSP amplitude peak amplitude (mV)')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    clear plot_id channel_id sweep_id cs

function fEPSP_heat_latency(hObject, eventdata, handles) % plot fEPSP latency as heat map
figure('Name',strcat('latency to peak of LFP burst.'),'NumberTitle','off')
imagesc(handles.data.latency'-handles.params.first_stim); figure(gcf) 
    colormap(flipud(jet)); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'time to peak(ms)')
    title('Post-stimulus latency to minima (ms)')
    box off
    %     caxis([0 200])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
function monosynaptic_zoom(hObject, eventdata, handles) % plot fEPSP latency as heat map
figure('Name',strcat('Monosynaptic LFP (mean +/- St.Dev.).'),'NumberTitle','off')
for plot_id=1:handles.params.no_channels
        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01); hold on
        ciplot((handles.data.mean_channels(:,plot_id)-handles.data.std_channels(:,plot_id)),...
               (handles.data.mean_channels(:,plot_id)+handles.data.std_channels(:,plot_id)),handles.data.tb,'b');
        plot(handles.data.tb,handles.data.mean_channels(:,plot_id),'b','LineWidth',2)
        axis([89 119 -0.13 0.13])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
%         text(93,0.06,num2str(plot_id),'FontWeight','bold','FontSize',12,'FontName', 'Arial')
        box off; axis off
end

function PlotStrongestChannel(hObject, eventdata, handles) % plot each rep on strongest channel
figure('Name',strcat('Activity on strongest channel, each repeat staggered.'),'NumberTitle','off')
    StrongestChannel=find(handles.data.max_amp==min(min(handles.data.max_amp)));
%     StrongestChannel=inputdlg('Which channel to plot?','Choose a channel...');
%     StrongestChannel =str2double(StrongestChannel)
    x_shift=000; 
    y_shift=.2;
    temp1  =  staggerplot(squeeze(handles.data.raw_data       (:,StrongestChannel,:)),x_shift,y_shift);
    temp2  =  staggerplot(squeeze(handles.data.filtered_lfp   (:,StrongestChannel,:)),x_shift,y_shift);
    temp3 =  staggerplot(squeeze(handles.data.filtered_spikes(:,StrongestChannel,:)),x_shift,y_shift);
    % plot raw
    subplot(1,3,1); hold on
        text(100, 2.4,strcat(num2str(handles.params.last_sweep),'x repeated stimulations: Activity on strongest channel (Ch.', num2str(StrongestChannel),').')) 

        text(100,2.2,'Unfiltered Data','Fontweight','bold')
        plot((1:size(temp1,1))./20,sgolayfilt(temp1,0,15),'b')                    
        axis([-20 1020 -0.1 Inf])
        xlabel('time (ms)')%,'Fontweight','bold')
        ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    % plot filtered LFP
    subplot(1,3,2); hold on
        text(100,2.2,'Filtered LFP','Fontweight','bold')
        plot((1:size(temp2,1))./20,sgolayfilt(temp2,0,15),'b');
        axis([-20 1020 -0.1 Inf])
        xlabel('time (ms)')%,'Fontweight','bold')
        ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    % plot filtered spikes
    subplot(1,3,3), hold on
        text(100,2.2,'Filtered spikes','Fontweight','bold')
        plot((1:size(temp3,1))./20,sgolayfilt(temp3,0,15),'b');
        axis([-20 1020 -0.1 Inf])
        xlabel('time (ms)')%,'Fontweight','bold')
    ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    clear x_shift y_shift temp1 temp2 temp3     
    
function PlotChosenChannel(hObject, eventdata, handles) % plot each rep on strongest channel
    ChosenChannel=inputdlg('Which channel to plot?','Choose a channel...');
    ChosenChannel =str2double(ChosenChannel);
    figure('Name',strcat('Channel ',num2str(ChosenChannel),', each repeat staggered.'),'NumberTitle','off')
    x_shift=000; 
    y_shift=.2;
    temp1  =  staggerplot(squeeze(handles.data.raw_data       (:,ChosenChannel,:)),x_shift,y_shift);
    temp2  =  staggerplot(squeeze(handles.data.filtered_lfp   (:,ChosenChannel,:)),x_shift,y_shift);
    temp3 =  staggerplot(squeeze(handles.data.filtered_spikes(:,ChosenChannel,:)),x_shift,y_shift);
    % plot raw
    subplot(1,3,1); hold on
        text(100, 2.4,strcat(num2str(handles.params.last_sweep),'x repeated stimulations: actitivy on Ch.', num2str(ChosenChannel),').')) 

        text(100,2.2,'Unfiltered Data','Fontweight','bold')
        plot((1:size(temp1,1))./20,sgolayfilt(temp1,0,15),'b')                    
        axis([-20 1020 -0.1 Inf])
        xlabel('time (ms)')%,'Fontweight','bold')
        ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    % plot filtered LFP
    subplot(1,3,2); hold on
        text(100,2.2,'Filtered LFP','Fontweight','bold')
        plot((1:size(temp2,1))./20,sgolayfilt(temp2,0,15),'b');
        axis([-20 1020 -0.1 Inf])
        xlabel('time (ms)')%,'Fontweight','bold')
        ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    % plot filtered spikes
    subplot(1,3,3), hold on
        text(100,2.2,'Filtered spikes','Fontweight','bold')
        plot((1:size(temp3,1))./20,sgolayfilt(temp3,0,15),'b');
        axis([-20 1020 -0.1 Inf])
        xlabel('time (ms)')%,'Fontweight','bold')
    ylabel('Amplitude (\muV)')%,'Fontweight','bold')
    clear x_shift y_shift temp1 temp2 temp3     

function mono_amp_heat(hObject, eventdata, handles) % plot fEPSP minima as heat map
figure('Name',strcat('Peak amplitude of 1st monosynaptic response.'),'NumberTitle','off')
imagesc(handles.data.max_amp_mono'); figure(gcf) 
cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat');
cmap=flipud(cmap.cmap);
    colormap(cmap);
    caxis([-50 50]); 
    c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Amplitude (\muV)')
    title('Peak amplitude of first monosynaptic response')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    clear plot_id channel_id sweep_id cs

function mono_latency_heat(hObject, eventdata, handles) % plot fEPSP latency as heat map
figure('Name',strcat('latency to peak of 1st monosynaptic LFP.'),'NumberTitle','off')
imagesc(handles.data.latency_mono'-handles.params.first_stim); figure(gcf) 
    colormap(flipud(hot)); c= colorbar;%('title','fEPSP minima')
    ylabel(c,'time to peak(ms)')
    title('Post-stimulus latency to minima (ms')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 20])

function mono_CSD(hObject, eventdata, handles) % plot fEPSP latency as heat map
%% figure('Name',strcat('mean CSD of 1st monosynaptic LFP.'),'NumberTitle','off')
figure('Name',strcat('Mean current source density.'),'NumberTitle','off')
% imagesc(CSD.grand_mean_monosynaptic'); 
   handles.CSD = evalin('base','CSD');

imgdata = interp2(handles.CSD.grand_mean_monosynaptic, 4, 'spline');
% imgdata=handles.CSD.grand_mean_monosynaptic
  imagesc(imgdata')   
   
    c= colorbar;%('title','fEPSP minima')
    ylabel(c,'Scaled current source density ')
    title('Current source density (mV.mm\^-2)','Interpreter','tex')
    box off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([-1 1]*max(max(abs(handles.CSD.grand_mean_monosynaptic))))
%     cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat');
%     cmap=cmap.cmap;
cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blueblackyellow.mat');

cmap=cmap.cmap; %red is source

colormap(cmap)

function CSD_selected_rep(hObject, eventdata, handles) % CSD choose a sweep
   handles.CSD = evalin('base','CSD');
  this_rep=inputdlg('Which rep to plot?','Choose a rep...');
  this_rep =str2double(this_rep);
  figure('Name',strcat('CSD, trial ',num2str(this_rep)),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
%         cla; 
        hold on,%box off; axis off
%         plot(handles.data.tb,squeeze(handles.CSD.csd_array(:,channel_id,this_rep)))
        if handles.CSD.sweep_sort.window_std(channel_id,this_rep)>handles.CSD.sweep_sort.detection_threshold
        plot(handles.data.tb,squeeze(handles.CSD.csd_array(:,channel_id,this_rep)),'r')
        else
        plot(handles.data.tb,squeeze(handles.CSD.csd_array(:,channel_id,this_rep)),'color',[0.3 0.3 0.3])
        end
        axis([0 1000 -300 300])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,0.1,num2str(channel_id),'FontWeight','bold')
end

function CSD_all_overlaid(hObject, eventdata, handles) % CSD all sweeps overlaid
handles.CSD = evalin('base','CSD');
figure('Name',strcat('CSD, each trial overlaid'),'NumberTitle','off')
for channel_id=1:handles.params.no_channels
       

        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);  
%         cla; 
        hold on,%box off; axis off
%         plot(handles.data.tb,squeeze(handles.data.filtered_lfp(:,channel_id,:)))
        plot(handles.data.tb,squeeze(handles.CSD.csd_array(:,channel_id,handles.data.sweep_sort.successful_sweeps)),'r')
        if ~isempty(handles.data.sweep_sort.failures)
            plot(handles.data.tb,squeeze(handles.CSD.csd_array(:,channel_id,handles.data.sweep_sort.failures)),'color',[0.3 0.3 0.3])
        end
        axis([0 1000 -300 300])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        text(25,0.1,num2str(channel_id),'FontWeight','bold')
end
 

% --------------------------------------------------------------------
function DisplayStatsCB_Callback(hObject, eventdata, handles)
displaystats(hObject, eventdata, handles)


% --- Executes on button press in WC_load_PB.
function WC_load_PB_Callback(hObject, eventdata, handles)
% hObject    handle to WC_load_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
