function [x_corr] = filteredLFP_xcorr(data,spikes,params)
params.flags.plot_online=1;
if params.flags.plot_online==1
   plotYN=1;
else
   plotYN=0;
end

%% choices for strongest channel
x_corr.x_corr_radius=[];
channel_index=1:64;    
%%% 1 choose mode stongest channel   
for sweep_id=1:data.sweep_sort.successful_sweeps
    meta.StrongestChannel(sweep_id)=...
    find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
         min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
         );
end;
meta.StrongestChannel(meta.StrongestChannel==0)=NaN;
chosen_channel=repmat(mode(meta.StrongestChannel),1,numel(data.sweep_sort.successful_sweeps));
    
% for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
%         this_sweep=data.sweep_sort.successful_sweeps(sweep_id);   
% %% 2 - choose biggest for each trial individually
% chosen_channel(sweep_id)=channel_index(...
%                           min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))==...
%                           min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
% 
% %%% 3 - choose earliest for each trial individually                        
% % chosen_channel(sweep_id)=round(median(channel_index(data.burst_timing.latency{this_sweep}==min(min(data.burst_timing.latency{this_sweep})))));              
% end;
x_corr.centreChannel=chosen_channel;
clear sweep_id this_sweep
%% find co-active channels
x_corr.x_corr_radius.threshold=0.7;
x_corr.x_corr_radius.activeYN= zeros(64,numel(data.sweep_sort.successful_sweeps));
x_corr.x_corr_radius.active_channels=cell(numel(data.sweep_sort.successful_sweeps),1);
for trial_id=1:numel(data.sweep_sort.successful_sweeps)
    this_trial=data.sweep_sort.successful_sweeps(trial_id);
    amp_temp=min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_trial),[],1); amp_temp=amp_temp/min(amp_temp);
    x_corr.x_corr_radius.activeYN(gt(amp_temp,x_corr.x_corr_radius.threshold),trial_id)=1;
    x_corr.x_corr_radius.active_channels{trial_id}  = channel_index(gt(amp_temp,x_corr.x_corr_radius.threshold));
    % remove auto-correlation at burst centre
    x_corr.x_corr_radius.active_channels{trial_id}(x_corr.x_corr_radius.active_channels{trial_id}==chosen_channel(trial_id))=[];
end
clear trial_id this_trial amp_temp 
%% run looping cross-correlogram
x_corr.params.windowsize=100;
    x_corr.params.windowsize_ms=x_corr.params.windowsize*params.Fs^-1*1000;
x_corr.params.noverlap=20;
    x_corr.params.noverlap_ms=x_corr.params.noverlap*params.Fs^-1*1000;
x_corr.params.maxlags=500;                                                % samples
x_corr.params.lags=(-x_corr.params.maxlags:x_corr.params.maxlags)';       % samples
    x_corr.params.maxlags_ms=x_corr.params.maxlags*params.Fs^-1*1000;     % ms
    x_corr.params.lags_ms=(x_corr.params.lags)*params.Fs^-1*1000;         % ms

t1=tic;
progbar = waitbar(0,'Initializing...','name','cross-correlogram progress');

for sweep_id=1:numel(data.sweep_sort.successful_sweeps)
    this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
    %%%%% update progressbar
        waitbar(this_sweep/numel(data.sweep_sort.successful_sweeps),progbar,...
                strcat(['Crunching sweep ' num2str(sweep_id)  '/' num2str(numel(data.sweep_sort.successful_sweeps)) '...'...
                ' elepased time ', num2str(round(toc(t1)/60)), ' mins.'] ))
    %%%%% Do X-corr
    signal_1=data.filtered_lfp(:,chosen_channel(sweep_id),this_sweep); 
    for channel_id=x_corr.x_corr_radius.active_channels{sweep_id}%1:64
        signal_2=data.filtered_lfp(:,channel_id,this_sweep);
        signal_1([1994:2040, 2394:2440,  2794:2840, 3194:3240, 3594:3640]) = 0;
        signal_2([1994:2040, 2394:2440,  2794:2840, 3194:3240, 3594:3640]) = 0;
        signal_1=zscore(signal_1);signal_2=zscore(signal_2);
        [x_corr.c{channel_id,this_sweep} ...
         x_corr.L...
         x_corr.T ]     =   corrgram(signal_1,signal_2,...
                                     x_corr.params.maxlags,...
                                     x_corr.params.windowsize,...
                                     x_corr.params.noverlap);        
        x_corr.c{channel_id,this_sweep}=smooth2a(x_corr.c{channel_id,this_sweep},20,0); %2D smoothing ...,lags , time
        x_corr.max_c{channel_id,this_sweep}=zeros(1,size(x_corr.c{channel_id,this_sweep},2));
        x_corr.max_c_lagindex{channel_id,this_sweep}=zeros(2,size(x_corr.c{channel_id,this_sweep},2));
        x_corr.actual_max_lag{channel_id,this_sweep}=zeros(1,size(x_corr.c{channel_id,this_sweep},2));
        
%         for idx= 1:size(x_corr.c{channel_id,this_sweep},2)
% 
%                  x_corr.max_c{channel_id,this_sweep}(1,idx)=max(x_corr.c{channel_id,this_sweep}(:,idx));
%                  [row column]=find(x_corr.c{channel_id,this_sweep}==...
%                                    x_corr.max_c{channel_id,this_sweep}(idx),1,'first');
%                  x_corr.max_c_lagindex{channel_id,this_sweep}(1 ,idx)=row ;
%                  x_corr.max_c_lagindex{channel_id,this_sweep}(2 ,idx)=column ;         
%                  x_corr.actual_max_lag{channel_id,this_sweep}(idx)=...
%                                    x_corr.params.lags(x_corr.max_c_lagindex{channel_id,this_sweep}(1 ,idx));
%                  clear row column
%         end
    end
    waitbar(this_sweep/numel(data.sweep_sort.successful_sweeps),...
            progbar,...
            ['Finished: Processed ',num2str(numel(data.sweep_sort.successful_sweeps)),' Sweeps in ', num2str(round(toc(t1)/60)), ' mins.']);
end
close(progbar)
clear signal_2 signal_1 progbar channel_id channel_index chosen_channel t1 this_sweep sweep_id
%% Mean x-correlogram
x_corr.mean_c=[];
temp=[]; temp2=[];
% mean by trial
for sweep_id=1:numel(data.sweep_sort.successful_sweeps)
    this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
    temp=x_corr.c(:,this_sweep);
    temp(cellfun(@numel,temp)==0)=[];
    for idx=1:numel(temp)
        test=temp{idx}; test(isnan(test))=0;
        if max(var(test))==0        
            test=zeros(size(test)); test(test==0)=NaN;
        end
            temp2(:,:,idx)=temp{idx};
    end; clear idx
    x_corr.mean_c.mean_c_by_sweep(:,:,sweep_id)=nanmean(temp2,3);
end
x_corr.mean_c.mean_c_by_sweep(isnan(x_corr.mean_c.mean_c_by_sweep))=0;
x_corr.mean_c.mean_c_grand = nanmean(x_corr.mean_c.mean_c_by_sweep,3);
x_corr.mean_c.SEM_c_grand = nansem(x_corr.mean_c.mean_c_by_sweep,3);
% x_corr.mean_c.mean_c_grand=smooth2a(x_corr.mean_c.mean_c_grand,0,10); %2D smoothing ...,lags , time
clear temp temp2 sweep_id plotYN chosen_channel this_sweep
%% plot mean 
load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat')

mean_c=figure;
subplot(2,1,1)
hold on
% plot a dotted line to indicate zero lag
%     % coordinates for mean mid-line line
%         tempx=(1:size(this_cell.c{plot_id},2))/(fs/windowsize);
%         tempy=repmat(lags((ceil(numel(lags)/2))),...
%                      1,size(this_cell.mean_c.mean_c,2));
%         tempz=this_cell.mean_c.mean_c((numel(lags)-1)/2,...
%               1:size(this_cell.mean_c.mean_c,2));
%     % threshold coordinates of peak points
%         tempa=this_cell.mean_c.max_c;
%         tempa(le(tempa,0.175))=NaN;
%         plot3(tempx,tempy,tempz,'k','LineWidth',0.5)
%         plot3(tempx,...
%               (this_cell.mean_c.actual_max_lag)./10,...
%                tempa,...
%                'k','LineWidth',0.5
%     clear  tempx tempy tempz tempa

mean_c_plot=surf(x_corr.T/params.Fs,...       
                 x_corr.L,...   
                 x_corr.mean_c.mean_c_grand);

set(mean_c_plot,...                                       
    'FaceColor','interp',...
	'EdgeColor','none',...
	'FaceLighting','phong'); 

shading interp
colormap(cmap)
axis([0 max(x_corr.T)/params.Fs...
      -x_corr.params.maxlags_ms x_corr.params.maxlags_ms -0.3 0.3])
% view(0,0)
caxis([-0.5 0.5])
xlabel('time window(ms)')
ylabel('lag time (ms)')
zlabel('cross-correlation coefficient')
clear tempx tempy tempz tempa
subplot(2,1,2); hold on
    plot(x_corr.T/params.Fs,x_corr.mean_c.mean_c_grand(ceil(numel(x_corr.L)/2),:),'b','LineWidth',1.5)
    ciplot((x_corr.mean_c.mean_c_grand(ceil(numel(x_corr.L)/2),:) - x_corr.mean_c.SEM_c_grand(ceil(numel(x_corr.L)/2),:)),...
           (x_corr.mean_c.mean_c_grand(ceil(numel(x_corr.L)/2),:) + x_corr.mean_c.SEM_c_grand(ceil(numel(x_corr.L)/2),:)),...
            x_corr.T/params.Fs,'b')
   for idx=1:size(x_corr.mean_c.mean_c_by_sweep,3)
        plot(x_corr.T/params.Fs,...
            (x_corr.mean_c.mean_c_by_sweep(ceil(numel(x_corr.L)/2),:,idx)),...
            ':k','LineWidth',1)
   end; clear idx
        axis([0 1 -0.2 1])
        