function [x_corr] = MUAxcorr(data,spikes,params)
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
x_corr.x_corr_radius.threshold=0.8;
x_corr.x_corr_radius.activeYN= zeros(64,numel(data.sweep_sort.successful_sweeps));
x_corr.x_corr_radius.active_channels=cell(numel(data.sweep_sort.successful_sweeps),1);
for trial_id=1:numel(data.sweep_sort.successful_sweeps)
    this_trial=data.sweep_sort.successful_sweeps(trial_id);
    amp_temp=min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_trial),[],1); amp_temp=amp_temp/min(amp_temp);
    x_corr.x_corr_radius.activeYN(gt(amp_temp,x_corr.x_corr_radius.threshold),trial_id)=1;
    x_corr.x_corr_radius.active_channels{trial_id}  = channel_index(gt(amp_temp,x_corr.x_corr_radius.threshold));
    % remove auto-correlation at burst centre
    x_corr.x_corr_radius.active_channels{trial_id}(x_corr.x_corr_radius.active_channels{trial_id}==chosen_channel(trial_id))=[]
end
clear trial_id this_trial amp_temp 
%% run looping cross-correlogram
x_corr.params.windowsize=2;
x_corr.params.noverlap=1;
x_corr.params.maxlags=50;                                                 % samples
x_corr.params.lags=(-x_corr.params.maxlags:x_corr.params.maxlags)';       % samples
    x_corr.params.maxlags_ms=x_corr.params.maxlags;                       % ms
    x_corr.params.lags_ms=(x_corr.params.lags);                           % ms

t1=tic;
progbar = waitbar(0,'Initializing...','name','cross-correlogram progress');

for sweep_id=1:numel(data.sweep_sort.successful_sweeps)
    this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
    %%%%% update progressbar
        waitbar(this_sweep/numel(data.sweep_sort.successful_sweeps),progbar,...
                strcat(['Crunching sweep ' num2str(sweep_id)  '/' num2str(numel(data.sweep_sort.successful_sweeps)) '...'...
                ' elepased time ', num2str(round(toc(t1)/60)), ' mins.'] ))
    %%%%% Do X-corr
    signal_1=spikes.PDF.PDF_trimmed{chosen_channel(sweep_id)}(sweep_id,:);
    for channel_id=x_corr.x_corr_radius.active_channels{sweep_id}%1:64
        signal_2=spikes.PDF.PDF_trimmed{channel_id}(sweep_id,:);
        signal_1(isnan(signal_1))=0;signal_2(isnan(signal_2))=0;
        [x_corr.c{channel_id,this_sweep} ...
         x_corr.L...
         x_corr.T ]     =   corrgram(signal_1,signal_2,...
                                     x_corr.params.maxlags,...
                                     x_corr.params.windowsize,...
                                     x_corr.params.noverlap);
        if range(x_corr.c{channel_id,this_sweep})==0
                 x_corr.c{channel_id,this_sweep}(x_corr.c{channel_id,this_sweep}==0)=NaN;
        end
%         x_corr.c{channel_id,this_sweep}=smooth2a(x_corr.c{channel_id,this_sweep},10, 1); %2D smoothing
        x_corr.max_c{channel_id,this_sweep}=zeros(1,size(x_corr.c{channel_id,this_sweep},2));
        x_corr.max_c_lagindex{channel_id,this_sweep}=zeros(2,size(x_corr.c{channel_id,this_sweep},2));
        x_corr.actual_max_lag{channel_id,this_sweep}=zeros(1,size(x_corr.c{channel_id,this_sweep},2));

        
        for idx= 1:size(x_corr.c{channel_id,this_sweep},2)

                 x_corr.max_c{channel_id,this_sweep}(1,idx)=max(x_corr.c{channel_id,this_sweep}(:,idx));
                 [row column]=find(x_corr.c{channel_id,this_sweep}==...
                                   x_corr.max_c{channel_id,this_sweep}(idx),1,'first');
                 if isempty(row)
                     row= NaN; column=NaN;
                 end
                 x_corr.max_c_lagindex{channel_id,this_sweep}(1 ,idx)=row ;
                 x_corr.max_c_lagindex{channel_id,this_sweep}(2 ,idx)=column ;         
                 try
                     x_corr.actual_max_lag{channel_id,this_sweep}(idx)=...
                                   x_corr.params.lags(x_corr.max_c_lagindex{channel_id,this_sweep}(1 ,idx));
                 catch
                     x_corr.actual_max_lag{channel_id,this_sweep}(idx)=NaN;
                 end
                 clear row column
        end; clear idx
    end
    waitbar(this_sweep/numel(data.sweep_sort.successful_sweeps),...
            progbar,...
            ['Finished: Processed ',num2str(numel(data.sweep_sort.successful_sweeps)),' Sweeps in ', num2str(round(toc(t1)/60)), ' mins.']);
end
disp(['Finished: Processed ',num2str(numel(data.sweep_sort.successful_sweeps)),' Sweeps in ', num2str(round(toc(t1)/60)), ' mins.'])
close(progbar);
clear signal_2 signal_1 progbar channel_id channel_index  t1 this_sweep sweep_id
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
% x_corr.mean_c.mean_c_by_sweep(isnan(x_corr.mean_c.mean_c_by_sweep))=0;
x_corr.mean_c.mean_c_grand = nanmean(x_corr.mean_c.mean_c_by_sweep,3);
x_corr.mean_c.mean_c_grand=smooth2a(x_corr.mean_c.mean_c_grand,10, 1); %2D smoothing
clear temp temp2 sweep_id plotYN chosen_channel this_sweep
%% plot mean 
load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat')

mean_c=figure;
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

mean_c_plot=surf(x_corr.T,...       
                 x_corr.L,...
                  x_corr.mean_c.mean_c_grand);

set(mean_c_plot,...                                       
    'FaceColor','interp',...
	'EdgeColor','none',...
	'FaceLighting','phong'); 

shading interp
colormap(cmap)
axis([0 999 -x_corr.params.maxlags_ms x_corr.params.maxlags_ms -1 1])
caxis([-1 1])
xlabel('time window(ms)')
ylabel('lag time (ms)')
zlabel('cross-correlation coefficient')

clear tempx tempy tempz tempa
