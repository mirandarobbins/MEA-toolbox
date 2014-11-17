%% Indexes exported data directory and iteratively loads files for measurements
%% sorts filenames into WT or KO
meta.archive.WT_filenames=cell(1);
meta.archive.KO_filenames=cell(1);
meta.archive.dir=pwd;
files=dir;
for idx=1:numel(files)
    filenames{idx,1}=files(idx).name;
    switch filenames{idx,1}(1)
        case 'W'
            if isempty (meta.archive.WT_filenames{1})
                meta.archive.WT_filenames{1}=filenames{idx,1};
            else
                meta.archive.WT_filenames{numel(meta.archive.WT_filenames)+1}=filenames{idx,1};
            end
        case 'K'
            if isempty (meta.archive.KO_filenames{1})
                meta.archive.KO_filenames{1}=filenames{idx,1};
            else
                meta.archive.KO_filenames{numel(meta.archive.KO_filenames)+1}=filenames{idx,1};
            end
        otherwise
    end
end
clear filenames files idx
%% WT load loop
tic;
no_files=numel(meta.archive.WT_filenames);
channel_index=1:64;    
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};
    % update progressbar
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    spikes = meta_spike_extract(data,params,4);
    spikes = meta_MUA64plot(data,params,spikes);
    toc
    %%%%%%%%%%% extract LFP/MUA histograms for strongest channels
    %%%%%%%%% choices for channel plotting
    temp_LFP=[]; temp_MUA=[]; temp_PDF=[];
    for sweep_id=data.sweep_sort.successful_sweeps
        meta.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
    end; clear sweep_id
    meta.StrongestChannel(meta.StrongestChannel==0)=NaN;
	
    for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
        this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
%         chosen_channel=meta.StrongestChannel(this_sweep);
    chosen_channel=mode(meta.StrongestChannel); %use mode?

    %%% 2 - choose biggest for each trial individually
%         chosen_channel=channel_index(...
%                                      min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))...
%                                      ==min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
    %%% 3 - choose earliest for each trial individually                        
%         chosen_channel=round(median(channel_index(data.burst_timing.latency{this_sweep}==min(min(data.burst_timing.latency{this_sweep})))));       
        temp_LFP(:,sweep_id)=downsample(data.filtered_lfp(:,chosen_channel,this_sweep),20);
        temp_LFP(:,std(temp_LFP)==0)=NaN; % get rid of zero-only traces
        temp_MUA(:,sweep_id)=downsample(data.filtered_spikes(:,chosen_channel,this_sweep),20);
        temp_MUA(:,std(temp_MUA)==0)=NaN; % get rid of zero-only traces
%%% for PDF
        temp_PDF(:,sweep_id)=spikes.PDF.PDF_trimmed{1,chosen_channel}(sweep_id,:);
%%% for spike histograms
%         temp_PDF(:,sweep_id)=spikes.PDF.timestamps{1,chosen_channel}(sweep_id,:);
    end
    temp_PDF(temp_PDF==0)=NaN;
    temp_PDF(:,std(temp_PDF)==0)=NaN; % get rid of zero-only traces
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP{file_idx}=temp_LFP;
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx}=temp_MUA;
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF{file_idx}=temp_PDF;
end
clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);
channel_index=1:64;    
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
    % update progressbar
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    spikes = meta_spike_extract(data,params,4);
    spikes = meta_MUA64plot(data,params,spikes);
    toc
    %%%%%%%%%%% extract LFP/MUA histograms for strongest channels
    %%%%%%%%% choices for channel plotting
    temp_LFP=[]; temp_MUA=[]; temp_PDF=[];
    for sweep_id=data.sweep_sort.successful_sweeps
        meta.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
    end; clear sweep_id
    meta.StrongestChannel(meta.StrongestChannel==0)=NaN;

    for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
        this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
%         chosen_channel=meta.StrongestChannel(this_sweep);
    chosen_channel=mode(meta.StrongestChannel); % use mode?
    %%% 2 - choose biggest for each trial individually
%         chosen_channel=channel_index(...
%                                      min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))...
%                                      ==min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
    %%% 3 - choose earliest for each trial individually                        
%         chosen_channel=round(median(channel_index(data.burst_timing.latency{this_sweep}==min(min(data.burst_timing.latency{this_sweep})))));       
        temp_LFP(:,sweep_id)=downsample(data.filtered_lfp(:,chosen_channel,this_sweep),20);
        temp_LFP(:,std(temp_LFP)==0)=NaN; % get rid of zero-only traces
        temp_MUA(:,sweep_id)=downsample(data.filtered_spikes(:,chosen_channel,this_sweep),20);
        temp_MUA(:,std(temp_MUA)==0)=NaN; % get rid of zero-only traces
%%% for PDF
        temp_PDF(:,sweep_id)=spikes.PDF.PDF_trimmed{1,chosen_channel}(sweep_id,:);
%%% for spike histograms
%         temp_PDF(:,sweep_id)=spikes.PDF.timestamps{1,chosen_channel}(sweep_id,:);
    end
    temp_PDF(temp_PDF==0)=NaN;
    temp_PDF(:,std(temp_PDF)==0)=NaN; % get rid of zero-only traces
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP{file_idx}=temp_LFP;
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx}=temp_MUA;
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF{file_idx}=temp_PDF;
end
clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% calculate population means... LFP, MUA and PDF
for file_idx=1:numel(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP)
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean(:,file_idx)=...
        nanmean(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP{file_idx},2);
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_SEM(:,file_idx)=...
        nansem(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP{file_idx},2);
    
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_pop_mean(:,file_idx)=...
        nanmean(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx},2);
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_pop_SEM(:,file_idx)=...
        nansem(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx},2);
    
    
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean(:,file_idx)=...
        nanmean(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF{file_idx},2);
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_SEM(:,file_idx)=...
        nansem(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF{file_idx},2);
end
meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean=...
    nanmean(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,2);
meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM=...
    nansem(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,2);

meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean=...
    nanmean(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,2);
meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM=...
    nansem(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,2);
for file_idx=1:numel(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP)
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean(:,file_idx)=...
        nanmean(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP{file_idx},2);
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_SEM(:,file_idx)=...
        nansem(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP{file_idx},2);
    
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_pop_mean(:,file_idx)=...
        nanmean(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx},2);
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_pop_SEM(:,file_idx)=...
        nansem(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx},2);
    
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean(:,file_idx)=...
        nanmean(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF{file_idx},2);
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_SEM(:,file_idx)=...
        nansem(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF{file_idx},2);
end

meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean=...
    nanmean(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,2);
meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM=...
    nansem(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,2);

meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean=...
    nanmean(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,2);
meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM=...
    nansem(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,2);

meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed=meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA{1};
for file_idx=2:numel(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP)
    meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed=...
        horzcat(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed,...
        meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx});
end
meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed=meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA{1};
for file_idx=2:numel(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP)
    meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed=...
        horzcat(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed,...
        meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA{file_idx});
end

%% Plot means

figure
% subplot(2,1,1); 
hold on
% plot(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,'b')
% plot(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,'r')
plot(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean,'b')
plot(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean,'r')
ciplot(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean+meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean-meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       1:1000,'b')
ciplot(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean+meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean-meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       1:1000,'r')
title('Mean LFP activity at burst centre')
xlabel('Time (ms)')   
ylabel('LFP amplitude (mV)')   

% subplot(3,2,3); hold on
% plot(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed,'b')
% % plot(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed,'r')
% axis([0 1000 -0.1 0.1])
% subplot(3,2,4); hold on
% % plot(meta.WT_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed,'b')
% plot(meta.KO_exported.mean_responses.strongest_channel.ByStrongest_MUA_collapsed,'r')
% axis([0 1000 -0.1 0.1])
%%
tb=-100:899;
figure
% subplot(2,1,2); 
hold on
% plot(tb,meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,'b')
% plot(tb,meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,'r')
plot(tb,meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean,'b','LineWidth',1.5)
plot(tb,meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean,'r','LineWidth',1.5)
upper=meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean+...
      meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower=meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean-...
      meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower(isnan(lower))=0; upper(isnan(upper))=0;

ciplot(lower,upper,tb,'b')

upper=meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean+...
      meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
  
lower=meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean-...
      meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower(isnan(lower))=0; upper(isnan(upper))=0;
ciplot(lower,upper,tb,'r')
% title('Mean spike density at burst centre')
xlabel('Post-stimulus time (ms)')   
ylabel('p(Spike/bin)')   
% axis([100 300 0 0.5])
