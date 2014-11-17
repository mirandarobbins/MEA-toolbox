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
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};
%     fieldname=active_file(1:numel(active_file)-4);
%     meta.WT_exported(file_idx)=struct(fieldname,[]);
    % update progressbar
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    % optional - re-run burst extraction
    [burst] = getburstshape(data,params,CSD);
    % extract burst amp, latency, window st.dev
    meta.WT_exported.burst.max_amplitude{file_idx}    =NaN(numel(data.sweep_sort.successful_sweeps),1);
    meta.WT_exported.burst.min_latency{file_idx}      =NaN(numel(data.sweep_sort.successful_sweeps),1);
    meta.WT_exported.burst.StrongestChannel{file_idx} =NaN(numel(data.sweep_sort.successful_sweeps),1);
    meta.WT_exported.burst.window_std_max{file_idx}   =data.sweep_sort.window_std_max;
    for trial_id=data.sweep_sort.successful_sweeps
        % strongest channel
        meta.WT_exported.burst.StrongestChannel{file_idx}(trial_id)=...
                find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,trial_id))==...
                     min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,trial_id)))...
                    );
        % largest burst amplitude
        meta.WT_exported.burst.max_amplitude{file_idx}(trial_id)=min(min(data.burst_timing.amp{trial_id}));
         % earliest peak amplitude time
%         meta.WT_exported.burst.min_latency{file_idx}(trial_id)=min(min(data.burst_timing.latency{trial_id}));        
        % time of max amplitude peak
        meta.WT_exported.burst.min_latency{file_idx}(trial_id)=data.burst_timing.latency{trial_id}(meta.WT_exported.burst.StrongestChannel{file_idx}(trial_id));
    end
    meta.WT_exported.burst.StrongestChannel{file_idx}(isnan(meta.WT_exported.burst.StrongestChannel{file_idx}))=[];
    meta.WT_exported.burst.trialtotrial_stability{file_idx}=numel(unique(meta.WT_exported.burst.StrongestChannel{file_idx}))./numel(meta.WT_exported.burst.StrongestChannel{file_idx});
    % extract burst progression
    meta.WT_exported.burst.burst_window_tb=burst.burst_window_tb(1:4000);
    meta.WT_exported.burst.no_ActiveChannels_mean(1:4000,file_idx)=burst.ActiveChannels.no_ActiveChannels_mean(1:4000);
    meta.WT_exported.burst.Location.peak_displacement_vector_mean(1:4000,file_idx)=burst.Location.peak_displacement_vector_mean(1:4000);
    meta.WT_exported.burst.Location.peak_velocity_mean(1:4000,file_idx)=burst.Location.peak_velocity_mean(1:4000);
end
clear data params burst spectro spikes CSD active_file file_idx no_files progbar trial_id
clear functions
clear hidden
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
%     fieldname=active_file(1:numel(active_file)-4);
%     meta.KO_exported(file_idx)=struct(fieldname,[]);
    % update progressbar
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    
    % optional - re-run burst extraction
    [burst] = getburstshape(data,params,CSD);
    % extract burst amp, latency, window st.dev
    meta.KO_exported.burst.max_amplitude{file_idx}    =NaN(numel(data.sweep_sort.successful_sweeps),1);
    meta.KO_exported.burst.min_latency{file_idx}      =NaN(numel(data.sweep_sort.successful_sweeps),1);
    meta.KO_exported.burst.StrongestChannel{file_idx} =NaN(numel(data.sweep_sort.successful_sweeps),1);
    meta.KO_exported.burst.window_std_max{file_idx}   =data.sweep_sort.window_std_max;
    for trial_id=data.sweep_sort.successful_sweeps
        % strongest channel
        meta.KO_exported.burst.StrongestChannel{file_idx}(trial_id)=...
                find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,trial_id))==...
                     min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,trial_id)))...
                    );
        % largest burst amplitude
        meta.KO_exported.burst.max_amplitude{file_idx}(trial_id)=min(min(data.burst_timing.amp{trial_id}));
        % earliest peak amplitude time
%         meta.KO_exported.burst.min_latency{file_idx}(trial_id)=min(min(data.burst_timing.latency{trial_id}));        
        % time of max amplitude peak
        meta.KO_exported.burst.min_latency{file_idx}(trial_id)=data.burst_timing.latency{trial_id}(meta.KO_exported.burst.StrongestChannel{file_idx}(trial_id));
    end
    meta.KO_exported.burst.StrongestChannel{file_idx}(isnan(meta.KO_exported.burst.StrongestChannel{file_idx}))=[];
    meta.KO_exported.burst.trialtotrial_stability{file_idx}=numel(unique(meta.KO_exported.burst.StrongestChannel{file_idx}))./numel(meta.KO_exported.burst.StrongestChannel{file_idx});
    % extract burst progression
    meta.KO_exported.burst.burst_window_tb=burst.burst_window_tb(1:4000);
    meta.KO_exported.burst.no_ActiveChannels_mean(1:4000,file_idx)=burst.ActiveChannels.no_ActiveChannels_mean(1:4000);
    meta.KO_exported.burst.Location.peak_displacement_vector_mean(1:4000,file_idx)=burst.Location.peak_displacement_vector_mean(1:4000);
    meta.KO_exported.burst.Location.peak_velocity_mean(1:4000,file_idx)=burst.Location.peak_velocity_mean(1:4000);
end
clear data params burst spectro spikes CSD active_file file_idx no_files progbar trial_id
clear functions
clear hidden
%% calculate population means - WT
% Peak burst amplitude
meta.WT_exported.burst.max_amplitude_mean=cellfun(@nanmean,meta.WT_exported.burst.max_amplitude);
meta.WT_exported.burst.max_amplitude_sem=cellfun(@nansem,meta.WT_exported.burst.max_amplitude);
    meta.WT_exported.burst.max_amplitude_grand_mean=nanmean(meta.WT_exported.burst.max_amplitude_mean);
    meta.WT_exported.burst.max_amplitude_grand_SEM=nansem(meta.WT_exported.burst.max_amplitude_mean);

% min burst latency
meta.WT_exported.burst.min_latency_mean=cellfun(@nanmean,meta.WT_exported.burst.min_latency);
meta.WT_exported.burst.min_latency_sem=cellfun(@nansem,meta.WT_exported.burst.min_latency);
    meta.WT_exported.burst.min_latency_grand_mean=nanmean(meta.WT_exported.burst.min_latency_mean);
    meta.WT_exported.burst.min_latency_grand_SEM=nansem(meta.WT_exported.burst.min_latency_mean);
    
% trial to trial stability
meta.WT_exported.burst.trialtotrial_stability_grand_mean = nanmean(cell2mat(meta.WT_exported.burst.trialtotrial_stability));
meta.WT_exported.burst.trialtotrial_stability_grand_SEM  = nansem(cell2mat(meta.WT_exported.burst.trialtotrial_stability));
% no active channels in burst 
meta.WT_exported.burst.no_ActiveChannels_mean(meta.WT_exported.burst.no_ActiveChannels_mean==0)=NaN;
meta.WT_exported.burst.no_ActiveChannels_popMean=nanmean(meta.WT_exported.burst.no_ActiveChannels_mean,2);
meta.WT_exported.burst.no_ActiveChannels_popSEM=nansem(meta.WT_exported.burst.no_ActiveChannels_mean,2);

% no active channels in burst... normalised by peak number
for file_idx=1:numel(meta.archive.WT_filenames)
    meta.WT_exported.burst.no_ActiveChannels_mean_norm(:,file_idx)=...
        meta.WT_exported.burst.no_ActiveChannels_mean(:,file_idx)/max(meta.WT_exported.burst.no_ActiveChannels_mean(:,file_idx));
end
meta.WT_exported.burst.no_ActiveChannels_norm_popMean=nanmean(meta.WT_exported.burst.no_ActiveChannels_mean_norm,2);
meta.WT_exported.burst.no_ActiveChannels_norm_popSEM=nansem(meta.WT_exported.burst.no_ActiveChannels_mean_norm,2);

% distance from origin
meta.WT_exported.burst.Location.peak_displacement_vector_popMean=nanmean(meta.WT_exported.burst.Location.peak_displacement_vector_mean,2);
meta.WT_exported.burst.Location.peak_displacement_vector_popSEM=nansem(meta.WT_exported.burst.Location.peak_displacement_vector_mean,2);

% burst progression speed
meta.WT_exported.burst.Location.peak_velocity_popMean=nanmean(downsample(meta.WT_exported.burst.Location.peak_velocity_mean,20),2);
meta.WT_exported.burst.Location.peak_velocity_popSEM=nansem(downsample(meta.WT_exported.burst.Location.peak_velocity_mean,20),2);

%% calculate population means - KO
% Peak burst amplitude
meta.KO_exported.burst.max_amplitude_mean=cellfun(@nanmean,meta.KO_exported.burst.max_amplitude);
meta.KO_exported.burst.max_amplitude_sem=cellfun(@nansem,meta.KO_exported.burst.max_amplitude);
    meta.KO_exported.burst.max_amplitude_grand_mean=nanmean(meta.KO_exported.burst.max_amplitude_mean);
    meta.KO_exported.burst.max_amplitude_grand_SEM=nansem(meta.KO_exported.burst.max_amplitude_mean);

% min burst latency
meta.KO_exported.burst.min_latency_mean=cellfun(@nanmean,meta.KO_exported.burst.min_latency);
meta.KO_exported.burst.min_latency_sem=cellfun(@nansem,meta.KO_exported.burst.min_latency);
    meta.KO_exported.burst.min_latency_grand_mean=nanmean(meta.KO_exported.burst.min_latency_mean);
    meta.KO_exported.burst.min_latency_grand_SEM=nansem(meta.KO_exported.burst.min_latency_mean);
    
% trial to trial stability
meta.KO_exported.burst.trialtotrial_stability_grand_mean = nanmean(cell2mat(meta.KO_exported.burst.trialtotrial_stability));
meta.KO_exported.burst.trialtotrial_stability_grand_SEM  = nansem(cell2mat(meta.KO_exported.burst.trialtotrial_stability));
% no active channels in burst 
meta.KO_exported.burst.no_ActiveChannels_mean(meta.KO_exported.burst.no_ActiveChannels_mean==0)=NaN;
meta.KO_exported.burst.no_ActiveChannels_popMean=nanmean(meta.KO_exported.burst.no_ActiveChannels_mean,2);
meta.KO_exported.burst.no_ActiveChannels_popSEM=nansem(meta.KO_exported.burst.no_ActiveChannels_mean,2);

% no active channels in burst... normalised by peak number
for file_idx=1:numel(meta.archive.KO_filenames)
    meta.KO_exported.burst.no_ActiveChannels_mean_norm(:,file_idx)=...
        meta.KO_exported.burst.no_ActiveChannels_mean(:,file_idx)/max(meta.KO_exported.burst.no_ActiveChannels_mean(:,file_idx));
end
meta.KO_exported.burst.no_ActiveChannels_norm_popMean=nanmean(meta.KO_exported.burst.no_ActiveChannels_mean_norm,2);
meta.KO_exported.burst.no_ActiveChannels_norm_popSEM=nansem(meta.KO_exported.burst.no_ActiveChannels_mean_norm,2);

% distance from origin
meta.KO_exported.burst.Location.peak_displacement_vector_popMean=nanmean(meta.KO_exported.burst.Location.peak_displacement_vector_mean,2);
meta.KO_exported.burst.Location.peak_displacement_vector_popSEM=nansem(meta.KO_exported.burst.Location.peak_displacement_vector_mean,2);

% burst progression speed
meta.KO_exported.burst.Location.peak_velocity_popMean=nanmean(downsample(meta.KO_exported.burst.Location.peak_velocity_mean,20),2);
meta.KO_exported.burst.Location.peak_velocity_popSEM=nansem(downsample(meta.KO_exported.burst.Location.peak_velocity_mean,20),2);

%% Plot means
area_factor=0.15^2;
distance_factor=150;
figure; 
subplot(2,1,1);
hold on% mean distance from burst origin

ciplot((meta.WT_exported.burst.Location.peak_displacement_vector_popMean+meta.WT_exported.burst.Location.peak_displacement_vector_popSEM)*distance_factor,...
       (meta.WT_exported.burst.Location.peak_displacement_vector_popMean-meta.WT_exported.burst.Location.peak_displacement_vector_popSEM)*distance_factor,...
        meta.WT_exported.burst.burst_window_tb+100,'b');

ciplot((meta.KO_exported.burst.Location.peak_displacement_vector_popMean+meta.KO_exported.burst.Location.peak_displacement_vector_popSEM)*distance_factor,...
       (meta.KO_exported.burst.Location.peak_displacement_vector_popMean-meta.KO_exported.burst.Location.peak_displacement_vector_popSEM)*distance_factor,...
        meta.KO_exported.burst.burst_window_tb+100,'r');
plot(meta.WT_exported.burst.burst_window_tb+100,(meta.WT_exported.burst.Location.peak_displacement_vector_popMean)*distance_factor,'b');
plot(meta.KO_exported.burst.burst_window_tb+100,(meta.KO_exported.burst.Location.peak_displacement_vector_popMean)*distance_factor,'r');
% plot(meta.WT_exported.burst.burst_window_tb,meta.WT_exported.burst.no_ActiveChannels_mean,'b')
% plot(meta.KO_exported.burst.burst_window_tb,meta.KO_exported.burst.no_ActiveChannels_mean,'r')
% xlabel('Post-stimulus time (ms)')
set(gca,'XTick',[])
ylabel({'Distance from';'starting position (\mum)'})



% figure; 
% hold on% no of channels active.... normalised
% ciplot(meta.WT_exported.burst.no_ActiveChannels_norm_popMean+meta.WT_exported.burst.no_ActiveChannels_norm_popSEM,...
%        meta.WT_exported.burst.no_ActiveChannels_norm_popMean-meta.WT_exported.burst.no_ActiveChannels_norm_popSEM,...
%        meta.WT_exported.burst.burst_window_tb+100,'b');
% 
% ciplot(meta.KO_exported.burst.no_ActiveChannels_norm_popMean+meta.KO_exported.burst.no_ActiveChannels_norm_popSEM,...
%        meta.KO_exported.burst.no_ActiveChannels_norm_popMean-meta.KO_exported.burst.no_ActiveChannels_norm_popSEM,...
%        meta.KO_exported.burst.burst_window_tb+100,'r');
% plot(meta.WT_exported.burst.burst_window_tb+100,meta.WT_exported.burst.no_ActiveChannels_norm_popMean,'b');
% plot(meta.KO_exported.burst.burst_window_tb+100,meta.KO_exported.burst.no_ActiveChannels_norm_popMean,'r');
% % plot(meta.WT_exported.burst.burst_window_tb,meta.WT_exported.burst.no_ActiveChannels_mean,'b')
% % plot(meta.KO_exported.burst.burst_window_tb,meta.KO_exported.burst.no_ActiveChannels_mean,'r')
% xlabel('Post-stimulus time (ms)')
% ylabel('Normalised coactive channels')

% figure; 
subplot(2,1,2);
hold on% no of channels active
upper=(meta.WT_exported.burst.no_ActiveChannels_popMean+meta.WT_exported.burst.no_ActiveChannels_popSEM)*area_factor;
lower=(meta.WT_exported.burst.no_ActiveChannels_popMean-meta.WT_exported.burst.no_ActiveChannels_popSEM)*area_factor;
ciplot(upper,...
       lower,...
       meta.WT_exported.burst.burst_window_tb+100,'b');
   
upper=(meta.WT_exported.burst.no_ActiveChannels_popMean+meta.WT_exported.burst.no_ActiveChannels_popSEM)*area_factor;
lower=(meta.WT_exported.burst.no_ActiveChannels_popMean-meta.WT_exported.burst.no_ActiveChannels_popSEM)*area_factor;
ciplot((meta.KO_exported.burst.no_ActiveChannels_popMean+meta.KO_exported.burst.no_ActiveChannels_popSEM)*area_factor,...
       (meta.KO_exported.burst.no_ActiveChannels_popMean-meta.KO_exported.burst.no_ActiveChannels_popSEM)*area_factor,...
        meta.KO_exported.burst.burst_window_tb+100,'r');
    
plot(meta.WT_exported.burst.burst_window_tb+100,meta.WT_exported.burst.no_ActiveChannels_popMean*area_factor,'b');
plot(meta.KO_exported.burst.burst_window_tb+100,meta.KO_exported.burst.no_ActiveChannels_popMean*area_factor,'r');
% plot(meta.WT_exported.burst.burst_window_tb,meta.WT_exported.burst.no_ActiveChannels_mean,'b')
% plot(meta.KO_exported.burst.burst_window_tb,meta.KO_exported.burst.no_ActiveChannels_mean,'r')
xlabel('Post-stimulus time (ms)')
% set(gca,'XTick',[])
ylabel({'Co-active';'area (mm^2)'})
axis([100 300 0 0.12])

%%
% subplot(2,2,4);
figure;hold on% mean burst speed
ciplot(meta.WT_exported.burst.Location.peak_velocity_popMean+meta.WT_exported.burst.Location.peak_velocity_popSEM,...
       meta.WT_exported.burst.Location.peak_velocity_popMean-meta.WT_exported.burst.Location.peak_velocity_popSEM,...
       200:399,'b');
% 
ciplot(meta.KO_exported.burst.Location.peak_velocity_popMean+meta.KO_exported.burst.Location.peak_velocity_popSEM,...
       meta.KO_exported.burst.Location.peak_velocity_popMean-meta.KO_exported.burst.Location.peak_velocity_popSEM,...
       200:399,'r');
plot(200:399,(meta.WT_exported.burst.Location.peak_velocity_popMean),'b');
plot(200:399,(meta.KO_exported.burst.Location.peak_velocity_popMean),'r');
xlabel('Post-stimulus time (ms)')
ylabel('Burst propagation speed (m/s)')

%% histogram
bins=0:0.05:0.5;
meta.WT_exported.burst.Location.peak_velocity_popMean_hist=histc(meta.WT_exported.burst.Location.peak_velocity_popMean,bins)./numel(meta.WT_exported.burst.Location.peak_velocity_popMean);
meta.KO_exported.burst.Location.peak_velocity_popMean_hist=histc(meta.KO_exported.burst.Location.peak_velocity_popMean,bins)./numel(meta.KO_exported.burst.Location.peak_velocity_popMean);
figure; hold on
plot(bins,meta.WT_exported.burst.Location.peak_velocity_popMean_hist,'b','LineWidth',1.5)
plot(bins,meta.KO_exported.burst.Location.peak_velocity_popMean_hist,'r','LineWidth',1.5)
xlabel('Propagation speed (m/s)')
ylabel('Fraction of time')
axis([0 0.5 0 Inf])