%% Indexes exported data directory and iteratively loads files for measurements
% WT loop
for file_idx=1:numel(meta.archive.WT_filenames)
    L4=meta.exported.WT_L4_LFP{file_idx};
    L23=meta.exported.WT_L23_LFP{file_idx};
    x_corr = pairwise_xcorr_LFP(L4,L23);

    meta.WT_exported.x_corr.mean_c_grand{file_idx}=x_corr.mean_c.mean_c_grand;
    meta.WT_exported.x_corr.SEM_c_grand{file_idx}=x_corr.mean_c.SEM_c_grand;
    meta.params=x_corr.params;
    meta.params.L=x_corr.L;
    meta.params.T=x_corr.T;
%     meta.params.amp_threshold=x_corr.x_corr_radius.threshold;
end

clear data params burst spectro spikes CSD x_corr
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% Indexes exported data directory and iteratively loads files for measurements
% KO loop
for file_idx=1:numel(meta.archive.KO_filenames)
    L4=meta.exported.KO_L4_LFP{file_idx};
    L23=meta.exported.KO_L23_LFP{file_idx};
    x_corr = pairwise_xcorr_LFP(L4,L23);

    meta.KO_exported.x_corr.mean_c_grand{file_idx}=x_corr.mean_c.mean_c_grand;
    meta.KO_exported.x_corr.SEM_c_grand{file_idx}=x_corr.mean_c.SEM_c_grand;
    meta.params=x_corr.params;
    meta.params.L=x_corr.L;
    meta.params.T=x_corr.T;
%     meta.params.amp_threshold=x_corr.x_corr_radius.threshold;
end

clear data params burst spectro spikes CSD x_corr
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden

%% calculate population means... 
for file_idx=numel(meta.archive.KO_filenames)
    temp(:,:,file_idx)=meta.WT_exported.x_corr.mean_c_grand{file_idx};
end
    meta.WT_exported.mean_responses.grand_mean_c=nanmean(temp,3);
    meta.WT_exported.mean_responses.grand_SEM_c=nansem(temp,3);
    clear temp
    
for file_idx=1:numel(meta.KO_exported.x_corr.mean_c_grand)
    temp(:,:,file_idx)=meta.KO_exported.x_corr.mean_c_grand{file_idx};
end
    meta.KO_exported.mean_responses.grand_mean_c=nanmean(temp,3);
    meta.KO_exported.mean_responses.grand_SEM_c=nansem(temp,3);
    clear temp

%% recalculate grand mean window
% take average between 0.11 and 0.26s (210->350 on timebase)
meta.WT_exported.mean_responses.grand_mean_window_redux = nanmean(meta.WT_exported.mean_responses.grand_mean_c(:,50:100),2);
meta.WT_exported.mean_responses.grand_sem_window_redux  = nanmean(meta.WT_exported.mean_responses.grand_SEM_c (:,50:100),2);
meta.KO_exported.mean_responses.grand_mean_window_redux = nanmean(meta.KO_exported.mean_responses.grand_mean_c(:,50:100),2);
meta.KO_exported.mean_responses.grand_sem_window_redux  = nanmean(meta.KO_exported.mean_responses.grand_SEM_c(:,50:100),2);

%% plot final layout
load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat')

meta.params.bin_size_in_ms=meta.params.maxlags/meta.params.windowsize;
c_lims    = [-0.1 .2];

lags_lims = 100/meta.params.bin_size_in_ms;
mean_c=figure;
subplot(2,3,1); hold on
imagesc(meta.params.T/20000-0.1,...
    meta.params.L/meta.params.bin_size_in_ms,...
    meta.WT_exported.mean_responses.grand_mean_c)

plot([-0.1,0.9],[0,0],':k')
axis([-0.1 0.9 -lags_lims lags_lims  c_lims(1) c_lims(2)])
caxis(c_lims)
xlabel('Time (s)')
ylabel('Time lag (ms)')

 colormap(hot)


subplot(2,3,2); hold on

imagesc(meta.params.T/20000-0.1,...
    meta.params.L/meta.params.bin_size_in_ms,...
    meta.KO_exported.mean_responses.grand_mean_c)
plot([-0.1,0.9],[0,0],':k')
axis([-0.1 0.9 -lags_lims lags_lims  c_lims(1) c_lims(2)])
caxis(c_lims)
xlabel('Time (s)')


subplot(2,6,5); hold on
plot(meta.params.L/meta.params.bin_size_in_ms,meta.WT_exported.mean_responses.grand_mean_window_redux,'b','LineWidth',1.2)
% plot(meta.params.window.lags'/20000*1000,meta.WT_exported.x_corr.window.mean_c.grand_mean_c,'b')
ciplot(meta.WT_exported.mean_responses.grand_mean_window_redux - meta.WT_exported.mean_responses.grand_sem_window_redux,...
    meta.WT_exported.mean_responses.grand_mean_window_redux + meta.WT_exported.mean_responses.grand_sem_window_redux,...
    meta.params.L/meta.params.bin_size_in_ms,'b')

plot(meta.params.L/meta.params.bin_size_in_ms,meta.KO_exported.mean_responses.grand_mean_window_redux,'r','LineWidth',1.2)
% plot(meta.params.window.lags'/20000*1000,meta.KO_exported.x_corr.window.mean_c.grand_mean_c,'r')
ciplot(meta.KO_exported.mean_responses.grand_mean_window_redux - meta.KO_exported.mean_responses.grand_sem_window_redux,...
    meta.KO_exported.mean_responses.grand_mean_window_redux + meta.KO_exported.mean_responses.grand_sem_window_redux,...
    meta.params.L/meta.params.bin_size_in_ms,'r')
plot([0,0],[-0,.25],':k')
axis([-50 50 -0.1 .3])
%         set(gca,'Xscale','log')
set(gca,'xtick',[])

% xlabel('Time lag (ms)')
ylabel('Correlation')
view(90,-90)


subplot(2,1,2); hold on
plot(meta.params.T/20000-0.1,...
    smooth_hist(meta.WT_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:)),'b','LineWidth',1.5)

plot(meta.params.T/20000-0.1,...
    smooth_hist(meta.KO_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:)),'r','LineWidth',1.5)

ciplot(smooth_hist(meta.WT_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) - meta.WT_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...
    smooth_hist(meta.WT_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) + meta.WT_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...
    meta.params.T/20000-0.1,'b')
ciplot(smooth_hist(meta.KO_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) - meta.KO_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...
    smooth_hist(meta.KO_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) + meta.KO_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...
    meta.params.T/20000-0.1,'r')
axis([-.1 0.9 0 0.5])
xlabel('time (s)')
ylabel('Zero-lag correlation')