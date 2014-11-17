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

progbar_MAIN = waitbar(0,'Initializing...','name','file loading progress');

for file_idx=8:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};
    % update progressbar
    waitbar(file_idx/no_files,progbar_MAIN,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
%     spikes = meta_spike_extract(data,params,4);
%     spikes = meta_MUA64plot(data,params,spikes);
    x_corr = filteredLFP_xcorr(data,spikes,params);
    toc
    meta.WT_exported.x_corr.mean_c_grand{file_idx}=x_corr.mean_c.mean_c_grand;
    meta.WT_exported.x_corr.SEM_c_grand{file_idx}=x_corr.mean_c.SEM_c_grand;
    meta.params=x_corr.params;
    meta.params.L=x_corr.L;
    meta.params.T=x_corr.T;
    meta.params.amp_threshold=x_corr.x_corr_radius.threshold;
end

clear data params burst spectro spikes CSD x_corr
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);

    channel_index=1:64;    

progbar_MAIN = waitbar(0,'Initializing...','name','file loading progress');

for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
    % update progressbar
    waitbar(file_idx/no_files,progbar_MAIN,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
%     spikes = meta_spike_extract(data,params,4);
%     spikes = meta_MUA64plot(data,params,spikes);
    x_corr = filteredLFP_xcorr(data,spikes,params);
    toc
    meta.KO_exported.x_corr.mean_c_grand{file_idx}=x_corr.mean_c.mean_c_grand;
    meta.KO_exported.x_corr.SEM_c_grand{file_idx}=x_corr.mean_c.SEM_c_grand;
    meta.params=x_corr.params;
    meta.params.L=x_corr.L;
    meta.params.T=x_corr.T;
    meta.params.amp_threshold=x_corr.x_corr_radius.threshold;
end

clear data params burst spectro spikes CSD x_corr
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden

%% calculate population means... 
for file_idx=1:14%numel(meta.WT_exported.x_corr.mean_c_grand)
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

%% plot mean 
load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat')
c_lims=[-0.3 0.3];
mean_c=figure;
subplot(3,1,1); hold on
mean_c_plot=surf(meta.params.T/20000,...       
                 meta.params.L,...   
                 meta.WT_exported.mean_responses.grand_mean_c);
set(mean_c_plot,...                                       
    'FaceColor','interp',...
	'EdgeColor','none',...
	'FaceLighting','phong'); 
    shading interp
    colormap(cmap)
    axis([0 max(meta.params.T)/20000 -meta.params.maxlags_ms meta.params.maxlags_ms c_lims(1) c_lims(2)])
    view(0,0)
    caxis(c_lims)
    xlabel('time window(ms)')
    ylabel('lag time (ms)')
    zlabel('cross-correlation coefficient')

    
subplot(3,1,2); hold on
mean_c_plot=surf(meta.params.T/20000,...       
                 meta.params.L,...   
                 meta.KO_exported.mean_responses.grand_mean_c);
set(mean_c_plot,...                                       
    'FaceColor','interp',...
	'EdgeColor','none',...
	'FaceLighting','phong'); 
    shading interp
    colormap(cmap)
    axis([0 max(meta.params.T)/20000 -meta.params.maxlags_ms meta.params.maxlags_ms c_lims(1) c_lims(2)])
    view(0,0)
    caxis(c_lims)
    xlabel('time window(ms)')
    ylabel('lag time (ms)')
    zlabel('cross-correlation coefficient')    
    

subplot(3,1,3); hold on
    plot(meta.params.T/20000,...
         meta.WT_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:),'b','LineWidth',1.5)
    
    plot(meta.params.T/20000,...
         meta.KO_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:),'r','LineWidth',1.5)
    
    ciplot((meta.WT_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) - meta.WT_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...
           (meta.WT_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) + meta.WT_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...  
           meta.params.T/20000,'b')
    ciplot((meta.KO_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) - meta.KO_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...
           (meta.KO_exported.mean_responses.grand_mean_c(ceil(numel(meta.params.L)/2),:) + meta.KO_exported.mean_responses.grand_SEM_c(ceil(numel(meta.params.L)/2),:)),...  
           meta.params.T/20000,'r')
       axis([0 1 0 0.8])
    xlabel('time window (ms)')
    ylabel('zero-lag correlation')
%    for  idx=1:numel(meta.WT_exported.x_corr.mean_c_grand)
%         plot(meta.params.T/20000,...
%             meta.WT_exported.x_corr.mean_c_grand{idx}(ceil(numel(meta.params.L)/2),:),...
%             ':b','LineWidth',1)
%    end
%    for  idx=1:numel(meta.KO_exported.x_corr.mean_c_grand)
%         plot(meta.params.T/20000,...
%             meta.KO_exported.x_corr.mean_c_grand{idx}(ceil(numel(meta.params.L)/2),:),...
%             ':r','LineWidth',1)
%    end;
%    clear idx

        