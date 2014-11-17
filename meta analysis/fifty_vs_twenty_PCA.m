start_file=pwd;
cd(strcat(start_file,'/20Hz stim/Data files'))
files=dir;
reanalyse_spikes=1;

%%
meta.archive.WT_filenames=cell(1);
meta.archive.KO_filenames=cell(1);
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
%% WT load loop - 20Hz
cd(strcat(start_file,'/20Hz stim/Data files'))
tic;
no_files=numel(meta.archive.WT_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};

    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    if reanalyse_spikes==1
        spikes = meta_spike_extract(data,params,2);
        spikes = meta_MUA64plot(data,params,spikes);
    end
    toc
    [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
    
    meta.WT_exported.twenty.PCAs.MUA.trajectories{file_idx}             =   PCAs.MUA.trajectories;
    meta.WT_exported.twenty.PCAs.MUA.trajectories_mean{file_idx}        =   PCAs.MUA.trajectories_mean;
    meta.WT_exported.twenty.PCAs.MUA.trajectories_error_max{file_idx}   =   nanmax(PCAs.MUA.trajectories_error,2);
    meta.WT_exported.twenty.PCAs.MUA.trajectories_error_mean{file_idx}  =   nanmean(PCAs.MUA.trajectories_error,2);
    
    meta.WT_exported.twenty.PCAs.LFP.trajectories{file_idx}             =   PCAs.LFP.trajectories;
    meta.WT_exported.twenty.PCAs.LFP.trajectories_mean{file_idx}        =   PCAs.LFP.trajectories_mean;
    meta.WT_exported.twenty.PCAs.LFP.trajectories_error_max{file_idx}   =   nanmax(PCAs.LFP.trajectories_error,2);
    meta.WT_exported.twenty.PCAs.LFP.trajectories_error_mean{file_idx}  =   nanmean(PCAs.LFP.trajectories_error,2);
    
    meta.WT_exported.twenty.PCAs.CSD.trajectories{file_idx}             =   PCAs.CSD.trajectories;
    meta.WT_exported.twenty.PCAs.CSD.trajectories_mean{file_idx}        =   PCAs.CSD.trajectories_mean;
    meta.WT_exported.twenty.PCAs.CSD.trajectories_error_max{file_idx}   =   nanmax(PCAs.CSD.trajectories_error,2);
    meta.WT_exported.twenty.PCAs.CSD.trajectories_error_mean{file_idx}  =   nanmean(PCAs.CSD.trajectories_error,2);
end
clear data params burst spectro spikes CSD 
clear active_file file_idx no_files  CSD burst data params spikes PCAs progbar
clear functions
clear hidden
%% WT load loop - 50Hz
cd(strcat(start_file,'/50Hz stim/Data files'))
tic;
no_files=numel(meta.archive.WT_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};

    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    if reanalyse_spikes==1
        spikes = meta_spike_extract(data,params,4);
        spikes = meta_MUA64plot(data,params,spikes);
    end
    toc
    [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
    
    meta.WT_exported.fifty.PCAs.MUA.trajectories{file_idx}             =   PCAs.MUA.trajectories;
    meta.WT_exported.fifty.PCAs.MUA.trajectories_mean{file_idx}        =   PCAs.MUA.trajectories_mean;
    meta.WT_exported.fifty.PCAs.MUA.trajectories_error_max{file_idx}   =   nanmax(PCAs.MUA.trajectories_error,2);
    meta.WT_exported.fifty.PCAs.MUA.trajectories_error_mean{file_idx}  =   nanmean(PCAs.MUA.trajectories_error,2);
    
    meta.WT_exported.fifty.PCAs.LFP.trajectories{file_idx}             =   PCAs.LFP.trajectories;
    meta.WT_exported.fifty.PCAs.LFP.trajectories_mean{file_idx}        =   PCAs.LFP.trajectories_mean;
    meta.WT_exported.fifty.PCAs.LFP.trajectories_error_max{file_idx}   =   nanmax(PCAs.LFP.trajectories_error,2);
    meta.WT_exported.fifty.PCAs.LFP.trajectories_error_mean{file_idx}  =   nanmean(PCAs.LFP.trajectories_error,2);
    
    meta.WT_exported.fifty.PCAs.CSD.trajectories{file_idx}             =   PCAs.CSD.trajectories;
    meta.WT_exported.fifty.PCAs.CSD.trajectories_mean{file_idx}        =   PCAs.CSD.trajectories_mean;
    meta.WT_exported.fifty.PCAs.CSD.trajectories_error_max{file_idx}   =   nanmax(PCAs.CSD.trajectories_error,2);
    meta.WT_exported.fifty.PCAs.CSD.trajectories_error_mean{file_idx}  =   nanmean(PCAs.CSD.trajectories_error,2);
end
clear data params burst spectro spikes CSD 
clear active_file file_idx no_files  CSD burst data params spikes PCAs progbar
clear functions
clear hidden
%% KO load loop - 20Hz
cd(strcat(start_file,'/20Hz stim/Data files'))
tic;
no_files=numel(meta.archive.KO_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};

    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    if reanalyse_spikes==1
       spikes = meta_spike_extract(data,params,2);
       spikes = meta_MUA64plot(data,params,spikes);
    end
    toc
    [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
    
    meta.KO_exported.twenty.PCAs.MUA.trajectories{file_idx}             =   PCAs.MUA.trajectories;
    meta.KO_exported.twenty.PCAs.MUA.trajectories_mean{file_idx}        =   PCAs.MUA.trajectories_mean;
    meta.KO_exported.twenty.PCAs.MUA.trajectories_error_max{file_idx}   =   nanmax(PCAs.MUA.trajectories_error,2);
    meta.KO_exported.twenty.PCAs.MUA.trajectories_error_mean{file_idx}  =   nanmean(PCAs.MUA.trajectories_error,2);
    
    meta.KO_exported.twenty.PCAs.LFP.trajectories{file_idx}             =   PCAs.LFP.trajectories;
    meta.KO_exported.twenty.PCAs.LFP.trajectories_mean{file_idx}        =   PCAs.LFP.trajectories_mean;
    meta.KO_exported.twenty.PCAs.LFP.trajectories_error_max{file_idx}   =   nanmax(PCAs.LFP.trajectories_error,2);
    meta.KO_exported.twenty.PCAs.LFP.trajectories_error_mean{file_idx}  =   nanmean(PCAs.LFP.trajectories_error,2);
    
    meta.KO_exported.twenty.PCAs.CSD.trajectories{file_idx}             =   PCAs.CSD.trajectories;
    meta.KO_exported.twenty.PCAs.CSD.trajectories_mean{file_idx}        =   PCAs.CSD.trajectories_mean;
    meta.KO_exported.twenty.PCAs.CSD.trajectories_error_max{file_idx}   =   nanmax(PCAs.CSD.trajectories_error,2);
    meta.KO_exported.twenty.PCAs.CSD.trajectories_error_mean{file_idx}  =   nanmean(PCAs.CSD.trajectories_error,2);
end
clear data params burst spectro spikes CSD 
clear active_file file_idx no_files  CSD burst data params spikes PCAs progbar
clear functions
clear hidden
%% KO load loop - 50Hz
cd(strcat(start_file,'/50Hz stim/Data files'))
tic;
no_files=numel(meta.archive.KO_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};

    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    if reanalyse_spikes==1
        spikes = meta_spike_extract(data,params,4);
        spikes = meta_MUA64plot(data,params,spikes);
    end
    toc
    [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
    
    meta.KO_exported.fifty.PCAs.MUA.trajectories{file_idx}             =   PCAs.MUA.trajectories;
    meta.KO_exported.fifty.PCAs.MUA.trajectories_mean{file_idx}        =   PCAs.MUA.trajectories_mean;
    meta.KO_exported.fifty.PCAs.MUA.trajectories_error_max{file_idx}   =   nanmax(PCAs.MUA.trajectories_error,2);
    meta.KO_exported.fifty.PCAs.MUA.trajectories_error_mean{file_idx}  =   nanmean(PCAs.MUA.trajectories_error,2);
    
    meta.KO_exported.fifty.PCAs.LFP.trajectories{file_idx}             =   PCAs.LFP.trajectories;
    meta.KO_exported.fifty.PCAs.LFP.trajectories_mean{file_idx}        =   PCAs.LFP.trajectories_mean;
    meta.KO_exported.fifty.PCAs.LFP.trajectories_error_max{file_idx}   =   nanmax(PCAs.LFP.trajectories_error,2);
    meta.KO_exported.fifty.PCAs.LFP.trajectories_error_mean{file_idx}  =   nanmean(PCAs.LFP.trajectories_error,2);
    
    meta.KO_exported.fifty.PCAs.CSD.trajectories{file_idx}             =   PCAs.CSD.trajectories;
    meta.KO_exported.fifty.PCAs.CSD.trajectories_mean{file_idx}        =   PCAs.CSD.trajectories_mean;
    meta.KO_exported.fifty.PCAs.CSD.trajectories_error_max{file_idx}   =   nanmax(PCAs.CSD.trajectories_error,2);
    meta.KO_exported.fifty.PCAs.CSD.trajectories_error_mean{file_idx}  =   nanmean(PCAs.CSD.trajectories_error,2);
end
clear data params burst spectro spikes CSD 
clear active_file file_idx no_files  CSD burst data params spikes PCAs progbar
clear functions
clear hidden
cd (start_file)
%% compare LFP PCAs for 20Hz and 50Hz
include_each_rep=0;
include_means=0;
    for file_idx=2%:numel(meta.archive.KO_filenames)
        figure; hold on      
%         sp_twenty=subplot(1,2,1); hold on      
%%%%%% Plot 20Hz stim
      if include_means==1
        [X Y Z] =tubeplot(meta.KO_exported.twenty.PCAs.LFP.trajectories_mean{file_idx}(:,1),...
                          meta.KO_exported.twenty.PCAs.LFP.trajectories_mean{file_idx}(:,2),...
                          meta.KO_exported.twenty.PCAs.LFP.trajectories_mean{file_idx}(:,3),...
                          meta.KO_exported.twenty.PCAs.LFP.trajectories_error_max{file_idx});
        fig_twenty_mean=surf(X,Y,Z);
        set(fig_twenty_mean,'EdgeColor','none','FaceColor',[0 1 0],'FaceLighting','phong')
%       set(fig_twenty_mean,'FaceColor','interp')
        alpha(fig_twenty_mean,0.3) 
      end
      if include_each_rep==1
         no_reps=numel(meta.KO_exported.twenty.PCAs.LFP.trajectories{file_idx});
         for trial_id=1:no_reps
             [X Y Z] =tubeplot(meta.KO_exported.twenty.PCAs.LFP.trajectories{file_idx}{trial_id}(:,1),...
                               meta.KO_exported.twenty.PCAs.LFP.trajectories{file_idx}{trial_id}(:,2),...
                               meta.KO_exported.twenty.PCAs.LFP.trajectories{file_idx}{trial_id}(:,3),...
                               0.002);
             KO_fig=surf(X,Y,Z);
             set(KO_fig,'EdgeColor','none','FaceColor',[0 1 0],'FaceLighting','phong')
%            set(KO_fig,'FaceColor','interp')
             alpha(KO_fig,0.5) 
         end
      end
%   axis([-0.1 2 -0.6 0.6 -0.3 0.3])
    grid on
    axis vis3d
      view(3)

    camlight
    xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')

%%%%%% Plot 50Hz stim
%     sp_fifty=subplot(1,2,2); hold on      
      if include_means==1
      [X Y Z] =  tubeplot(meta.KO_exported.fifty.PCAs.LFP.trajectories_mean{file_idx}(:,1),...
                          meta.KO_exported.fifty.PCAs.LFP.trajectories_mean{file_idx}(:,2),...
                          meta.KO_exported.fifty.PCAs.LFP.trajectories_mean{file_idx}(:,3),...
                          meta.KO_exported.fifty.PCAs.LFP.trajectories_error_mean{file_idx});
        fig_fifty_mean=surf(X,Y,Z);
        set(fig_fifty_mean,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
%       set(fig_fifty_mean,'FaceColor','interp')
        alpha(fig_fifty_mean,0.5) 
      end
     
      if include_each_rep==1
         no_reps=numel(meta.KO_exported.fifty.PCAs.LFP.trajectories{file_idx});
         for trial_id=1:no_reps
             [X Y Z] =tubeplot(meta.KO_exported.fifty.PCAs.LFP.trajectories{file_idx}{trial_id}(:,1),...
                               meta.KO_exported.fifty.PCAs.LFP.trajectories{file_idx}{trial_id}(:,2),...
                               meta.KO_exported.fifty.PCAs.LFP.trajectories{file_idx}{trial_id}(:,3),...
                               0.002);
             KO_fig=surf(X,Y,Z);
             set(KO_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
%            set(KO_fig,'FaceColor','interp')
             alpha(KO_fig,0.5) 
         end
      end
      
      
%   axis([-0.1 2 -0.6 0.6 -0.3 0.3])
    grid on
    axis vis3d
  view(3)
    camlight
    xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%     hlink = linkprop([sp_twenty;sp_fifty], {'CameraPosition','CameraUpVector'});

    end
