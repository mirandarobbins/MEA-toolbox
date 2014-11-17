%% Indexes exported data directory and iteratively loads files for measurements
reanalyse_spikes=0;
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
    if reanalyse_spikes==1
        spikes = AD_spike_extract(data,params);
        spikes=meta_spike_extract(data,params,3);
        spikes = meta_MUA64plot(data,params,spikes);
    end
    toc
    [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
    clear active_file CSD burst data params spikes
    meta.WT_exported.PCAs.MUA.percent_var_explained{file_idx}    =   PCAs.MUA.percent_var_explained;
    meta.WT_exported.PCAs.MUA.trajectories{file_idx}             =   PCAs.MUA.trajectories;
    meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}        =   PCAs.MUA.trajectories_mean;
    meta.WT_exported.PCAs.MUA.trajectories_error_max{file_idx}   =   nanmax(PCAs.MUA.trajectories_error,2);
    meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}  =   nanmean(PCAs.MUA.trajectories_error,2);
    
    meta.WT_exported.PCAs.LFP.percent_var_explained{file_idx}    =   PCAs.LFP.percent_var_explained;
    meta.WT_exported.PCAs.LFP.trajectories{file_idx}             =   PCAs.LFP.trajectories;
    meta.WT_exported.PCAs.LFP.trajectories_mean{file_idx}        =   PCAs.LFP.trajectories_mean;
    meta.WT_exported.PCAs.LFP.trajectories_error_max{file_idx}   =   nanmax(PCAs.LFP.trajectories_error,2);
    meta.WT_exported.PCAs.LFP.trajectories_error_mean{file_idx}  =   nanmean(PCAs.LFP.trajectories_error,2);
    
    meta.WT_exported.PCAs.CSD.percent_var_explained{file_idx}    =   PCAs.CSD.percent_var_explained;
    meta.WT_exported.PCAs.CSD.trajectories{file_idx}             =   PCAs.CSD.trajectories;
    meta.WT_exported.PCAs.CSD.trajectories_mean{file_idx}        =   PCAs.CSD.trajectories_mean;
    meta.WT_exported.PCAs.CSD.trajectories_error_max{file_idx}   =   nanmax(PCAs.CSD.trajectories_error,2);
    meta.WT_exported.PCAs.CSD.trajectories_error_mean{file_idx}  =   nanmean(PCAs.CSD.trajectories_error,2);
    clear PCAs
end
clear file_idx no_files
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
    if reanalyse_spikes==1
        %         spikes = AD_spike_extract(data,params);
        spikes=meta_spike_extract(data,params,3);
        spikes = meta_MUA64plot(data,params,spikes);
    end
    toc
    [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
    clear active_file CSD burst data params spikes
    meta.KO_exported.PCAs.MUA.percent_var_explained{file_idx}    =   PCAs.MUA.percent_var_explained;
    meta.KO_exported.PCAs.MUA.trajectories{file_idx}             =   PCAs.MUA.trajectories;
    meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}        =   PCAs.MUA.trajectories_mean;
    meta.KO_exported.PCAs.MUA.trajectories_error_max{file_idx}   =   nanmax(PCAs.MUA.trajectories_error,2);
    meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx}  =   nanmean(PCAs.MUA.trajectories_error,2);
    
    meta.KO_exported.PCAs.LFP.percent_var_explained{file_idx}    =   PCAs.LFP.percent_var_explained;
    meta.KO_exported.PCAs.LFP.trajectories{file_idx}             =   PCAs.LFP.trajectories;
    meta.KO_exported.PCAs.LFP.trajectories_mean{file_idx}        =   PCAs.LFP.trajectories_mean;
    meta.KO_exported.PCAs.LFP.trajectories_error_max{file_idx}   =   nanmax(PCAs.LFP.trajectories_error,2);
    meta.KO_exported.PCAs.LFP.trajectories_error_mean{file_idx}  =   nanmean(PCAs.LFP.trajectories_error,2);
    
    meta.KO_exported.PCAs.CSD.percent_var_explained{file_idx}    =   PCAs.CSD.percent_var_explained;
    meta.KO_exported.PCAs.CSD.trajectories{file_idx}             =   PCAs.CSD.trajectories;
    meta.KO_exported.PCAs.CSD.trajectories_mean{file_idx}        =   PCAs.CSD.trajectories_mean;
    meta.KO_exported.PCAs.CSD.trajectories_error_max{file_idx}   =   nanmax(PCAs.CSD.trajectories_error,2);
    meta.KO_exported.PCAs.CSD.trajectories_error_mean{file_idx}  =   nanmean(PCAs.CSD.trajectories_error,2);
    clear PCAs
end
clear active_file file_idx no_files
clear functions
clear hidden
%% extract pop data
for file_idx=1:numel(meta.archive.WT_filenames)
    meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative{file_idx}=cumsum(meta.WT_exported.PCAs.MUA.trajectories_error_max{file_idx});
    meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative{file_idx}=cumsum(meta.WT_exported.PCAs.LFP.trajectories_error_max{file_idx});
    meta.WT_exported.PCAs.CSD.trajectories_error_max_cumulative{file_idx}=cumsum(meta.WT_exported.PCAs.CSD.trajectories_error_max{file_idx});
    
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx});
    meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.WT_exported.PCAs.LFP.trajectories_error_mean{file_idx});
    meta.WT_exported.PCAs.CSD.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.WT_exported.PCAs.CSD.trajectories_error_mean{file_idx});
end

for file_idx=1:numel(meta.archive.KO_filenames)
    meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.MUA.trajectories_error_max{file_idx});
    meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.LFP.trajectories_error_max{file_idx});
    meta.KO_exported.PCAs.CSD.trajectories_error_max_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.CSD.trajectories_error_max{file_idx});
    
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx});
    meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.LFP.trajectories_error_mean{file_idx});
    meta.KO_exported.PCAs.CSD.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.CSD.trajectories_error_mean{file_idx});
end
% grand average errors in first 3 dims...
%MUA
meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean  = nanmean(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative),2);
meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_SEM   = nansem(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative),2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean = nanmean(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative),2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM  = nansem(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative),2);

meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean  = nanmean(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative),2);
meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_SEM   = nansem(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative),2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean = nanmean(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative),2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM  = nansem(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative),2);

%LFP
meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean  = nanmean(cell2mat(meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative),2);
meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_SEM   = nansem(cell2mat(meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative),2);
meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean = nanmean(cell2mat(meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative),2);
meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_SEM  = nansem(cell2mat(meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative),2);

meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean  = nanmean(cell2mat(meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative),2);
meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_SEM   = nansem(cell2mat(meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative),2);
meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean = nanmean(cell2mat(meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative),2);
meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_SEM  = nansem(cell2mat(meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative),2);

%CSD
meta.WT_exported.PCAs.CSD.trajectories_error_max_cumulative_grand_mean  = nanmean(cell2mat(meta.WT_exported.PCAs.CSD.trajectories_error_max_cumulative),2);
meta.WT_exported.PCAs.CSD.trajectories_error_max_cumulative_grand_SEM   = nansem(cell2mat(meta.WT_exported.PCAs.CSD.trajectories_error_max_cumulative),2);
meta.WT_exported.PCAs.CSD.trajectories_error_mean_cumulative_grand_mean = nanmean(cell2mat(meta.WT_exported.PCAs.CSD.trajectories_error_mean_cumulative),2);
meta.WT_exported.PCAs.CSD.trajectories_error_mean_cumulative_grand_SEM  = nansem(cell2mat(meta.WT_exported.PCAs.CSD.trajectories_error_mean_cumulative),2);

meta.KO_exported.PCAs.CSD.trajectories_error_max_cumulative_grand_mean  = nanmean(cell2mat(meta.KO_exported.PCAs.CSD.trajectories_error_max_cumulative),2);
meta.KO_exported.PCAs.CSD.trajectories_error_max_cumulative_grand_SEM   = nansem(cell2mat(meta.KO_exported.PCAs.CSD.trajectories_error_max_cumulative),2);
meta.KO_exported.PCAs.CSD.trajectories_error_mean_cumulative_grand_mean = nanmean(cell2mat(meta.KO_exported.PCAs.CSD.trajectories_error_mean_cumulative),2);
meta.KO_exported.PCAs.CSD.trajectories_error_mean_cumulative_grand_SEM  = nansem(cell2mat(meta.KO_exported.PCAs.CSD.trajectories_error_mean_cumulative),2);

%% inspect variance explained by each PCA - pseudo-pareto plot
% MUA
meta.WT_exported.PCAs.MUA.percent_var_explained_mean     = nanmean(cell2mat(meta.WT_exported.PCAs.MUA.percent_var_explained),2);
meta.WT_exported.PCAs.MUA.percent_var_explained_SEM      = nansem(cell2mat(meta.WT_exported.PCAs.MUA.percent_var_explained),2);
meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean = nanmean(cumsum(cell2mat(meta.WT_exported.PCAs.MUA.percent_var_explained)),2);
meta.WT_exported.PCAs.MUA.percent_var_explained_cum_SEM  = nansem(cumsum(cell2mat(meta.WT_exported.PCAs.MUA.percent_var_explained)),2);

meta.KO_exported.PCAs.MUA.percent_var_explained_mean     = nanmean(cell2mat(meta.KO_exported.PCAs.MUA.percent_var_explained),2);
meta.KO_exported.PCAs.MUA.percent_var_explained_SEM      = nansem(cell2mat(meta.KO_exported.PCAs.MUA.percent_var_explained),2);
meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean = nanmean(cumsum(cell2mat(meta.KO_exported.PCAs.MUA.percent_var_explained)),2);
meta.KO_exported.PCAs.MUA.percent_var_explained_cum_SEM  = nansem(cumsum(cell2mat(meta.KO_exported.PCAs.MUA.percent_var_explained)),2);
% LFP
meta.WT_exported.PCAs.LFP.percent_var_explained_mean     = nanmean(cell2mat(meta.WT_exported.PCAs.LFP.percent_var_explained),2);
meta.WT_exported.PCAs.LFP.percent_var_explained_SEM      = nansem(cell2mat(meta.WT_exported.PCAs.LFP.percent_var_explained),2);
meta.WT_exported.PCAs.LFP.percent_var_explained_cum_mean = nanmean(cumsum(cell2mat(meta.WT_exported.PCAs.LFP.percent_var_explained)),2);
meta.WT_exported.PCAs.LFP.percent_var_explained_cum_SEM  = nansem(cumsum(cell2mat(meta.WT_exported.PCAs.LFP.percent_var_explained)),2);

meta.KO_exported.PCAs.LFP.percent_var_explained_mean     = nanmean(cell2mat(meta.KO_exported.PCAs.LFP.percent_var_explained),2);
meta.KO_exported.PCAs.LFP.percent_var_explained_SEM      = nansem(cell2mat(meta.KO_exported.PCAs.LFP.percent_var_explained),2);
meta.KO_exported.PCAs.LFP.percent_var_explained_cum_mean = nanmean(cumsum(cell2mat(meta.KO_exported.PCAs.LFP.percent_var_explained)),2);
meta.KO_exported.PCAs.LFP.percent_var_explained_cum_SEM  = nansem(cumsum(cell2mat(meta.KO_exported.PCAs.LFP.percent_var_explained)),2);

figure; hold on
subplot(1,2,1); hold on
%plot absolute
ciplot(meta.WT_exported.PCAs.MUA.percent_var_explained_mean+meta.WT_exported.PCAs.MUA.percent_var_explained_SEM,...
    meta.WT_exported.PCAs.MUA.percent_var_explained_mean-meta.WT_exported.PCAs.MUA.percent_var_explained_SEM,1:64,'b')
plot(meta.WT_exported.PCAs.MUA.percent_var_explained_mean,'b','LineWidth',1.5)
ciplot(meta.KO_exported.PCAs.MUA.percent_var_explained_mean+meta.KO_exported.PCAs.MUA.percent_var_explained_SEM,...
    meta.KO_exported.PCAs.MUA.percent_var_explained_mean-meta.KO_exported.PCAs.MUA.percent_var_explained_SEM,1:64,'r')
plot(meta.KO_exported.PCAs.MUA.percent_var_explained_mean,'r','LineWidth',1.5)


%plot cumulative
ciplot(meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean+meta.WT_exported.PCAs.MUA.percent_var_explained_cum_SEM,...
    meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean-meta.WT_exported.PCAs.MUA.percent_var_explained_cum_SEM,1:64,'b')
plot(meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean,':b','LineWidth',1)

ciplot(meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean+meta.KO_exported.PCAs.MUA.percent_var_explained_cum_SEM,...
    meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean-meta.KO_exported.PCAs.MUA.percent_var_explained_cum_SEM,1:64,'r')
plot(meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean,':r','LineWidth',1)
% text(7, 10, 'MUA','FontSize', 18)
axis([1 5 0 100])
xlabel('Principle component no.')
ylabel('% Variance explained.')

subplot(1,2,2); hold on
%plot absolute
ciplot(meta.WT_exported.PCAs.LFP.percent_var_explained_mean+meta.WT_exported.PCAs.LFP.percent_var_explained_SEM,...
    meta.WT_exported.PCAs.LFP.percent_var_explained_mean-meta.WT_exported.PCAs.LFP.percent_var_explained_SEM,1:64,'b')
plot(meta.WT_exported.PCAs.LFP.percent_var_explained_mean,'b','LineWidth',1.5)
ciplot(meta.KO_exported.PCAs.LFP.percent_var_explained_mean+meta.KO_exported.PCAs.LFP.percent_var_explained_SEM,...
    meta.KO_exported.PCAs.LFP.percent_var_explained_mean-meta.KO_exported.PCAs.LFP.percent_var_explained_SEM,1:64,'r')
plot(meta.KO_exported.PCAs.LFP.percent_var_explained_mean,'r','LineWidth',1.5)


%plot cumulative
ciplot(meta.WT_exported.PCAs.LFP.percent_var_explained_cum_mean+meta.WT_exported.PCAs.LFP.percent_var_explained_cum_SEM,...
    meta.WT_exported.PCAs.LFP.percent_var_explained_cum_mean-meta.WT_exported.PCAs.LFP.percent_var_explained_cum_SEM,1:64,'b')
plot(meta.WT_exported.PCAs.LFP.percent_var_explained_cum_mean,':b','LineWidth',1)

ciplot(meta.KO_exported.PCAs.LFP.percent_var_explained_cum_mean+meta.KO_exported.PCAs.LFP.percent_var_explained_cum_SEM,...
    meta.KO_exported.PCAs.LFP.percent_var_explained_cum_mean-meta.KO_exported.PCAs.LFP.percent_var_explained_cum_SEM,1:64,'r')
plot(meta.KO_exported.PCAs.LFP.percent_var_explained_cum_mean,':r','LineWidth',1)
text(7, 10, 'LFP','FontSize', 18)
axis([0 10 0 100])
xlabel('Principle component no.')
ylabel('% Variance explained.')
%% plot max error
tb=0:999;
figure; hold on
subplot(1,2,1), hold on
for file_idx=1:numel(meta.archive.WT_filenames)
    plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative{file_idx}./...
        max(meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative{file_idx}),'b')
end
for file_idx=1:numel(meta.archive.KO_filenames)
    plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative{file_idx}./...
        max(meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative{file_idx}),'r')
end
ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_SEM,...
    meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_SEM,...
    tb,'b')
plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean,'b','LineWidth',1.5)


ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_SEM,...
    meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_SEM,...
    tb,'r')
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_max_cumulative_grand_mean,'r','LineWidth',1.5)


% text(850, 0.06, 'MUA','FontSize', 18)
xlabel('Time (ms)')
ylabel('Cumulative variance in window')
% axis([200 800 0 1])

subplot(1,2,2), hold on
% figure, hold on
% for file_idx=1:numel(meta.archive.WT_filenames)
%     plot(tb,meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative{file_idx}./max(meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative{file_idx}),'b')
% end
% for file_idx=1:numel(meta.archive.KO_filenames)
%     plot(tb,meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative{file_idx}./max(meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative{file_idx}),'r')
% end
ciplot(meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean+meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_SEM,...
    meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean-meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_SEM,...
    tb,'b')
plot(tb,meta.WT_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean,'b','LineWidth',1.5)


ciplot(meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean+meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_SEM,...
    meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean-meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_SEM,...
    tb,'r')
plot(tb,meta.KO_exported.PCAs.LFP.trajectories_error_max_cumulative_grand_mean,'r','LineWidth',1.5)

text(850, 0.06, 'LFP','FontSize', 18)
xlabel('Time (ms)')
ylabel('Cumulative error')
% axis([200 800 0 1])
%% plot mean error
tb=-100:899;
figure; hold on
subplot(1,2,1), hold on
% for file_idx=1:numel(meta.archive.WT_filenames)
%     plot(200:1000,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx},'b')%./max(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})
% end
% for file_idx=1:numel(meta.archive.KO_filenames)
%     plot(200:1000,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx},'r')%./max(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})
% end
ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    tb,'b')
plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean,'b','LineWidth',1.5)


ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    tb,'r')
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean,'r','LineWidth',1.5)


% text(800, 2, 'MUA','FontSize', 18)
xlabel('Post stimulus time (ms)')
ylabel('Cumulative variance (A.U.)')
% axis([100 300 0 20])

subplot(1,2,2), hold on
% for file_idx=1:numel(meta.archive.WT_filenames)
%     plot(200:1000,meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative{file_idx}./max(meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative{file_idx}),'b')
% end
% for file_idx=1:numel(meta.archive.KO_filenames)
%     plot(200:1000,meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative{file_idx}./max(meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative{file_idx}),'r')
% end
ciplot(meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean+meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_SEM,...
    meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean-meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_SEM,...
    tb,'b')
plot(tb,meta.WT_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean,'b','LineWidth',1.5)


ciplot(meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean+meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_SEM,...
    meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean-meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_SEM,...
    tb,'r')
plot(tb,meta.KO_exported.PCAs.LFP.trajectories_error_mean_cumulative_grand_mean,'r','LineWidth',1.5)

text(800, 2, 'LFP','FontSize', 18)
xlabel('Time (ms)')
ylabel('Cumulative error')
axis([100 1000 0 20])
%% plot mean trajectories - MUA
figure; hold on
for file_idx=1:numel(meta.archive.WT_filenames)
    [X Y Z] =tubeplot(meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
        meta.WT_exported.PCAs.MUA.trajectories_error_max{file_idx}/10);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(WT_fig,'FaceColor','interp')
    alpha(WT_fig,0.3)
end

for file_idx=1:numel(meta.archive.KO_filenames)
    [X Y Z] =tubeplot(meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
        meta.KO_exported.PCAs.MUA.trajectories_error_max{file_idx}/10);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(KO_fig,'FaceColor','interp')
end
grid on
view(3)
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%% plot mean trajectories - LFP
win=200:1000;
figure; hold on
% sp1=subplot(1,2,1);hold on
cmap=hsv(numel(meta.archive.WT_filenames));
for file_idx=1:numel(meta.archive.WT_filenames)
    [X Y Z] =tubeplot(meta.WT_exported.PCAs.LFP.trajectories_mean{file_idx}(win,1),...
        meta.WT_exported.PCAs.LFP.trajectories_mean{file_idx}(win,2),...
        meta.WT_exported.PCAs.LFP.trajectories_mean{file_idx}(win,3),...
        meta.WT_exported.PCAs.LFP.trajectories_error_max{file_idx}/10);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(WT_fig,'FaceColor','interp')
end
% axis([-0.5 0.5 -0.3 0.2 -0.1 0.1])
grid on
% view(3)
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')

% sp2=subplot(1,2,2);hold on
for file_idx=1:numel(meta.archive.KO_filenames)
    [X Y Z] =tubeplot(meta.KO_exported.PCAs.LFP.trajectories_mean{file_idx}(win,1),...
        meta.KO_exported.PCAs.LFP.trajectories_mean{file_idx}(win,2),...
        meta.KO_exported.PCAs.LFP.trajectories_mean{file_idx}(win,3),...
        meta.KO_exported.PCAs.LFP.trajectories_error_max{file_idx}/10);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(KO_fig,'FaceColor','interp')
end

camlight
axis vis3d

xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
% axis([-0.5 0.5 -0.3 0.2 -0.1 0.1])
% hlink = linkprop([sp1;sp2], {'CameraPosition','PlotBoxAspectRatio','CameraUpVector'});

% hlinkprops = getappdata(sp1);
% addprop(hlink,'PlotBoxAspectRatio')
%% plot mean trajectories - CSD
figure; hold on
for file_idx=1:numel(meta.archive.WT_filenames)
    [X Y Z] =tubeplot(meta.WT_exported.PCAs.CSD.trajectories_mean{file_idx}(:,1),...
        meta.WT_exported.PCAs.CSD.trajectories_mean{file_idx}(:,2),...
        meta.WT_exported.PCAs.CSD.trajectories_mean{file_idx}(:,3),...
        meta.WT_exported.PCAs.CSD.trajectories_error_max{file_idx}/10);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(WT_fig,'FaceColor',cmap(file_idx,:));
    % set(WT_fig,'FaceColor','interp')
    % alpha(WT_fig,0.3)
end

for file_idx=1:numel(meta.archive.KO_filenames)
    [X Y Z] =tubeplot(meta.KO_exported.PCAs.CSD.trajectories_mean{file_idx}(:,1),...
        meta.KO_exported.PCAs.CSD.trajectories_mean{file_idx}(:,2),...
        meta.KO_exported.PCAs.CSD.trajectories_mean{file_idx}(:,3),...
        meta.KO_exported.PCAs.CSD.trajectories_error_max{file_idx}/10);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    
    % set(KO_fig,'FaceColor',cmap(file_idx,:));
    % set(KO_fig,'FaceColor','interp')
end
grid on
view(3)
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%% plot mean trajectories as subplot - MUA
figure; hold on
sp1=subplot(1,2,1);hold on
title('WT (N=15) Mean Burst Trajectories')
cmap=hsv(numel(meta.archive.WT_filenames));
for file_idx=1:numel(meta.archive.WT_filenames)
    [X Y Z] =tubeplot(meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
        meta.WT_exported.PCAs.MUA.trajectories_error_max{file_idx}/10);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(WT_fig,'FaceColor',cmap(file_idx,:));
    % set(WT_fig,'FaceColor','interp')
    % alpha(WT_fig,0.3)
end
% axis([-0.1 2 -0.6 0.6 -0.3 0.3])
axis([-2 15 -5 10 -5 10])
grid on
% view(3)
view([-23 38])
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')

sp2=subplot(1,2,2);hold on
title('KO (N=9) Mean Burst Trajectories')
cmap=hsv(numel(meta.archive.KO_filenames));
for file_idx=1:numel(meta.archive.KO_filenames)
    [X Y Z] =tubeplot(meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
        meta.KO_exported.PCAs.MUA.trajectories_error_max{file_idx}/10);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(KO_fig,'FaceColor',cmap(file_idx,:));
    % set(KO_fig,'FaceColor','interp')
end
grid on
% view(3)
view([-23 38])
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')


% axis([-0.1 2 -0.6 0.6 -0.3 0.3])
axis([-2 15 -5 10 -5 10])
% hlink = linkprop([sp1;sp2], {'CameraPosition','CameraUpVector','axis'});
