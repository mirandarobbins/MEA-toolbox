%% Indexes exported data directory and iteratively loads files for measurements
reanalyse_spikes=0;
tb=-100:899;
load pushover_parameters.mat
cname = gethostname;
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
tElapsed = tic;
no_files=numel(meta.archive.WT_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};
    try
        %     fieldname=active_file(1:numel(active_file)-4);
        %     meta.WT_exported(file_idx)=struct(fieldname,[]);
        % update progressbar
        waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...loading'] ))
        load(active_file,'-mat');
        if reanalyse_spikes==1
            waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...extracting spikes'] ))
            %             spikes = AD_spike_extract(data,params);
            spikes=meta_spike_extract(data,params,3);
            waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...processing MUA'] ))
            spikes = meta_MUA64plot(data,params,spikes);
        end
        waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...Processing PCAs'] ))
        [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
        clear active_file CSD burst data params spikes
        meta.WT_exported.PCAs.MUA.percent_var_explained{file_idx}       =   PCAs.MUA.percent_var_explained;
        meta.WT_exported.PCAs.MUA.trajectories{file_idx}                =   PCAs.MUA.trajectories;
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}           =   PCAs.MUA.trajectories_mean;
        meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV{file_idx}  =   nanmean(PCAs.MUA.trajectories_error_CV,2);
        meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}     =   nanmean(PCAs.MUA.trajectories_error,2);
        meta.WT_exported.PCAs.MUA.trajectories_mean_all_PCs{file_idx}   =   PCAs.MUA.trajectories_mean_all_PCs;
        meta.WT_exported.PCAs.MUA.trajectories_error_all_PCs{file_idx}  =   PCAs.MUA.trajectories_error_all_PCs;
        clear PCAs
        %%%%% pushover messaging
        tElapsed = toc(uint64(tElapsed));
        title = ['Execution Success on ', cname];
        message = sprintf('Execution success: %s\nTime Elapsed: %.1fs', 'WT batch analysis', tElapsed);
    catch exception
        tElapsed = toc(uint64(tElapsed));
        title = ['Execution Error on ', cname];
        message = sprintf('Error executing %s!\n---\n%s\n---\nTime Elapsed: %.1fs', 'WT batch analysis', exception.message, tElapsed);
        cprintf('red', getReport(exception));
    end
end
cprintf('blue', '\n---\nExecution Complete: %s\nTime Elapsed: %.1fs\n---\n', 'WT  batch analysis', tElapsed);
post_params = {'token', API_TOKEN,...    % API token
    'user', USER_TOKEN,...    % user's ID token
    'message', message,...    % message to push
    'title', title};          % message title in notification bar
urlread('https://api.pushover.net/1/messages.json', 'Post', post_params);
clear file_idx no_files
clear functions
clear hidden
%% KO load loop
tElapsed = tic;
no_files=numel(meta.archive.KO_filenames);
% progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
    try
        
        %     fieldname=active_file(1:numel(active_file)-4);
        %     meta.KO_exported(file_idx)=struct(fieldname,[]);
        % update progressbar
        waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...loading'] ))
        load(active_file,'-mat');
        if reanalyse_spikes==1
            waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...extracting spikes'] ))
            %             spikes = AD_spike_extract(data,params);
            spikes=meta_spike_extract(data,params,3);
            waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...processing MUA'] ))
            spikes = meta_MUA64plot(data,params,spikes);
        end
        waitbar(file_idx/no_files,progbar,strcat(['Processing file ',num2str(file_idx),'/',num2str(no_files),'...Processing PCAs'] ))
        [PCAs] = PCA_individual_trials(data,spikes,CSD,0);
        clear active_file CSD burst data params spikes
        meta.KO_exported.PCAs.MUA.percent_var_explained{file_idx}       =   PCAs.MUA.percent_var_explained;
        meta.KO_exported.PCAs.MUA.trajectories{file_idx}                =   PCAs.MUA.trajectories;
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}           =   PCAs.MUA.trajectories_mean;
        meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV{file_idx}  =   nanmean(PCAs.MUA.trajectories_error_CV,2);
        meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx}     =   nanmean(PCAs.MUA.trajectories_error,2);
        
        meta.KO_exported.PCAs.MUA.trajectories_mean_all_PCs{file_idx}   =   PCAs.MUA.trajectories_mean_all_PCs;
        meta.KO_exported.PCAs.MUA.trajectories_error_all_PCs{file_idx}  =   PCAs.MUA.trajectories_error_all_PCs;
        
        clear PCAs
        %%%%% pushover messaging
        tElapsed = toc(uint64(tElapsed));
        title = ['Execution Success on ', cname];
        message = sprintf('Execution success: %s\nTime Elapsed: %.1fs', 'KO batch analysis', tElapsed);
    catch exception
        tElapsed = toc(uint64(tElapsed));
        title = ['Execution Error on ', cname];
        message = sprintf('Error executing %s!\n---\n%s\n---\nTime Elapsed: %.1fs', 'KO batch analysis', exception.message, tElapsed);
        cprintf('red', getReport(exception));
    end
end
cprintf('blue', '\n---\nExecution Complete: %s\nTime Elapsed: %.1fs\n---\n', 'KO batch analysis', tElapsed);
post_params = {'token', API_TOKEN,...    % API token
    'user', USER_TOKEN,...    % user's ID token
    'message', message,...    % message to push
    'title', title};          % message title in notification bar
urlread('https://api.pushover.net/1/messages.json', 'Post', post_params);
clear active_file file_idx no_files
clear functions
clear hidden
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

temp_WT=cell2mat(meta.WT_exported.PCAs.MUA.percent_var_explained); temp_WT=sum(temp_WT(1:3,:));
temp_KO=cell2mat(meta.KO_exported.PCAs.MUA.percent_var_explained); temp_KO=sum(temp_KO(1:3,:));
[h,p]=ttest2(temp_WT,temp_KO) 
figure; hold on
%plot absolute
% plot(meta.WT_exported.PCAs.MUA.percent_var_explained_mean,'b','LineWidth',1.5)
% plot(meta.KO_exported.PCAs.MUA.percent_var_explained_mean,'r','LineWidth',1.5)
temp=horzcat(meta.WT_exported.PCAs.MUA.percent_var_explained_mean,meta.KO_exported.PCAs.MUA.percent_var_explained_mean);
temp_map=[0.3 0.3 1;1 0.3 0.3];
bar_h=bar(temp(1:10,:),1.5);
set(bar_h,'linewidth',0.001,'edgecolor','none'); colormap(temp_map)
errorbar((1:10)-0.14,meta.WT_exported.PCAs.MUA.percent_var_explained_mean(1:10),...
        meta.WT_exported.PCAs.MUA.percent_var_explained_SEM(1:10),'b','LineStyle','none','LineWidth',2)
errorbar((1:10)+0.14,meta.KO_exported.PCAs.MUA.percent_var_explained_mean(1:10),...
        meta.KO_exported.PCAs.MUA.percent_var_explained_SEM(1:10),'r','LineStyle','none','LineWidth',2)

% ciplot(meta.WT_exported.PCAs.MUA.percent_var_explained_mean+meta.WT_exported.PCAs.MUA.percent_var_explained_SEM,...
%        meta.WT_exported.PCAs.MUA.percent_var_explained_mean-meta.WT_exported.PCAs.MUA.percent_var_explained_SEM,1:64,'b')
% ciplot(meta.KO_exported.PCAs.MUA.percent_var_explained_mean+meta.KO_exported.PCAs.MUA.percent_var_explained_SEM,...
%        meta.KO_exported.PCAs.MUA.percent_var_explained_mean-meta.KO_exported.PCAs.MUA.percent_var_explained_SEM,1:64,'r')



%plot cumulative
ciplot(meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean+meta.WT_exported.PCAs.MUA.percent_var_explained_cum_SEM,...
    meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean-meta.WT_exported.PCAs.MUA.percent_var_explained_cum_SEM,1:64,'b')
plot(meta.WT_exported.PCAs.MUA.percent_var_explained_cum_mean,'b','LineWidth',1.5);

ciplot(meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean+meta.KO_exported.PCAs.MUA.percent_var_explained_cum_SEM,...
    meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean-meta.KO_exported.PCAs.MUA.percent_var_explained_cum_SEM,1:64,'r')
plot(meta.KO_exported.PCAs.MUA.percent_var_explained_cum_mean,'r','LineWidth',1.5);
% text(7, 10, 'MUA','FontSize', 18)
% axis([0.7 5.3 0 100])
axis([0.7 10.3 0 100])
xlabel('Principal Component no.')
ylabel({'% Variance explained'})
%% plot mean trajectories - MUA
figure; hold on
for file_idx=1:numel(meta.archive.WT_filenames)
    [X Y Z] =tubeplot(meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
        meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}/10);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(WT_fig,'FaceColor','interp')
    alpha(WT_fig,0.5)
end

for file_idx=1:numel(meta.archive.KO_filenames)
    [X Y Z] =tubeplot(meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
        meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx}/10);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    % set(KO_fig,'FaceColor','interp')
    % alpha(KO_fig,0.5)
end
grid off
% box on
view(3)
axis square
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
axis([-0 2 -Inf Inf -Inf Inf])
%% extract pop data for time-resolved variability

meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean          =  nanmean(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean),2);
meta.WT_exported.PCAs.MUA.trajectories_error_grand_SEM           =  nansem (cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean),2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean  =  nanmean(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV),2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM   =  nansem (cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV),2);

meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean          =  nanmean(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean),2);
meta.KO_exported.PCAs.MUA.trajectories_error_grand_SEM           =  nansem (cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean),2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean  =  nanmean(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV),2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM   =  nansem (cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV),2);

temp_WT=cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean); temp_WT(temp_WT==0)=NaN;   temp_WT=nanmean(temp_WT(260:300,:));
temp_KO=cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean); temp_KO(temp_KO==0)=NaN;   temp_KO=nanmean(temp_KO(260:300,:));
[h,p]=ttest2(temp_WT,temp_KO) ;
temp_WT=nansem(temp_WT);temp_KO=nansem(temp_KO);

% Plot mean±SEM trajectory error
figure; hold on 
    plot(tb,cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean),'b','LineWidth',1)
    plot(tb,cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean),'r','LineWidth',1)
    plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean,'b','LineWidth',1.5)
    plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean,'r','LineWidth',1.5)
    ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_grand_SEM,...
        meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_grand_SEM,...
        tb,'b')
    ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_grand_SEM,...
        meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_grand_SEM,...
        tb,'r')
    %    plot([100 300],[0,0],'k','linewidth',1.5)
%     axis([100 300 0 0.3])
    xlabel('Post-stimulus time (ms)')
    ylabel('Trial-to-trial variability (A.U.)')

% % Plot mean trajectory error as solid bar
% figure; hold on
%     ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean,...
%         -1*meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean,...
%         tb,'b')
%     ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean,...
%         -1*meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean,...
%         tb,'r')
%     %    plot([100 300],[0,0],'k','linewidth',1.5)
%     axis([100 300 0 0.15])
%     xlabel('Post-stimulus time (ms)')
%     ylabel('Trial-to-trial variability (A.U.)')
%% extract pop data for time-resolved variability cumulative
top_PCs =3; % no top PCs to include in error estimate
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative=[];                    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm=[];               meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs=[];            meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean=[];       meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm=[];       meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean=[];  meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean=[]; meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean=[];

for file_idx=1:numel(meta.archive.WT_filenames)
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}                   = cumsum(meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx});
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx}              = meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}...
                                                                                                    ./max(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})*100;
    
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}           = cumsum(meta.WT_exported.PCAs.MUA.trajectories_error_all_PCs{file_idx}(:,1:top_PCs));
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean(:,file_idx)    = nanmean(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx},2);

        
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm{file_idx}      = meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}(:,1:top_PCs)...
                                                                                                    ./repmat(max(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}(:,1:top_PCs)),...
                                                                                                        length(meta.WT_exported.PCAs.MUA.trajectories_error_grand_mean),1)*100;
    
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean(:,file_idx)  = nanmean(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm{file_idx},2);
             
end
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean=nanmean(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean,2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_sem=nansem(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean,2);

meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_grand_mean=nanmean(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean,2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_grand_sem=nansem(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean,2);

for file_idx=1:numel(meta.archive.KO_filenames)
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}                   = cumsum(meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx});
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx}              = meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}...
                                                                                                    ./max(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})*100;
    
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}           = cumsum(meta.KO_exported.PCAs.MUA.trajectories_error_all_PCs{file_idx});
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean(:,file_idx)    = nanmean(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}(:,1:top_PCs),2);

        
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm{file_idx}      = meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}(:,1:top_PCs)...
                                                                                                    ./repmat(max(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs{file_idx}(:,1:top_PCs)),...
                                                                                                        length(meta.KO_exported.PCAs.MUA.trajectories_error_grand_mean),1)*100;
    
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean(:,file_idx)  = nanmean(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm{file_idx},2);
             
end
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean=nanmean(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean,2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_sem=nansem(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean,2);

meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_grand_mean=nanmean(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean,2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_grand_sem=nansem(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_norm_mean,2);

figure; hold on
plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean,'b','LineWidth',1.5)
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean,'r','LineWidth',1.5)
ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_sem,...
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_sem,...
    tb,'b');

ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_sem,...
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_grand_sem,...
    tb,'r');

plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean,'b')
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean,'r')
temp_WT=meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean(300,:);
temp_KO=meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_all_PCs_mean(300,:);
temp_WT(temp_WT==0)=[];temp_KO(temp_KO==0)=[];
[h,p]=ttest2(temp_WT,temp_KO)

%%
temp=cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative);
temp(:,std(temp)==0)=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean      = nanmean(temp,2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM       = nansem(temp,2);

temp=cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm);
temp(:,std(temp)==0)=[];
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean = nanmean(temp,2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_SEM  = nansem(temp,2);

for file_idx=1:numel(meta.archive.KO_filenames)
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}=cumsum(meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx});
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}=meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx};
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx}=meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx}...
        ./max(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})*100;
end
temp=cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative);
temp(:,std(temp)==0)=[];
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean      = nanmean(temp,2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM       = nansem(temp,2);
temp=cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm);
temp(:,std(temp)==0)=[];

meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean = nanmean(temp,2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_SEM  = nansem(temp,2);

temp_WT=cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative); 
temp_KO=cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative); 

figure; hold on;plot(temp_WT,'b');plot(temp_KO,'r')
temp_WT=temp_WT(250,:);temp_WT(temp_WT==0)=NaN;  
temp_KO=temp_KO(250,:);temp_KO(temp_KO==0)=NaN;   
[h,p]=ttest2(temp_WT,temp_KO) 
% temp_WT=nansem(temp_WT);temp_KO=nansem(temp_KO);

%%%%%
% (1) plot mean accumulation of variance across top 3 PC's
figure; hold on
%     for file_idx=1:numel(meta.archive.WT_filenames)
%         plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx},'b')%./max(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})
%     end
%     for file_idx=1:numel(meta.archive.KO_filenames)
%         plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx},'r')%./max(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative{file_idx})
%     end
ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    tb,'b')
plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean,'b','LineWidth',1.5)
ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_SEM,...
    tb,'r')
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_grand_mean,'r','LineWidth',1.5)
xlabel('Post-stimulus time (ms)')
ylabel({'Cumulative trial-trial';'variability (A.U.)'})
axis([100 300 0 50])

%%%%%
% (2) plot mean accumulation of variance across top 3 PC's
figure; hold on
%     for file_idx=1:numel(meta.archive.WT_filenames)
%         plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx},'b')%./max(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx})
%     end
%     for file_idx=1:numel(meta.archive.KO_filenames)
%         plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx},'r')%./max(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm{file_idx})
%     end
ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_SEM,...
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_SEM,...
    tb,'b')
plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean,'b','LineWidth',1.5)
ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_SEM,...
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_SEM,...
    tb,'r')
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_cumulative_norm_grand_mean,'r','LineWidth',1.5)
xlabel('Post stimulus time (ms)')
ylabel('Cumulative variance (A.U.)')
axis([100 300 0 100])

%% inspect CV vs time
meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean = nanmean(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV),2);
meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM = nansem(cell2mat(meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV),2);

meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean = nanmean(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV),2);
meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM = nansem(cell2mat(meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV),2);
figure; hold on
% for file_idx=1:numel(meta.archive.WT_filenames)
%     plot(meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV{file_idx},'b');
% end
% for file_idx=1:numel(meta.archive.KO_filenames)
%     plot(meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV{file_idx},'r');
% end
plot(tb,meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean,'b','LineWidth',1.5)
plot(tb,meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean,'r','LineWidth',1.5)
ciplot(meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean-meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM,...
    meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean+meta.WT_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM,...
    tb,'b')

ciplot(meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean-meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM,...
    meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_mean+meta.KO_exported.PCAs.MUA.trajectories_error_mean_CV_grand_SEM,...
    tb,'r')
%% wrap up

notify('I`m finished, come and take a look!!!')