%% prepare mean trajectories - MUA each cell mean only
% one colour per recording
figure; hold on
sp1=subplot(1,2,1);hold on
title('WT (N=15) Mean Burst Trajectories')
cmap=hsv(numel(meta.archive.WT_filenames));
clear X Y Z
for file_idx=1:numel(meta.archive.WT_filenames)
 [X{file_idx} Y{file_idx} Z{file_idx}] =tubeplot(meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
                                                 meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
                                                 meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
                                                 0.005);%meta.WT_exported.PCAs.MUA.trajectories_error_max{file_idx}/10
end
for file_idx=1:numel(meta.archive.WT_filenames)
    WT_fig=surf(X{file_idx},Y{file_idx},Z{file_idx});
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    set(WT_fig,'FaceColor',cmap(file_idx,:));
% set(WT_fig,'FaceColor','interp')
% alpha(WT_fig,0.3) 
end
axis([-0.1 2 -0.6 0.6 -0.3 0.3])
grid on
% view(3)
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
title('KO (N=15) Mean Burst Trajectories')
cmap=hsv(numel(meta.archive.KO_filenames));
clear X Y Z
sp1=subplot(1,2,2);hold on

for file_idx=1:numel(meta.archive.KO_filenames)
 [X{file_idx} Y{file_idx} Z{file_idx}] =tubeplot(meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
                                                 meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
                                                 meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
                                                 0.005);%meta.KO_exported.PCAs.MUA.trajectories_error_max{file_idx}/10
end
for file_idx=1:numel(meta.archive.KO_filenames)
    KO_fig=surf(X{file_idx},Y{file_idx},Z{file_idx});
    set(KO_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
    set(KO_fig,'FaceColor',cmap(file_idx,:));
% set(KO_fig,'FaceColor','interp')
% alpha(KO_fig,0.3) 
end
axis([-0.1 2 -0.6 0.6 -0.3 0.3])
grid on
% view(3)
axis vis3d
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%% plot moving means - movie
M = sweep2movie_PCA3D_multitrace(X,Y,Z,[1 500],'axlims','auto');
last_dir=pwd;
if  ~isdir('movies')
    mkdir('movies')
end
target_dir=strcat(last_dir,'/movies');
cd movies
movie2avi(M, ['PCA_average.avi'], 'compression', 'None')
cd(last_dir)
close(gcf)


clear M movie trial_id input last_dir target_dir
disp( 'finished coding AVI files!')
%% prepare mean trajectories as subplot - MUA one cell with means
% One giant spaghetti-like mess
fig_plot=figure; 
set(fig_plot,'pos',[11 21 1268 958],'color','w')
hold on

% sp1=subplot(1,2,1);hold on
title('WT & KO')

for file_idx=1:numel(meta.WT_exported.PCAs.MUA.trajectories)
    cmap=hsv(numel(meta.WT_exported.PCAs.MUA.trajectories));
    no_trials=numel(meta.WT_exported.PCAs.MUA.trajectories{file_idx});
    clear X Y Z
%%%Individual trials
    for trial_id=1:no_trials
        [X{trial_id} ...
         Y{trial_id} ...
         Z{trial_id}] = tubeplot(...
                       meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(:,1),...
                       meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(:,2),...
                       meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(:,3),...
                       0.01);
    end
%%% MEAN over trials 
    [X{no_trials+1} ...
     Y{no_trials+1} ...
     Z{no_trials+1}] = tubeplot(...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
           meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}/20);
    for trial_id=1:no_trials
        WT_fig=surf(X{trial_id},Y{trial_id},Z{trial_id});
        set(WT_fig,'EdgeColor','none','FaceColor',[0.8 0.8 1],'FaceLighting','phong')
%         alpha(WT_fig,0.6) 
    end
     WT_fig=surf(X{no_trials+1},Y{no_trials+1},Z{no_trials+1});
    % set(WT_fig,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
    alpha(WT_fig,1) 

end

for file_idx=1:numel(meta.KO_exported.PCAs.MUA.trajectories)
    cmap=hsv(numel(meta.KO_exported.PCAs.MUA.trajectories));
    no_trials=numel(meta.KO_exported.PCAs.MUA.trajectories{file_idx});
    clear X Y Z
%%%Individual trials
    for trial_id=1:no_trials
        [X{trial_id} ...
         Y{trial_id} ...
         Z{trial_id}] = tubeplot(...
           meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(:,1),...
           meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(:,2),...
           meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(:,3),...
           0.01);
    end
%%%% MEAN over trials 
    [X{no_trials+1} ...
     Y{no_trials+1} ...
     Z{no_trials+1}] = tubeplot(...
           meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1),...
           meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,2),...
           meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(:,3),...
           meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx}/20);
    for trial_id=1:no_trials
        KO_fig=surf(X{trial_id},Y{trial_id},Z{trial_id});
        set(KO_fig,'EdgeColor','none','FaceColor',[1 0.8 0.8],'FaceLighting','phong')
        cmap(file_idx,:)% set(KO_fig,'FaceColor','interp')
%         alpha(KO_fig,0.6) 
    end
    KO_fig=surf(X{no_trials+1},Y{no_trials+1},Z{no_trials+1});
%     set(KO_fig,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')
    alpha(KO_fig,1) 

end

% axis([-0.1 1.5 -0.6 0.5 -0.4 0.3])
grid on
view(1,30)
axis vis3d
axis off
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
zoom(1.5)
set(fig_plot,'NextPlot','replacechildren');
%% Rotate and make movie

for frame_id = 1:180
    view(frame_id,30)
    set(fig_plot,'pos',[11 21 1268 958],'color','w')

    M(frame_id) = getframe(fig_plot);

end

last_dir=pwd;
if  ~isdir('movies')
    mkdir('movies')
end
target_dir=strcat(last_dir,'/movies');
cd movies
movie2avi(M, ['PCA_all_rotate.avi'], 'compression', 'None')
cd(last_dir)
close(gcf)


clear M movie trial_id input last_dir target_dir
disp( 'finished coding AVI files!')

