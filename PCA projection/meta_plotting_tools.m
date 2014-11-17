%% colour-coded by time
    plot_until=[200 250];

    
    fig_plot=figure; 
cameratoolbar('Show');
    cameratoolbar('SetMode','orbit','SetCoordSys','none');
set(fig_plot,'pos',[11 21 1268 958],'color','w')
hold on
time_rainbow=jet(numel(meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1)));    
for file_idx=3%:numel(meta.WT_exported.PCAs.MUA.trajectories)
%     subplot(round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),file_idx)
            hold on
            cmap=hsv(numel(meta.WT_exported.PCAs.MUA.trajectories));
    no_trials=numel(meta.WT_exported.PCAs.MUA.trajectories{file_idx});
    clear X Y Z
%%%Individual trials
    for trial_id=1:3%no_trials
        [X{trial_id} ...
         Y{trial_id} ...
         Z{trial_id}] = tubeplot(...
                       meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_until(1):plot_until(2),1),...
                       meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_until(1):plot_until(2),2),...
                       meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_until(1):plot_until(2),3),...
                       0.01);
    end
%%% MEAN over trials 
    [X{no_trials+1} ...
     Y{no_trials+1} ...
     Z{no_trials+1}] = tubeplot(...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),1),...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),2),...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),3),...
           meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}(plot_until(1):plot_until(2))/5);
    for trial_id=1:no_trials
        WT_fig=surf(X{trial_id},Y{trial_id},Z{trial_id});
        set(WT_fig,'EdgeColor','none','FaceColor',[0.8 0.8 1],'FaceLighting','phong')
%         alpha(WT_fig,0.6) 
    end
     WT_fig=surf(X{no_trials+1},Y{no_trials+1},Z{no_trials+1});
    % set(WT_fig,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
    alpha(WT_fig,1) 
axis([-0.5 2 -1 1.5 -1 1.5])
grid on
view(1,30)
axis vis3d
axis off
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
zoom(1.5)
end

%%
figure;
for file_idx=1:numel(meta.KO_exported.PCAs.MUA.trajectories)
    subplot(round(numel(meta.KO_exported.PCAs.MUA.trajectories)^0.5),round(numel(meta.KO_exported.PCAs.MUA.trajectories)^0.5),file_idx)
            hold on
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

    
axis([-0.5 2 -1 1.5 -1 1.5])
grid on
view(1,30)
axis vis3d
axis off
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
zoom(1.5)

    
end

fig_plot=figure; 
set(fig_plot,'pos',[11 21 1268 958],'color','w')
hold on

% sp1=subplot(1,2,1);hold on
% title('WT & KO')
    
for file_idx=1:numel(meta.WT_exported.PCAs.MUA.trajectories)
    subplot(round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),file_idx)
            hold on
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
axis([-0.5 2 -1 1.5 -1 1.5])
grid on
view(1,30)
axis vis3d
axis off
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
zoom(1.5)
end


%% coloured by time

fig_plot=figure; 
cameratoolbar('Show');
    cameratoolbar('SetMode','orbit','SetCoordSys','none');
set(fig_plot,'pos',[11 21 1268 958],'color','w')
hold on
time_rainbow=jet(numel(meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(:,1)));    
for file_idx=2%:numel(meta.WT_exported.PCAs.MUA.trajectories)
%     subplot(round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),file_idx)
            hold on
%             cmap=hsv(numel(meta.WT_exported.PCAs.MUA.trajectories));
    no_trials=numel(meta.WT_exported.PCAs.MUA.trajectories{file_idx});
    clear X Y Z
    
%%%Individual trials
    for trial_id=1:no_trials
         X{trial_id} = meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_until(1):plot_until(2),1);
         Y{trial_id} = meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_until(1):plot_until(2),2);
         Z{trial_id} = meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_until(1):plot_until(2),3);
    end
% MEAN over trials 
         X{no_trials+1} =  meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),1);
         Y{no_trials+1} =  meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),2);
         Z{no_trials+1} =  meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),3);
    

    for trial_id=1:no_trials
        WT_fig=clinep(X{trial_id},...
                      Y{trial_id},...
                      Z{trial_id},...
                      0:numel(X{trial_id})-1,3); %
    end
    [X{no_trials+1} ...
     Y{no_trials+1} ...
     Z{no_trials+1}] = tubeplot(...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),1),...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),2),...
           meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_until(1):plot_until(2),3),...
           meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}(plot_until(1):plot_until(2))/2);
    
     WT_fig=surf(X{no_trials+1},Y{no_trials+1},Z{no_trials+1});
    set(WT_fig,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
    alpha(WT_fig,1) 
axis([-0.5 2 -1 1.5 -1 1.5])
caxis(plot_until)%numel(X{trial_id})])
colormap jet 
grid on
view(1,30)
axis vis3d
axis off
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
zoom(1.5)
end
