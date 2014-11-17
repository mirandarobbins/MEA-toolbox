%% Plot to highlight WT
plot_window_wide=[200 1000];
plot_window=[278 280];
highlight_time=279%plot_window(2);
divide_factor=2 % scale factor for mean± error plots

plot_width_all=0.02;
plot_width_highlight=plot_width_all*2;
%%% WT plot
WT_fig_plot=figure;
cameratoolbar('Show');
cameratoolbar('SetMode','orbit','SetCoordSys','none');
set(WT_fig_plot,'pos',[11 21 1268 958],'color','w')
% subplot(1,2,1)
hold on
for file_idx=9%:numel(meta.WT_exported.PCAs.MUA.trajectories)
    %     subplot(round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),round(numel(meta.WT_exported.PCAs.MUA.trajectories)^0.5),file_idx)
    hold on
    cmap=hsv(numel(meta.WT_exported.PCAs.MUA.trajectories));
    no_trials=numel(meta.WT_exported.PCAs.MUA.trajectories{file_idx});
    clear X Y Z X_all Y_all Z_all plot_targets
    
    %%% find coordinates for straight line plots
    for trial_id=1:no_trials
        plot_targets(trial_id,1)=meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(highlight_time,1)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,1);
        plot_targets(trial_id,2)=meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(highlight_time,2)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2);
        plot_targets(trial_id,3)=meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(highlight_time,3)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,3);
    end
    plot_targets(no_trials+1,1)=meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(highlight_time,1)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,1);
    plot_targets(no_trials+1,2)=meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(highlight_time,2)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,2);
    plot_targets(no_trials+1,3)=meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(highlight_time,3)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,3);
    
    %%%Individual trials - highlight
    for trial_id=1:no_trials
        [X{trial_id} ...
            Y{trial_id} ...
            Z{trial_id}] = tubeplot(...
            meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window(1):plot_window(2),1)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,1),...
            meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window(1):plot_window(2),2)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window(1):plot_window(2),3)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            plot_width_highlight);
    end
    %%% MEAN over trials - highlight
    [X{no_trials+1} ...
        Y{no_trials+1} ...
        Z{no_trials+1}] = tubeplot(...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window(1):plot_window(2),1)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,1),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window(1):plot_window(2),2)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,2),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window(1):plot_window(2),3)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,3),...
        0.02)%meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}(plot_window(1):plot_window(2))/5);
    
    %%%Individual trials - all
    for trial_id=1:no_trials
        [X_all{trial_id} ...
            Y_all{trial_id} ...
            Z_all{trial_id}] = tubeplot(...
            meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window_wide(1):plot_window_wide(2),1)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,1),...
            meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window_wide(1):plot_window_wide(2),2)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window_wide(1):plot_window_wide(2),3)-meta.WT_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            plot_width_all);
    end
    %%% MEAN over trials - all
    [X_all{no_trials+1} ...
        Y_all{no_trials+1} ...
        Z_all{no_trials+1}] = tubeplot(...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window_wide(1):plot_window_wide(2),1)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,1),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window_wide(1):plot_window_wide(2),2)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,2),...
        meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window_wide(1):plot_window_wide(2),3)-meta.WT_exported.PCAs.MUA.trajectories_mean{file_idx}(1,3),...
        meta.WT_exported.PCAs.MUA.trajectories_error_mean{file_idx}(plot_window_wide(1):plot_window_wide(2))/divide_factor);%0.018) 
    
        
    %%%% Plot window wide %%%%
    for trial_id=1:no_trials
        fig_X_all=surf(X_all{trial_id},Y_all{trial_id},Z_all{trial_id});
        set(fig_X_all,'EdgeColor','none','FaceColor',[0.9 0.9 1],'FaceLighting','phong')
        alpha(fig_X_all,0.1)
    end
    fig_X_all_mean=surf(X_all{no_trials+1},Y_all{no_trials+1},Z_all{no_trials+1});
    % set(fig_X_all_mean,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(fig_X_all_mean,'EdgeColor','none','FaceColor',[0.2 0.2 1],'FaceLighting','phong')
    alpha(fig_X_all_mean,0.4)
    
    %%%% Plot highlight %%%%
    for trial_id=1:no_trials
        fig_X_highlight=surf(X{trial_id},Y{trial_id},Z{trial_id});
        set(fig_X_highlight,'EdgeColor','none','FaceColor',[0.2 0.2 1],'FaceLighting','phong')
        alpha(fig_X_highlight,0.8)
    end
    fig_X_highlight_mean=surf(X{no_trials+1},Y{no_trials+1},Z{no_trials+1});
    % set(fig_X_highlight_mean,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(fig_X_highlight_mean,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
    alpha(fig_X_highlight_mean,1)
    
    %%%% Dotted coordinate lines %%%%
    for trial_id=1:no_trials
        plot3([plot_targets(trial_id,1) plot_targets(no_trials+1,1)],...
            [plot_targets(trial_id,2) plot_targets(no_trials+1,2)],...
            [plot_targets(trial_id,3) plot_targets(no_trials+1,3)],...
            '-b','LineWidth',2)
    end
%     axis([-0.5 6 -1.5 2.5 -1 2.5])
    axis([0.5 5 -1.8 1.5 -1 1])

    grid on
    view(-37,58)
    axis vis3d
%     axis off
    camlight
    xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%     zoom(1.5)
end

% plot3([-1. -1.],[0 2],[0 0],'k','LineWidth',4)

% plot3([-1. -1.],[0 0],[0 1],'k','LineWidth',4)
%% Plot to highlight KO
plot_window_wide=[200 400];
plot_window=[250 252];
highlight_time=251%plot_window(2);
divide_factor=2 % scale factor for mean± error plots

plot_width_all=0.02;
plot_width_highlight=plot_width_all*1.5;
%%% KO plot
KO_fig_plot=figure;
cameratoolbar('Show');
cameratoolbar('SetMode','orbit','SetCoordSys','none');
set(KO_fig_plot,'pos',[11 21 1268 958],'color','w')
% subplot(1,2,2)
hold on
for file_idx=5%:numel(meta.KO_exported.PCAs.MUA.trajectories)
    %     subplot(round(numel(meta.KO_exported.PCAs.MUA.trajectories)^0.5),round(numel(meta.KO_exported.PCAs.MUA.trajectories)^0.5),file_idx)
    hold on
    cmap=hsv(numel(meta.KO_exported.PCAs.MUA.trajectories));
    no_trials=numel(meta.KO_exported.PCAs.MUA.trajectories{file_idx});
    clear X Y Z X_all Y_all Z_all plot_targets
    
    %%% find coordinates for straight line plots
    for trial_id=1:no_trials
        plot_targets(trial_id,1)=meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(highlight_time,1)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,1);
        plot_targets(trial_id,2)=meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(highlight_time,2)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2);
        plot_targets(trial_id,3)=meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(highlight_time,3)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,3);
    end
    plot_targets(no_trials+1,1)=meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(highlight_time,1)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,1);
    plot_targets(no_trials+1,2)=meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(highlight_time,2)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,2);
    plot_targets(no_trials+1,3)=meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(highlight_time,3)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,3);
    
    %%%Individual trials - highlight
    for trial_id=1:no_trials
        [X{trial_id} ...
            Y{trial_id} ...
            Z{trial_id}] = tubeplot(...
            meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window(1):plot_window(2),1)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,1),...
            meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window(1):plot_window(2),2)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window(1):plot_window(2),3)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            plot_width_highlight);
    end
    %%% MEAN over trials - highlight
    [X{no_trials+1} ...
        Y{no_trials+1} ...
        Z{no_trials+1}] = tubeplot(...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window(1):plot_window(2),1)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,1),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window(1):plot_window(2),2)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,2),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window(1):plot_window(2),3)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,3),...
        0.0005)%meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx}(plot_window(1):plot_window(2))/5);
    
    %%%Individual trials - all
    for trial_id=1:no_trials
        [X_all{trial_id} ...
            Y_all{trial_id} ...
            Z_all{trial_id}] = tubeplot(...
            meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window_wide(1):plot_window_wide(2),1)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,1),...
            meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window_wide(1):plot_window_wide(2),2)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(plot_window_wide(1):plot_window_wide(2),3)-meta.KO_exported.PCAs.MUA.trajectories{file_idx}{trial_id}(1,2),...
            plot_width_all);
    end
    %%% MEAN over trials - all
    [X_all{no_trials+1} ...
        Y_all{no_trials+1} ...
        Z_all{no_trials+1}] = tubeplot(...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window_wide(1):plot_window_wide(2),1)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,1),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window_wide(1):plot_window_wide(2),2)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,2),...
        meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(plot_window_wide(1):plot_window_wide(2),3)-meta.KO_exported.PCAs.MUA.trajectories_mean{file_idx}(1,3),...
        meta.KO_exported.PCAs.MUA.trajectories_error_mean{file_idx}(plot_window_wide(1):plot_window_wide(2))/divide_factor);%0.018) 

    
        
    %%%% Plot window wide %%%%
    for trial_id=1:no_trials
        fig_X_all=surf(X_all{trial_id},Y_all{trial_id},Z_all{trial_id});
        set(fig_X_all,'EdgeColor','none','FaceColor',[0.9 0.9 1],'FaceLighting','phong')
        alpha(fig_X_all,0.1)
    end
    fig_X_all_mean=surf(X_all{no_trials+1},Y_all{no_trials+1},Z_all{no_trials+1});
    % set(fig_X_all_mean,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(fig_X_all_mean,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')
    alpha(fig_X_all_mean,0.4)
    
    %%%% Plot highlight %%%%
    for trial_id=1:no_trials
        fig_X_highlight=surf(X{trial_id},Y{trial_id},Z{trial_id});
        set(fig_X_highlight,'EdgeColor','none','FaceColor',[1 0.2 0.2],'FaceLighting','phong')
        alpha(fig_X_highlight,0.8)
    end
    fig_X_highlight_mean=surf(X{no_trials+1},Y{no_trials+1},Z{no_trials+1});
    % set(fig_X_highlight_mean,'EdgeColor','none','FaceColor',cmap(file_idx,:),'FaceLighting','phong')
    set(fig_X_highlight_mean,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')
    alpha(fig_X_highlight_mean,1)
    
    %%%% Dotted coordinate lines %%%%
    for trial_id=1:no_trials
        plot3([plot_targets(trial_id,1) plot_targets(no_trials+1,1)],...
            [plot_targets(trial_id,2) plot_targets(no_trials+1,2)],...
            [plot_targets(trial_id,3) plot_targets(no_trials+1,3)],...
            'r','LineWidth',2)
    end
    axis([0.5 5 -1.8 1.5 -1 1])
    grid on
    view(-37,58)
    box on
    axis vis3d
    axis off
    camlight
    xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%     zoom(1.5)
end
