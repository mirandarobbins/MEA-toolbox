function [PCAs] = PCA_individual_trials_MUAonly(data_in,plotYN)
% if nargin<4    
% data_opt=2; %1=MUA, 2=LFP, 3=CSD
% end
for data_opt=1%:3;
    % time_lims=100:200;
    time_lims=1:1000;
    no_trials=size(data_in{1},1);
    n_time_lims=numel(time_lims);
    clear temp temp2 COEFF SCORE
    for trial_id= 1:no_trials
        for channel_id=1:64
            switch data_opt
                case 1
                    temp{trial_id}(channel_id,:)=data_in{channel_id}(trial_id,time_lims);
    %                 temp{trial_id}(channel_id,:)=temp{trial_id}(channel_id,:)./max(temp{trial_id}(channel_id,:));
                case 2
                    temp2=downsample(data.filtered_lfp(:,channel_id,trial_id),20);
                    temp{trial_id}(channel_id,:)= sgolayfilt(temp2(time_lims),3,15);
                case 3 
                    temp2=downsample(CSD.csd_array(:,channel_id,trial_id),20);
                    temp{trial_id}(channel_id,:)= sgolayfilt(temp2(time_lims),3,15);
            end
        end
    end

    temp=cell2mat(temp);temp(isnan(temp))=0;
    % temp=temp./max(max(temp));
    input=temp'; 
    clear temp temp2 trial_id channel_id
%% PCA overlaid
    clear PC1 PC1_mean PC1_SEM PC2 PC2_mean PC2_SEM PC3 PC3_mean PC3_SEM
    [COEFF,SCORE,var_explained]=princomp(input);
    percent_var_explained = 100*var_explained/sum(var_explained);
    for trial_id=1:no_trials
        clear X Y Z %PC1 PC2 PC3
        PC1(:,trial_id)=SCORE(1:n_time_lims,1);
        PC2(:,trial_id)=SCORE(1:n_time_lims,2);
        PC3(:,trial_id)=SCORE(1:n_time_lims,3);
        SCORE(1:n_time_lims,:)=[];
    end
    PC1_mean=nanmean(PC1,2);PC1_SEM=nansem(PC1,2);
    PC2_mean=nanmean(PC2,2);PC2_SEM=nansem(PC2,2);
    PC3_mean=nanmean(PC3,2);PC3_SEM=nansem(PC3,2);
    error=mean(horzcat(PC1_SEM,PC2_SEM,PC3_SEM),2);
    PC_temp=cell(1,no_trials);
    for trial_id=1:no_trials
        PC_temp{trial_id}=horzcat(PC1(:,trial_id),PC2(:,trial_id),PC3(:,trial_id));
    end
    switch data_opt
        case 1
              PCAs.MUA.trajectories=PC_temp;
              PCAs.MUA.trajectories_mean=horzcat(PC1_mean,PC2_mean,PC3_mean);
              PCAs.MUA.trajectories_error=horzcat(PC1_SEM,PC2_SEM,PC3_SEM);
              PCAs.MUA.percent_var_explained=percent_var_explained;
        case 2
              PCAs.LFP.trajectories=PC_temp;
              PCAs.LFP.trajectories_mean=horzcat(PC1_mean,PC2_mean,PC3_mean);
              PCAs.LFP.trajectories_error=horzcat(PC1_SEM,PC2_SEM,PC3_SEM);
              PCAs.LFP.percent_var_explained=percent_var_explained;
        case 3
              PCAs.CSD.trajectories=PC_temp;
              PCAs.CSD.trajectories_mean=horzcat(PC1_mean,PC2_mean,PC3_mean);
              PCAs.CSD.trajectories_error=horzcat(PC1_SEM,PC2_SEM,PC3_SEM);
              PCAs.CSD.percent_var_explained=percent_var_explained;
    end
    clear PC1 PC1_mean PC1_SEM PC2 PC2_mean PC2_SEM PC3 PC3_mean PC3_SEM X Y Z error COEFF SCORE WT_fig channel_id cmap PC_temp
    clear input n_time_lims temp temp2 time_lims trial_id percent_var_explained var_explained
end
%% optional plotting
if plotYN==1
    figure;
    set(gcf,'pos',[10 300 1200 600])
    sp1=subplot(1,2,1);hold on
    cmap=hsv(no_trials);
    for trial_id=1:no_trials
        switch data_opt
            case 1
                 [X Y Z] =tubeplot(PCAs.MUA.trajectories{trial_id}(:,1),...
                                   PCAs.MUA.trajectories{trial_id}(:,2),...
                                   PCAs.MUA.trajectories{trial_id}(:,3),0.01);
            case 2
                 [X Y Z] =tubeplot(PCAs.LFP.trajectories{trial_id}(:,1),...
                                   PCAs.LFP.trajectories{trial_id}(:,2),...
                                   PCAs.LFP.trajectories{trial_id}(:,3),0.002);
            case 3 
                 [X Y Z] =tubeplot(PCAs.CSD.trajectories{trial_id}(:,1),...
                                   PCAs.CSD.trajectories{trial_id}(:,2),...
                                   PCAs.CSD.trajectories{trial_id}(:,3),2);
        end
        WT_fig=surf(X,Y,Z);
            set(WT_fig,'EdgeColor','none','FaceColor',cmap(trial_id,:),'FaceLighting','phong')
%           set(WT_fig,'EdgeColor','none','FaceColor',[0.9 0.9 0.9],'FaceLighting','phong')
%           set(WT_fig,'FaceColor','interp')
            grid on
            view(3)
            axis vis3d
            camlight
            xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
    end
    sp2=subplot(1,2,2);hold on
    data_opt=1;
    switch data_opt
            case 1 
                [X Y Z] =tubeplot(PCAs.MUA.trajectories_mean(:,1),...
                                  PCAs.MUA.trajectories_mean(:,2),...
                                  PCAs.MUA.trajectories_mean(:,3),...
                                  PCAs.MUA.trajectories_error(:,1));
                disp('plotting MUA')
            case 2 
                [X Y Z] =tubeplot(PCAs.LFP.trajectories_mean(:,1),...
                                  PCAs.LFP.trajectories_mean(:,2),...
                                  PCAs.LFP.trajectories_mean(:,3),...
                                  PCAs.LFP.trajectories_error(:,1));
                disp('plotting LFP')
            case 3 
                [X Y Z] =tubeplot(PCAs.CSD.trajectories_mean(:,1),...
                                  PCAs.CSD.trajectories_mean(:,2),...
                                  PCAs.CSD.trajectories_mean(:,3),...
                                  PCAs.CSD.trajectories_error(:,1));
                disp('plotting CSD')
    end
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
%   set(WT_fig,'FaceColor','interp')
    grid on
    view(3)
    axis vis3d
    camlight
    xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%   axis([-0.4 0.4 -0.4 0.4 -0.1 0.1])
%     hlink = linkprop([sp1;sp2], {'CameraPosition','CameraUpVector'});
end
%% clean up
clear data_opt plotYN
assignin('base', 'PCAs', PCAs);
