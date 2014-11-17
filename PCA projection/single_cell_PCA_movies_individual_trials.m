    [PCAs] = PCA_individual_trials(data,spikes,CSD,1) ;

    %%
    last_dir=pwd;
    if  ~isdir('movies')
        mkdir('movies')
    end
    cd movies
    target_dir=strcat(last_dir,'/movies');



    X=[];Y=[];Z=[];
    for trial_id=1:size(PCAs.LFP.trajectories,2)
         [X{trial_id}...
          Y{trial_id}...
          Z{trial_id}] =tubeplot(PCAs.LFP.trajectories{1,trial_id}(:,1),...
                                 PCAs.LFP.trajectories{1,trial_id}(:,2),...
                                 PCAs.LFP.trajectories{1,trial_id}(:,3),...
                                 0.001);
    end
    M = sweep2movie_PCA3D_multitrace(X,Y,Z,[1 800],'axlims','auto');
    movie2avi(M, ['PCA_each_trial_LFP.avi'], 'compression', 'None')
    close(gcf)
    clear X Y Z M movie trial_id

    X=[];Y=[];Z=[];
     for trial_id=1:size(PCAs.CSD.trajectories,2)
         [X{trial_id}...
          Y{trial_id}...
          Z{trial_id}] =tubeplot(PCAs.CSD.trajectories{1,trial_id}(:,1),...
                                 PCAs.CSD.trajectories{1,trial_id}(:,2),...
                                 PCAs.CSD.trajectories{1,trial_id}(:,3),...
                                 0.5);
     end
    M = sweep2movie_PCA3D_multitrace(X,Y,Z,[1 800],'axlims','auto');
    movie2avi(M,['PCA_each_trial_CSD.avi'], 'compression', 'None')
    close(gcf)
    clear X Y Z M movie trial_id


    X=[];Y=[];Z=[];
     for trial_id=1:size(PCAs.MUA.trajectories,2)
         [X{trial_id}...
          Y{trial_id}...
          Z{trial_id}] =tubeplot(PCAs.MUA.trajectories{1,trial_id}(:,1),...
                                 PCAs.MUA.trajectories{1,trial_id}(:,2),...
                                 PCAs.MUA.trajectories{1,trial_id}(:,3),...
                                 0.002);
     end
    M = sweep2movie_PCA3D_multitrace(X,Y,Z,[1 800],'axlims','auto');
    movie2avi(M, ['PCA_each_trial_MUA.avi'], 'compression', 'None')
    close(gcf)
    clear X Y Z M movie trial_id

    cd(last_dir) 
    clear  X Y Z M movie trial_id input last_dir target_dir
    disp( 'finished coding AVI files!')