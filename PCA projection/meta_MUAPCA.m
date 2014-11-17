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
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_id=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_id};
%     fieldname=active_file(1:numel(active_file)-4);
%     meta.WT_exported(file_id)=struct(fieldname,[]);
    % update progressbar
    waitbar(file_id/no_files,progbar,strcat(['Loading file ',num2str(file_id),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    spikes = meta_spike_extract(data,params,3);
    spikes = meta_MUA64plot(data,params,spikes);
    meta.WT_exported.mean_PDF{file_id}=spikes.PDF.mean_PDF;
    meta.WT_exported.each_trial_PDF{file_id}=spikes.PDF.PDF_trimmed;
end
    clear data params burst spectro spikes CSD 
clear active_file file_id no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_id=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_id};
%     fieldname=active_file(1:numel(active_file)-4);
%     meta.KO_exported(file_id)=struct(fieldname,[]);
    % update progressbar
    waitbar(file_id/no_files,progbar,strcat(['Loading file ',num2str(file_id),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    spikes = meta_spike_extract(data,params,3);
    spikes = meta_MUA64plot(data,params,spikes);
    meta.KO_exported.mean_PDF{file_id}=spikes.PDF.mean_PDF;
    meta.KO_exported.each_trial_PDF{file_id}=spikes.PDF.PDF_trimmed;
end
    clear data params burst spectro spikes CSD 
clear active_file file_id no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% Pre process
clear temp
meta.PCA=[];
meta.no_WT=numel(meta.archive.WT_filenames);
meta.no_KO=numel(meta.archive.KO_filenames);
meta.PCA.time_lims=200:1000;
meta.PCA.no_time_steps=numel(meta.PCA.time_lims);
for file_id=1:meta.no_WT
    temp{file_id}=meta.WT_exported.mean_PDF{file_id}(:,meta.PCA.time_lims)./...
                  max(max(meta.WT_exported.mean_PDF{file_id}(:,meta.PCA.time_lims)));
    temp{file_id}(temp{file_id}==NaN)=0;
end; clear file_id
input=cell2mat(temp)';
clear temp

for file_id=1:meta.no_KO
    temp{file_id}=meta.KO_exported.mean_PDF{file_id}(:,meta.PCA.time_lims)./...
                  max(max(meta.KO_exported.mean_PDF{file_id}(:,meta.PCA.time_lims)));
    temp{file_id}(temp{file_id}==NaN)=0;              
end; clear file_id
input2=cell2mat(temp)';

meta.PCA.input=vertcat(input,input2);
meta.PCA.input(isnan(meta.PCA.input))=0;
input(isnan(input))=0;
input2(isnan(input2))=0;
clear temp a file_id
% clear input input2 
%% for WT and KO combined
% compute PCA, sepetate back out again and plot mean traces overlaid
[meta.PCA.COEFF,...
 meta.PCA.SCORE,...
 meta.PCA.var_explained] = princomp(meta.PCA.input);
meta.PCA.percent_var_explained = 100*meta.PCA.var_explained/sum(meta.PCA.var_explained);

figure; pareto(meta.PCA.percent_var_explained )
figure; hold on

for file_id=1:meta.no_WT
    for PC_id=1:64
        meta.PCA.WT_PCA.PCs{PC_id}(:,file_id)=meta.PCA.SCORE(1:meta.PCA.no_time_steps,PC_id);
    end
    [X Y Z] =tubeplot(meta.PCA.SCORE(1:meta.PCA.no_time_steps,1),...
                      meta.PCA.SCORE(1:meta.PCA.no_time_steps,2),...
                      meta.PCA.SCORE(1:meta.PCA.no_time_steps,3),0.01);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
    meta.PCA.SCORE(1:meta.PCA.no_time_steps,:)=[];
end

for file_id=1:meta.no_KO
    for PC_id=1:64
        meta.PCA.KO_PCA.PCs{PC_id}(:,file_id)=meta.PCA.SCORE(1:meta.PCA.no_time_steps,PC_id);
    end
    [X Y Z] =tubeplot(meta.PCA.SCORE(1:meta.PCA.no_time_steps,1),...
                      meta.PCA.SCORE(1:meta.PCA.no_time_steps,2),...
                      meta.PCA.SCORE(1:meta.PCA.no_time_steps,3),0.01);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')
    meta.PCA.SCORE(1:meta.PCA.no_time_steps,:)=[];
end
axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
grid on
view(3)
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
clear X Y Z WT_fig KO_fig file_id PC_id
%% for WT and KO individually
% compute PCA, sepetate back out again and plot mean traces overlaid
[meta.PCA.COEFF_WT,...
 meta.PCA.SCORE_WT,...
 meta.PCA.var_explained_WT] = princomp(input);
 meta.PCA.percent_var_explained_WT = 100*meta.PCA.var_explained_WT/sum(meta.PCA.var_explained_WT);

 [meta.PCA.COEFF_KO,...
 meta.PCA.SCORE_KO,...
 meta.PCA.var_explained_KO] = princomp(input2);
 meta.PCA.percent_var_explained_KO = 100*meta.PCA.var_explained_KO/sum(meta.PCA.var_explained_KO);



figure; hold on
for file_id=1:meta.no_WT
    for PC_id=1:64
        meta.PCA.WT_PCA.PCs{PC_id}(:,file_id)=meta.PCA.SCORE_WT(1:meta.PCA.no_time_steps,PC_id);
    end
    [X Y Z] =tubeplot(meta.PCA.SCORE_WT(1:meta.PCA.no_time_steps,1),...
                      meta.PCA.SCORE_WT(1:meta.PCA.no_time_steps,2),...
                      meta.PCA.SCORE_WT(1:meta.PCA.no_time_steps,3),0.01);
    WT_fig=surf(X,Y,Z);
    set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
    meta.PCA.SCORE_WT(1:meta.PCA.no_time_steps,:)=[];
end

for file_id=1:meta.no_KO
    for PC_id=1:64
        meta.PCA.KO_PCA.PCs{PC_id}(:,file_id)=meta.PCA.SCORE_KO(1:meta.PCA.no_time_steps,PC_id);
    end
    [X Y Z] =tubeplot(meta.PCA.SCORE_KO(1:meta.PCA.no_time_steps,1),...
                      meta.PCA.SCORE_KO(1:meta.PCA.no_time_steps,2),...
                      meta.PCA.SCORE_KO(1:meta.PCA.no_time_steps,3),0.01);
    KO_fig=surf(X,Y,Z);
    set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')
    meta.PCA.SCORE_KO(1:meta.PCA.no_time_steps,:)=[];
end
axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
grid on
view(3)
camlight
xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
clear X Y Z WT_fig KO_fig file_id PC_id

figure;
subplot(1,2,1)
pareto(meta.PCA.var_explained_WT)
subplot(1,2,2)
pareto(meta.PCA.var_explained_KO)

 
%%  - plot pareto (variance explained)
clear WT_fig WT_fig_child
h_fig=figure;
% temp=horzcat(percent_explained_WT,percent_explained_KO,percent_explained_WTKO);
hold on
WT_fig=pareto(meta.PCA.percent_var_explained_WT);
WT_fig_child(1)=get(WT_fig(1),'Children');
    set(WT_fig_child(1),'facealpha',0.99,'facecolor','b');
    set(WT_fig(2),'color','b','LineWidth',1.5);
KO_fig=pareto(meta.PCA.percent_var_explained_KO)
    KO_fig_child=get(KO_fig,'Children');
    set(KO_fig_child{1},'facealpha',0.99,'facecolor','r');
    set(KO_fig(2),'color','r','LineWidth',1.5);
% WTKO_fig=pareto(percent_explained_WTKO)
%     WTKO_fig_child=get(WTKO_fig,'Children');
%     set(WTKO_fig_child{1},'facealpha',1,'facecolor','g');
%     set(WTKO_fig_child{2},'color','g');

% bar(temp,'stacked')
% colormap([0 0 1;1 0 0; 0 1 0])
axis([0 10 0 100])
all_h=get(h_fig,'Children');
set(all_h,'box','off','FontSize',14)
xlabel('Principle Component no.')
ylabel('% Variance explained')
%% PCA individual trial subplots
% [COEFF,SCORE]=princomp(input);
% % SCORE=smooth2a(SCORE,1,1);
% figure; hold on
% % for file_id=1:size(input,1)/numel(time_lims)
% for file_id=1:no_WT
% subplot(ceil(no_WT^0.5),ceil(no_WT^0.5),file_id);hold on
% 
% [X Y Z] =tubeplot(SCORE(1:n_time_steps,1),...
%          SCORE(1:n_time_steps,2),...
%          SCORE(1:n_time_steps,3),0.01);
% WT_fig=surf(X,Y,Z);
% set(WT_fig,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')
%   SCORE(1:n_time_steps,:)=[];
% axis([0 2 -1 1 -1 1])
% grid on
% view(3)
% % axis vis3d
% camlight
% end
% figure; hold on
% 
% for file_id=1:no_KO
% subplot(ceil(no_KO^0.5),ceil(no_KO^0.5),file_id);hold on
% 
% [X Y Z] =tubeplot(SCORE(1:n_time_steps,1),...
%          SCORE(1:n_time_steps,2),...
%          SCORE(1:n_time_steps,3),0.01);
% KO_fig=surf(X,Y,Z);
% set(KO_fig,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')
%   SCORE(1:n_time_steps,:)=[];
% axis([0 2 -1 1 -1 1])
% grid on
% view(3)
% % axis vis3d
% camlight
% end
% xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3');

%% Mean responses of all PCs
for PC_id=1:64
    meta.PCA.WT_PCA.PC_mean(:,PC_id)=nanmean(meta.PCA.WT_PCA.PCs{PC_id},2);
    meta.PCA.WT_PCA.PC_SEM(:,PC_id) =nansem(meta.PCA.WT_PCA.PCs{PC_id},2);
    
    meta.PCA.KO_PCA.PC_mean(:,PC_id)=nanmean(meta.PCA.KO_PCA.PCs{PC_id},2);
    meta.PCA.KO_PCA.PC_SEM(:,PC_id) =nansem(meta.PCA.KO_PCA.PCs{PC_id},2);
end
figure; hold on
plot3(meta.PCA.WT_PCA.PC_mean(:,1),meta.PCA.WT_PCA.PC_mean(:,2),meta.PCA.WT_PCA.PC_mean(:,3),    'b','LineWidth',1.5);
plot3(meta.PCA.KO_PCA.PC_mean(:,1),meta.PCA.KO_PCA.PC_mean(:,2),meta.PCA.KO_PCA.PC_mean(:,3),    'r','LineWidth',1.5);
grid on
clear PC_id
%% calculate expansion coefficients in 1st 3 PCs 
for file_id=1:meta.no_WT
    for tb_id=1:meta.PCA.no_time_steps
        tstep=meta.PCA.time_lims(tb_id);
        pos1= [meta.PCA.WT_PCA.PCs{1}(tstep,file_id),...
               meta.PCA.WT_PCA.PCs{2}(tstep,file_id),...
               meta.PCA.WT_PCA.PCs{2}(tstep,file_id)];
        pos2= [meta.PCA.WT_PCA.PCs{1}(tstep+1,file_id),...
               meta.PCA.WT_PCA.PCs{2}(tstep+1,file_id),...
               meta.PCA.WT_PCA.PCs{2}(tstep+1,file_id)];
        vec = pos2-pos1;
        meta.PCA.WT_PCvector(tstep,file_id)=abs(sum(vec.^2)^0.5);  
    end
end

for file_id=1:meta.no_KO
    for tb_id=1:meta.PCA.no_time_steps
        tstep=meta.PCA.time_lims(tb_id);
        pos1= [meta.PCA.KO_PCA.PCs{1}(tstep,file_id),...
               meta.PCA.KO_PCA.PCs{2}(tstep,file_id),...
               meta.PCA.KO_PCA.PCs{2}(tstep,file_id)];
        pos2= [meta.PCA.KO_PCA.PCs{1}(tstep+1,file_id),...
               meta.PCA.KO_PCA.PCs{2}(tstep+1,file_id),...
               meta.PCA.KO_PCA.PCs{2}(tstep+1,file_id)];
        vec = pos2-pos1;
        meta.PCA.KO_PCvector(tstep,file_id)=abs(sum(vec.^2)^0.5);  
    end
end

WT=cumsum(meta.PCA.WT_PCvector)./max(cumsum(meta.PCA.WT_PCvector),2);
KO=cumsum(meta.PCA.KO_PCvector)./max(cumsum(meta.PCA.KO_PCvector),2);

figure; 
subplot(2,1,1);hold on
plot(meta.PCA.WT_PCvector,'b')
plot(meta.PCA.KO_PCvector,'r')
subplot(2,1,2);hold on
plot(WT,'b')
plot(KO,'r')
clear file_id tb_id tstep pos1 pos2 vec 








