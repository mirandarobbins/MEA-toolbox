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
for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};
%     fieldname=active_file(1:numel(active_file)-4);
%     meta.WT_exported(file_idx)=struct(fieldname,[]);
    % update progressbar
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    temp=mean(CSD.csd_array(:,:,data.sweep_sort.successful_sweeps),3);
    meta.WT_exported.mean_PDF{file_idx}=temp';
    
    
end
    clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
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
    toc
    temp=mean(CSD.csd_array(:,:,data.sweep_sort.successful_sweeps),3);
    meta.KO_exported.mean_PDF{file_idx}=temp';
    
    
end
    clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% Pre process
clear input
no_WT=numel(meta.archive.WT_filenames);
no_KO=numel(meta.archive.KO_filenames);
time_lims=200:1000;
n_time_lims=numel(time_lims);
for file_idx=1:no_WT
    a=downsample(meta.WT_exported.mean_PDF{file_idx}',20);
    a=a(time_lims,:)';
temp{file_idx}=a;
end; clear file_idx 

input=cell2mat(temp)';
clear temp
for file_idx=1:no_KO

    a=downsample(meta.KO_exported.mean_PDF{file_idx}',20);
    a=a(time_lims,:)';
temp{file_idx}=a;
end; clear file_idx
input2=cell2mat(temp)';

input=vertcat(input,input2);

clear temp input2 a

%% PCA plot 3D
[COEFF,SCORE]=princomp(input);
% SCORE=smooth2a(SCORE,1,1);
figure; hold on
% subplot(1,2,1);hold on
% for trial_id=1:size(input,1)/numel(time_lims)
for trial_id=1:no_WT

plot3(SCORE(1:n_time_lims,1),...
      SCORE(1:n_time_lims,2),...
      SCORE(1:n_time_lims,3),'b','LineWidth',1.5)
SCORE(1:n_time_lims,:)=[];
end
grid on
view(3)
shading interp
camlight

% subplot(1,2,2);hold on
for trial_id=1:no_KO
plot3(SCORE(1:n_time_lims,1),...
      SCORE(1:n_time_lims,2),...
      SCORE(1:n_time_lims,3),'r','LineWidth',1.5)
SCORE(1:n_time_lims,:)=[];
end;
grid on
view(3)
shading interp
camlight

%% PCA plot 2D
[COEFF,SCORE]=princomp(input);
% SCORE=smooth2a(SCORE,1,1);
figure; hold on
% subplot(1,2,1);hold on
% for trial_id=1:size(input,1)/numel(time_lims)
for trial_id=1:no_WT

plot(SCORE(1:n_time_lims,1),...
      SCORE(1:n_time_lims,2),'b','LineWidth',1)
SCORE(1:n_time_lims,:)=[];
end
grid on
view(3)
shading interp
camlight

% subplot(1,2,2);hold on
for trial_id=1:no_KO
plot(SCORE(1:n_time_lims,1),...
      SCORE(1:n_time_lims,2),'r','LineWidth',1)
SCORE(1:n_time_lims,:)=[];
end;
grid on
view(2)
shading interp
camlight