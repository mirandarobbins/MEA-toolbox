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
channel_index=1:64;    
progbar = waitbar(0,'Initializing...','name','file loading progress');

for file_idx=1:numel(meta.archive.WT_filenames)
    active_file=meta.archive.WT_filenames{file_idx};
    % update progressbar
%     waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    [data params spectro] = spectroexploreChronux64_AllChannels(data,params);
    toc
%%% 1 choose mode stongest channel   
    for sweep_id=data.sweep_sort.successful_sweeps
        meta.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
    end; clear sweep_id
    meta.StrongestChannel(meta.StrongestChannel==0)=NaN;
    chosen_channel=mode(meta.StrongestChannel);
    meta.WT_exported.StrongestChannel{file_idx}=meta.StrongestChannel;
    meta.WT_exported.chosen_channel{file_idx}=chosen_channel;

meta.WT_exported.spectro.average_window_power{file_idx}=spectro.average_window_power{chosen_channel};
meta.WT_exported.spectro.average_window_power_mean(:,file_idx)=spectro.average_window_power_mean{chosen_channel};
% 
%     for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
%         this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
% %%% 2 - choose biggest for each trial individually
% %         chosen_channel=channel_index(...
% %                                      min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))...
% %                                      ==min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
% %%% 3 - choose earliest for each trial individually                        
%     end
    meta.WT_exported.spectro.F=spectro.F;
    
    meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx}           = spectro.window_power.gamma_centre_freq;
    meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx}     = spectro.window_power.gamma_centre_freq_power;
    meta.WT_exported.spectro.window_power.gamma_centre_freq_foundYN{file_idx}   = spectro.window_power.gamma_centre_freq_foundYN;
    meta.WT_exported.spectro.window_power.gamma_centre_freq_mean{file_idx}      = spectro.window_power.gamma_centre_freq_mean;
    meta.WT_exported.spectro.window_power.gamma_centre_freq_SEM{file_idx}       = spectro.window_power.gamma_centre_freq_SEM;
end

clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);
channel_index=1:64;    
progbar = waitbar(0,'Initializing...','name','file loading progress');

for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
    % update progressbar
%     waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    [data params spectro] = spectroexploreChronux64_AllChannels(data,params);
    toc
%%% 1 choose mode stongest channel   
    for sweep_id=data.sweep_sort.successful_sweeps
        meta.StrongestChannel(sweep_id)=...
        find(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id))==...
                min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,sweep_id)))...
            );
    end; clear sweep_id
    meta.StrongestChannel(meta.StrongestChannel==0)=NaN;
    chosen_channel=mode(meta.StrongestChannel);
    meta.KO_exported.StrongestChannel{file_idx}=meta.StrongestChannel;
    meta.KO_exported.chosen_channel{file_idx}=chosen_channel;
    
meta.KO_exported.spectro.average_window_power{file_idx}=spectro.average_window_power{chosen_channel};
meta.KO_exported.spectro.average_window_power_mean(:,file_idx)=spectro.average_window_power_mean{chosen_channel};
% 
%     for sweep_id=1:numel(data.sweep_sort.successful_sweeps);
%         this_sweep=data.sweep_sort.successful_sweeps(sweep_id);
% %%% 2 - choose biggest for each trial individually
% %         chosen_channel=channel_index(...
% %                                      min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))...
% %                                      ==min(min(data.filtered_lfp(params.search_win_samples(1):params.search_win_samples(2),:,this_sweep))));
% %%% 3 - choose earliest for each trial individually                        
%     end
    meta.KO_exported.spectro.F=spectro.F;

    meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx}           = spectro.window_power.gamma_centre_freq;
    meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx}     = spectro.window_power.gamma_centre_freq_power;
    meta.KO_exported.spectro.window_power.gamma_centre_freq_foundYN{file_idx}   = spectro.window_power.gamma_centre_freq_foundYN;
    meta.KO_exported.spectro.window_power.gamma_centre_freq_mean{file_idx}      = spectro.window_power.gamma_centre_freq_mean;
    meta.KO_exported.spectro.window_power.gamma_centre_freq_SEM{file_idx}       = spectro.window_power.gamma_centre_freq_SEM;
end


clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% calculate population means... 
meta.WT_exported.spectro.average_window_power_grand_mean=nanmean(meta.WT_exported.spectro.average_window_power_mean,2);
meta.WT_exported.spectro.average_window_power_grand_SEM=nansem(meta.WT_exported.spectro.average_window_power_mean,2);
file_idx=1;
temp_freq  =reshape(meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx},numel(meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx}),1);
temp_power =reshape(meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx},numel(meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx}),1);
for file_idx=2:numel(meta.archive.WT_filenames)
%     clear temp_freq2 temp_power2  
        temp_freq2  =reshape(meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx},numel(meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx}),1);
        temp_power2 =reshape(meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx},numel(meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx}),1);
        temp_freq=vertcat(temp_freq,temp_freq2);
        temp_power=vertcat(temp_power,temp_power2);
end
meta.WT_exported.spectro.window_power.gamma_centre_freq_all_binned=temp_freq;% WT=horzcat(temp_freq,temp_power);
meta.WT_exported.spectro.window_power.gamma_centre_freq_power_all_binned=temp_power;
for file_idx=1:numel(meta.archive.WT_filenames)
meta.WT_exported.spectro.window_power.gamma_centre_freq_pop_mean(file_idx)=meta.WT_exported.spectro.window_power.gamma_centre_freq_mean{file_idx}(meta.WT_exported.chosen_channel{file_idx});
meta.WT_exported.spectro.window_power.gamma_centre_freq_power_pop_mean(file_idx)=meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx}(meta.WT_exported.chosen_channel{file_idx});
end

meta.KO_exported.spectro.average_window_power_grand_mean=nanmean(meta.KO_exported.spectro.average_window_power_mean,2);
meta.KO_exported.spectro.average_window_power_grand_SEM=nansem(meta.KO_exported.spectro.average_window_power_mean,2);
file_idx=1;
temp_freq  =reshape(meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx},numel(meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx}),1);
temp_power =reshape(meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx},numel(meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx}),1);
for file_idx=2:numel(meta.archive.KO_filenames)
%     clear temp_freq2 temp_power2  
        temp_freq2  =reshape(meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx},numel(meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx}),1);
        temp_power2 =reshape(meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx},numel(meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx}),1);
        temp_freq=vertcat(temp_freq,temp_freq2);
        temp_power=vertcat(temp_power,temp_power2);
end
meta.KO_exported.spectro.window_power.gamma_centre_freq_all_binned=temp_freq;% KO=horzcat(temp_freq,temp_power);
meta.KO_exported.spectro.window_power.gamma_centre_freq_power_all_binned=temp_power;
for file_idx=1:numel(meta.archive.KO_filenames)
meta.KO_exported.spectro.window_power.gamma_centre_freq_pop_mean(file_idx)=meta.KO_exported.spectro.window_power.gamma_centre_freq_mean{file_idx}(meta.KO_exported.chosen_channel{file_idx});
meta.KO_exported.spectro.window_power.gamma_centre_freq_power_pop_mean(file_idx)=meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx}(meta.KO_exported.chosen_channel{file_idx});
end

clear file_idx temp_freq temp_freq2 temp_power temp_power2
%% plot
figure; 
subplot(1,2,1);hold on
for file_idx=1:numel(meta.archive.WT_filenames)
% plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power{file_idx},'b')
plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power_mean(:,file_idx),'-b','LineWidth',0.5);
end
% % plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power_grand_mean,'b','LineWidth',2);
% upper=meta.WT_exported.spectro.average_window_power_grand_mean+meta.WT_exported.spectro.average_window_power_grand_SEM;
% lower=meta.WT_exported.spectro.average_window_power_grand_mean-meta.WT_exported.spectro.average_window_power_grand_SEM;
% ciplot(lower,upper,meta.WT_exported.spectro.F,'b')
set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([0 200 -75 -20])
xlabel('Frequency')
ylabel('Spectral Power (dB)')
subplot(1,2,2);hold on
for file_idx=1:numel(meta.archive.KO_filenames)
%     plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power{file_idx},'r')
    plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power_mean(:,file_idx),'-r','LineWidth',0.5);

end
% plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power_grand_mean,'r','LineWidth',2);
% upper=meta.KO_exported.spectro.average_window_power_grand_mean+meta.KO_exported.spectro.average_window_power_grand_SEM;
% lower=meta.KO_exported.spectro.average_window_power_grand_mean-meta.KO_exported.spectro.average_window_power_grand_SEM;
% ciplot(lower,upper,meta.KO_exported.spectro.F,'r')
set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([0 200 -75 -20])
xlabel('Frequency')
ylabel('Spectral Power (dB)')
%% plot overlaid each rep
figure; hold on
for file_idx=1:numel(meta.archive.WT_filenames)
plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power{file_idx},'b')
% plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power_mean(:,file_idx),'-b','LineWidth',0.5);
end
% plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power_grand_mean,'b','LineWidth',2);
% upper=meta.WT_exported.spectro.average_window_power_grand_mean+meta.WT_exported.spectro.average_window_power_grand_SEM;
% lower=meta.WT_exported.spectro.average_window_power_grand_mean-meta.WT_exported.spectro.average_window_power_grand_SEM;
% ciplot(lower,upper,meta.WT_exported.spectro.F,'b')
for file_idx=1:numel(meta.archive.KO_filenames)
    plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power{file_idx},'r')
%     plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power_mean(:,file_idx),'-r','LineWidth',0.5);

end
% plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power_grand_mean,'r','LineWidth',2);
% upper=meta.KO_exported.spectro.average_window_power_grand_mean+meta.KO_exported.spectro.average_window_power_grand_SEM;
% lower=meta.KO_exported.spectro.average_window_power_grand_mean-meta.KO_exported.spectro.average_window_power_grand_SEM;
% ciplot(lower,upper,meta.KO_exported.spectro.F,'r')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% axis([0 200 -75 -20])
xlabel('Frequency')
ylabel('Spectral Power (dB)')
%% plot overlaid averages
figure; hold on
% for file_idx=1:numel(meta.archive.WT_filenames)
% plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power{file_idx},'b')
% plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power_mean(:,file_idx),'-b','LineWidth',0.5);
% end
plot(meta.WT_exported.spectro.F,meta.WT_exported.spectro.average_window_power_grand_mean,'b','LineWidth',2);
upper=meta.WT_exported.spectro.average_window_power_grand_mean+meta.WT_exported.spectro.average_window_power_grand_SEM;
lower=meta.WT_exported.spectro.average_window_power_grand_mean-meta.WT_exported.spectro.average_window_power_grand_SEM;
ciplot(lower,upper,meta.WT_exported.spectro.F,'b')
% for file_idx=1:numel(meta.archive.KO_filenames)
%     plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power{file_idx},'r')
%     plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power_mean(:,file_idx),'-r','LineWidth',0.5);

% end
plot(meta.KO_exported.spectro.F,meta.KO_exported.spectro.average_window_power_grand_mean,'r','LineWidth',2);
upper=meta.KO_exported.spectro.average_window_power_grand_mean+meta.KO_exported.spectro.average_window_power_grand_SEM;
lower=meta.KO_exported.spectro.average_window_power_grand_mean-meta.KO_exported.spectro.average_window_power_grand_SEM;
ciplot(lower,upper,meta.KO_exported.spectro.F,'r')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([0 200 -75 -20])
xlabel('Frequency')
ylabel('Spectral Power (dB)')

%%
figure; hold on


for file_idx=1:numel(meta.archive.WT_filenames)
        temp_freq  =reshape(meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx},numel(meta.WT_exported.spectro.window_power.gamma_centre_freq{file_idx}),1);
        temp_power =reshape(meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx},numel(meta.WT_exported.spectro.window_power.gamma_centre_freq_power{file_idx}),1);
        scatter(temp_freq,temp_power,'b')
        clear temp_freq temp_power
end
for file_idx=1:numel(meta.archive.KO_filenames)
        temp_freq  =reshape(meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx},numel(meta.KO_exported.spectro.window_power.gamma_centre_freq{file_idx}),1);
        temp_power =reshape(meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx},numel(meta.KO_exported.spectro.window_power.gamma_centre_freq_power{file_idx}),1);
        scatter(temp_freq,temp_power,'r')
        clear temp_freq temp_power
end







