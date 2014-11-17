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
    
    % update progressbar & load file to workspace
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    close all
    [spikes]                     = meta_spike_extract(data,params,4);                    % blank peri-stim only
    [spikes]                     = meta_MUA64plot(data,params,spikes);                   % check kernel width of 1ms      
    [data params spectro]        = spectroexploreChronux64_AllChannels(data,params);
    [data params spikes spectro] = MEAChannel_coherence(data,spikes,params,spectro);
    [data params spikes spectro] = shared_coherence(data,spikes,params,spectro);
    
    
    % choose mode strongest channel   
    chosen_channel=spectro.channeltoanalyse;
    
%%%%% LFPLFP
    % extract mean LFPLFP coherence/phase across co-active channel radius
    meta.WT_exported.spectro.f_LFPLFP=spectro.coherence_LFPLFP.f;
    meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean{file_idx}    =  spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean;
    meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean{file_idx}        =  spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean;
    
    % extract mean LFPLFP coherence/phase at peak gamma freq
    meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean{file_idx}               =  spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean;
    meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean{file_idx}                   =  spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean;

%%%%% LFPspike
    % extract spike-LFP coherence/phase on strongest channel
    meta.WT_exported.spectro.f_LFPspike=spectro.coherence_LFPspike.f;
    for trial_id=data.sweep_sort.successful_sweeps
        temp_C(:,trial_id)   = spectro.coherence_LFPspike.C_mean(:,chosen_channel(trial_id));
        temp_phi(:,trial_id) = spectro.coherence_LFPspike.phi_mean(:,chosen_channel(trial_id));
    end
    temp_C(:,var(temp_C)==0)=NaN; temp_phi(:,var(temp_phi)==0)=NaN;
    meta.WT_exported.spectro.LFPspike.C_mean{file_idx}      =  nanmean(temp_C,2);
    meta.WT_exported.spectro.LFPspike.phi_mean{file_idx}    =  nanmean(temp_phi,2);
   
    % extract mean LFPspike coherence/phase across co-active channel radius
    meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean{file_idx}    =  spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean;
    meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean{file_idx}        =  spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean;
    
    % extract mean LFPspike coherence/phase at peak gamma freq
    meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean{file_idx}               =  spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean;
    meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean{file_idx}                   =  spectro.coherence_radius.LFPspike.peak_gamma_phase_mean;
end
% clear data params burst spectro spikes CSD 
% clear active_file channel_index chosen_channel file_idx no_files progbar 
% clear functions
% clear hidden
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);
channel_index=1:64;    
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
    
    % update progressbar & load file to workspace
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    close all
    [data params spectro] = spectroexploreChronux64_AllChannels(data,params);
    [data params spikes spectro] = MEAChannel_coherence(data,spikes,params,spectro);
    [data params spikes spectro] = shared_coherence(data,spikes,params,spectro);
    
    
    % choose mode strongest channel   
    chosen_channel=spectro.channeltoanalyse;
    
%%%%% LFPLFP
    % extract mean LFPLFP coherence/phase across co-active channel radius
    meta.KO_exported.spectro.f_LFPLFP=spectro.coherence_LFPLFP.f;
    meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean{file_idx}    =  spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean;
    meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean{file_idx}        =  spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean;
    
    % extract mean LFPLFP coherence/phase at peak gamma freq
    meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean{file_idx}               =  spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean;
    meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean{file_idx}                   =  spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean;

%%%%% LFPspike
    % extract spike-LFP coherence/phase on strongest channel
    meta.KO_exported.spectro.f_LFPspike=spectro.coherence_LFPspike.f;
    for trial_id=data.sweep_sort.successful_sweeps
        temp_C(:,trial_id)   = spectro.coherence_LFPspike.C_mean(:,chosen_channel(trial_id));
        temp_phi(:,trial_id) = spectro.coherence_LFPspike.phi_mean(:,chosen_channel(trial_id));
    end
    temp_C(:,var(temp_C)==0)=NaN; temp_phi(:,var(temp_phi)==0)=NaN;
    meta.KO_exported.spectro.LFPspike.C_mean{file_idx}      =  nanmean(temp_C,2);
    meta.KO_exported.spectro.LFPspike.phi_mean{file_idx}    =  nanmean(temp_phi,2);
   
    % extract mean LFPspike coherence/phase across co-active channel radius
    meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean{file_idx}    =  spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean;
    meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean{file_idx}        =  spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean;
    
    % extract mean LFPspike coherence/phase at peak gamma freq
    meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean{file_idx}               =  spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean;
    meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean{file_idx}                   =  spectro.coherence_radius.LFPspike.peak_gamma_phase_mean;
end
clear data params burst spectro spikes CSD temp_C temp_phi
clear active_file channel_index chosen_channel file_idx no_files progbar 
clear functions
clear hidden
%% Population stats (WT)
%%%%%%LFP-LFP coherence
% LFPLFP spectra over co-active radius
meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean  =  nanmean(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean'));
meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_SEM   =  nansem(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean'));
meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean      =  nanmean(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean'));
meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_SEM       =  nansem(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean'));

% LFPLFP coherence/phase at peak gamma
meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_grand_mean             =  nanmean(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean));
meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_grand_SEM              =  nansem(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean));
meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean_cat               =  cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean);
meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean_cat                   =  cell2mat(meta.WT_exported.spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean)*180/pi;
%%%%%%LFPspike coherence
% LFPspike spectra over co-active radius
meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean  =  nanmean(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean'));
meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_SEM   =  nansem(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean'));
meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean      =  nanmean(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean'));
meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_SEM       =  nansem(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean'));

% LFPspike coherence/phase at peak gamma
meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_grand_mean             =  nanmean(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean));
meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_grand_SEM              =  nansem(cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean));
meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean_cat               =  cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean);
meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean_cat                   =  cell2mat(meta.WT_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean)*180/pi;

% LFPspike spectra for strongest channel
meta.WT_exported.spectro.LFPspike.C_grand_mean = nanmean(cell2mat(meta.WT_exported.spectro.LFPspike.C_mean),2);
meta.WT_exported.spectro.LFPspike.C_grand_SEM  = nansem(cell2mat(meta.WT_exported.spectro.LFPspike.C_mean),2);
meta.WT_exported.spectro.LFPspike.phi_grand_mean = nanmean(cell2mat(meta.WT_exported.spectro.LFPspike.phi_mean),2);
meta.WT_exported.spectro.LFPspike.phi_grand_SEM  = nansem(cell2mat(meta.WT_exported.spectro.LFPspike.phi_mean),2);



%% Population stats (KO)
%%%%%%LFP-LFP coherence
% LFPLFP spectra over co-active radius
meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean  =  nanmean(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean'));
meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_SEM   =  nansem(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean'));
meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean      =  nanmean(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean'));
meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_SEM       =  nansem(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean'));

% LFPLFP coherence/phase at peak gamma
meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_grand_mean             =  nanmean(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean));
meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_grand_SEM              =  nansem(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean));
meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean_cat               =  cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_coherence_mean);
meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean_cat                   =  cell2mat(meta.KO_exported.spectro.coherence_radius.LFPLFP.peak_gamma_phase_mean)*180/pi;
%%%%%%LFPspike coherence
% LFPspike spectra over co-active radius
meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean  =  nanmean(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean'));
meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_SEM   =  nansem(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_mean'));
meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean      =  nanmean(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean'));
meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_SEM       =  nansem(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_mean'));

% LFPspike coherence/phase at peak gamma
meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_grand_mean             =  nanmean(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean));
meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_grand_SEM              =  nansem(cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean));
meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean_cat               =  cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_coherence_mean);
meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean_cat                   =  cell2mat(meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean)*180/pi;

% LFPspike spectra for strongest channel
meta.KO_exported.spectro.LFPspike.C_grand_mean = nanmean(cell2mat(meta.KO_exported.spectro.LFPspike.C_mean),2);
meta.KO_exported.spectro.LFPspike.C_grand_SEM  = nansem(cell2mat(meta.KO_exported.spectro.LFPspike.C_mean),2);
meta.KO_exported.spectro.LFPspike.phi_grand_mean = nanmean(cell2mat(meta.KO_exported.spectro.LFPspike.phi_mean),2);
meta.KO_exported.spectro.LFPspike.phi_grand_SEM  = nansem(cell2mat(meta.KO_exported.spectro.LFPspike.phi_mean),2);


%% LFPspike phase histograms
figure
wind_rose2(meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean_cat,ones(numel(meta.KO_exported.spectro.coherence_radius.LFPspike.peak_gamma_phase_mean_cat)))

%% Plot LFPLFP spectra over co-active radius
% 1- coherence
figure('name','Mean LFP-LFP coherence over active radius'); 
subplot(1,2,1);hold on
    plot(meta.WT_exported.spectro.f_LFPLFP,smooth_hist(meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean),'b')
    upper=meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean+...
          meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_SEM;
    lower=meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean-...
          meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_SEM;
    ciplot(smooth_hist(lower),smooth_hist(upper),meta.WT_exported.spectro.f_LFPLFP,'b')
%     plot(meta.WT_exported.spectro.f_LFPLFP,cell2mat((meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean)')','b')

    plot(meta.KO_exported.spectro.f_LFPLFP,smooth_hist(meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean),'r')
    upper=meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean+...
          meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_SEM;
    lower=meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_mean-...
          meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_grand_SEM;
    ciplot(smooth_hist(lower),smooth_hist(upper),meta.WT_exported.spectro.f_LFPLFP,'r')
%     plot(meta.KO_exported.spectro.f_LFPLFP,cell2mat((meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_coherence_mean)')','r')
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
    axis([10 200 0 1])
    xlabel('Frequency (Hz)')
    ylabel('Coherence magnitude')
% 2- phase
subplot(1,2,2);hold on
    plot(meta.WT_exported.spectro.f_LFPLFP,smooth_hist(meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean*180/pi),'b')
    plot(meta.KO_exported.spectro.f_LFPLFP,smooth_hist(meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean*180/pi),'r')
    upper=meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean+meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_SEM;
    lower=meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean-meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_SEM;
    ciplot(smooth_hist(lower*180/pi),smooth_hist(upper*180/pi),meta.WT_exported.spectro.f_LFPLFP,'b')
%     plot(meta.WT_exported.spectro.f_LFPLFP,cell2mat((meta.WT_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean)')'*180/pi,'b')

    upper=meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean+meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_SEM;
    lower=meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_mean-meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_grand_SEM;
    ciplot(smooth_hist(lower*180/pi),smooth_hist(upper*180/pi),meta.WT_exported.spectro.f_LFPLFP,'r')
%     plot(meta.KO_exported.spectro.f_LFPLFP,cell2mat((meta.KO_exported.spectro.coherence_radius.LFPLFP.radius_average.LFPLFP_phase_mean)')'*180/pi,'r')
    axis([10 200 -90 90])
    xlabel('Frequency (Hz)')
    ylabel('Coherence phase (degrees)')


%% Plot LFPspike spectra over co-active radius
% 1- coherence
figure('name','Mean spike-field coherence over active radius'); 
subplot(1,2,1);hold on
plot(meta.WT_exported.spectro.f_LFPspike,(meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean),'b')
plot(meta.KO_exported.spectro.f_LFPspike,(meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean),'r')
upper=meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean+meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_SEM;
lower=meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean-meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_SEM;
ciplot((lower),(upper),meta.WT_exported.spectro.f_LFPspike,'b')

upper=meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean+meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_SEM;
lower=meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_mean-meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_coherence_grand_SEM;
ciplot((lower),(upper),meta.WT_exported.spectro.f_LFPspike,'r')
axis([10 200 0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence magnitude')
% 2- phase
% figure; 
subplot(1,2,2);hold on
plot(meta.WT_exported.spectro.f_LFPspike,(meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean*180/pi),'b')
plot(meta.KO_exported.spectro.f_LFPspike,(meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean*180/pi),'r')
upper=meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean+meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_SEM;
lower=meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean-meta.WT_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_SEM;
ciplot((lower*180/pi),(upper*180/pi),meta.WT_exported.spectro.f_LFPspike,'b')

upper=meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean+meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_SEM;
lower=meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_mean-meta.KO_exported.spectro.coherence_radius.LFPspike.radius_average.LFPspike_phase_grand_SEM;
ciplot((lower*180/pi),(upper*180/pi),meta.WT_exported.spectro.f_LFPspike,'r')
axis([10 200 -90 90])
xlabel('Frequency (Hz)')
ylabel('Coherence phase (degrees)')
%% Plot LFP-spike on strongest channel
figure('name','Mean spike-field coherence at burst centre'); 
subplot(1,2,1);hold on
plot(meta.WT_exported.spectro.f_LFPspike,meta.WT_exported.spectro.LFPspike.C_grand_mean,'b')
plot(meta.KO_exported.spectro.f_LFPspike,meta.KO_exported.spectro.LFPspike.C_grand_mean,'r')

upper=meta.WT_exported.spectro.LFPspike.C_grand_mean+meta.WT_exported.spectro.LFPspike.C_grand_SEM;
lower=meta.WT_exported.spectro.LFPspike.C_grand_mean-meta.WT_exported.spectro.LFPspike.C_grand_SEM;
    ciplot((lower),(upper),meta.WT_exported.spectro.f_LFPspike,'b')
upper=meta.KO_exported.spectro.LFPspike.C_grand_mean+meta.KO_exported.spectro.LFPspike.C_grand_SEM;
lower=meta.KO_exported.spectro.LFPspike.C_grand_mean-meta.KO_exported.spectro.LFPspike.C_grand_SEM;
    ciplot((lower),(upper),meta.KO_exported.spectro.f_LFPspike,'r')
axis([10 200 0 1])
xlabel('Frequency (Hz)')
ylabel('Coherence magnitude')
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
set(gca,'color',[0 1 0])
% figure; 
subplot(1,2,2);hold on
plot(meta.WT_exported.spectro.f_LFPspike,(meta.WT_exported.spectro.LFPspike.phi_grand_mean*180/pi),'b')
plot(meta.KO_exported.spectro.f_LFPspike,(meta.KO_exported.spectro.LFPspike.phi_grand_mean*180/pi),'r')

upper=meta.WT_exported.spectro.LFPspike.phi_grand_mean+meta.WT_exported.spectro.LFPspike.phi_grand_SEM;
lower=meta.WT_exported.spectro.LFPspike.phi_grand_mean-meta.WT_exported.spectro.LFPspike.phi_grand_SEM;
ciplot((lower*180/pi),(upper*180/pi),meta.WT_exported.spectro.f_LFPspike,'b')

upper=meta.KO_exported.spectro.LFPspike.phi_grand_mean+meta.KO_exported.spectro.LFPspike.phi_grand_SEM;
lower=meta.KO_exported.spectro.LFPspike.phi_grand_mean-meta.KO_exported.spectro.LFPspike.phi_grand_SEM;
ciplot((lower*180/pi),(upper*180/pi),meta.KO_exported.spectro.f_LFPspike,'r')
axis([10 200 0 180])
xlabel('Frequency (Hz)')
ylabel('Coherence phase (degrees)')
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
set(gca,'color',[0 1 0])