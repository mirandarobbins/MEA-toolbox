%% Indexes exported data directory and iteratively loads files for measurements
% meta.archive.WT_channels_chosen=[]
% meta.archive.KO_channels_chosen=[]
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
    L4=meta.archive.WT_channels_chosen(file_idx,1);
    L23=meta.archive.WT_channels_chosen(file_idx,2);
    % update progressbar & load file to workspace
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    close all
    [spikes]                     = meta_spike_extract(data,params,4);                    % blank peri-stim only
    [spikes]                     = meta_MUA64plot(data,params,spikes);                   % check kernel width of 1ms      
%     [data params spikes spectro] = MEAChannel_coherence(data,spikes,params,spectro);
    [data params spikes spectro] = pairwise_coherence(data,spikes,params,L4,L23);
        
    
    
%%%%% LFPLFP
    meta.WT_exported.spectro.LFPLFP_L4_signal{file_idx}          = spectro.coherence_LFPLFP.input1;
    meta.WT_exported.spectro.LFPLFP_L23_signal{file_idx}         = spectro.coherence_LFPLFP.input2;
    meta.WT_exported.spectro.f_LFPLFP                            = spectro.coherence_LFPLFP.f;
    meta.WT_exported.spectro.LFPLFP_coherence(file_idx,:)        = spectro.coherence_LFPLFP.C_mean(:,L23);
    meta.WT_exported.spectro.LFPLFP_phase(file_idx,:)            = spectro.coherence_LFPLFP.phi_mean(:,L23);
    meta.WT_exported.spectro.LFP_spectrum_L4(file_idx,:)         = spectro.coherence_LFPLFP.S1_mean;
    meta.WT_exported.spectro.LFP_spectrum_L23(file_idx,:)        = spectro.coherence_LFPLFP.S2_mean;
%%%%% LFPspike
    % extract spike-LFP coherence/phase on strongest channel
    meta.WT_exported.spectro.f_LFPspike=spectro.coherence_LFPspike.f;
    for trial_id=data.sweep_sort.successful_sweeps
        temp_C_L4(:,trial_id)   = spectro.coherence_LFPspike.C_mean(:,L4);
        temp_phi_L4(:,trial_id) = spectro.coherence_LFPspike.phi_mean(:,L4);
        temp_C_L23(:,trial_id)   = spectro.coherence_LFPspike.C_mean(:,L23);
        temp_phi_L23(:,trial_id) = spectro.coherence_LFPspike.phi_mean(:,L23);
    end
    temp_C_L4(:,var(temp_C_L4)==0)=NaN; temp_phi_L4(:,var(temp_phi_L4)==0)=NaN;
    temp_C_L23(:,var(temp_C_L23)==0)=NaN; temp_phi_L23(:,var(temp_phi_L23)==0)=NaN;
    meta.WT_exported.spectro.LFPspike_coherence_L4(file_idx,:)      =  nanmean(temp_C_L4,2);
    meta.WT_exported.spectro.LFPspike_phase_L4(file_idx,:)          =  nanmean(temp_phi_L4,2);
    meta.WT_exported.spectro.LFPspike_coherence_L23(file_idx,:)     =  nanmean(temp_C_L23,2);
    meta.WT_exported.spectro.LFPspike_phase_L23(file_idx,:)         =  nanmean(temp_phi_L23,2);
   
    figure; hold on
    plot(meta.WT_exported.spectro.f_LFPLFP,meta.WT_exported.spectro.LFPLFP_coherence)
    set(gca,'XScale','log');axis([5 200 -Inf 1])


end
% clear data params burst spectro spikes CSD 
% clear active_file channel_index chosen_channel file_idx no_files progbar 
% clear functions
% clear hidden

% collate individual spectra
meta.WT_exported.spectro.LFP_spectrum_L4_mean    = nanmean(meta.WT_exported.spectro.LFP_spectrum_L4);
meta.WT_exported.spectro.LFP_spectrum_L4_SEM     = nansem(meta.WT_exported.spectro.LFP_spectrum_L4);
meta.WT_exported.spectro.LFP_spectrum_L23_mean   = nanmean(meta.WT_exported.spectro.LFP_spectrum_L23);
meta.WT_exported.spectro.LFP_spectrum_L23_SEM    = nansem(meta.WT_exported.spectro.LFP_spectrum_L23);

% collate LFP-LPF coherence 
meta.WT_exported.spectro.LFPLFP_coherence_mean   = nanmean(meta.WT_exported.spectro.LFPLFP_coherence);
meta.WT_exported.spectro.LFPLFP_coherence_SEM    = nansem(meta.WT_exported.spectro.LFPLFP_coherence);
meta.WT_exported.spectro.LFPLFP_phase_mean       = nanmean(meta.WT_exported.spectro.LFPLFP_phase);
meta.WT_exported.spectro.LFPLFP_phase_SEM        = nansem(meta.WT_exported.spectro.LFPLFP_phase);
% collate spike-field coherence 
meta.WT_exported.spectro.LFPspike_coherence_L4_mean = nanmean(meta.WT_exported.spectro.LFPspike_coherence_L4);
meta.WT_exported.spectro.LFPspike_coherence_L4_SEM  = nansem(meta.WT_exported.spectro.LFPspike_coherence_L4);
meta.WT_exported.spectro.LFPspike_phase_L4_mean     = nanmean(meta.WT_exported.spectro.LFPspike_phase_L4);
meta.WT_exported.spectro.LFPspike_phase_L4_SEM      = nansem(meta.WT_exported.spectro.LFPspike_phase_L4);

meta.WT_exported.spectro.LFPspike_coherence_L23_mean = nanmean(meta.WT_exported.spectro.LFPspike_coherence_L23);
meta.WT_exported.spectro.LFPspike_coherence_L23_SEM  = nansem(meta.WT_exported.spectro.LFPspike_coherence_L23);
meta.WT_exported.spectro.LFPspike_phase_L23_mean     = nanmean(meta.WT_exported.spectro.LFPspike_phase_L23);
meta.WT_exported.spectro.LFPspike_phase_L23_SEM      = nansem(meta.WT_exported.spectro.LFPspike_phase_L23);
%% KO load loop
tic;
no_files=numel(meta.archive.KO_filenames);
channel_index=1:64;    
progbar = waitbar(0,'Initializing...','name','file loading progress');
for file_idx=1:numel(meta.archive.KO_filenames)
    active_file=meta.archive.KO_filenames{file_idx};
    L4=meta.archive.KO_channels_chosen(file_idx,1);
    L23=meta.archive.KO_channels_chosen(file_idx,2);
    % update progressbar & load file to workspace
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    toc
    close all
    [spikes]                     = meta_spike_extract(data,params,4);                    % blank peri-stim only
    [spikes]                     = meta_MUA64plot(data,params,spikes);                   % check kernel width of 1ms      
%     [data params spikes spectro] = MEAChannel_coherence(data,spikes,params,spectro);
    [data params spikes spectro] = pairwise_coherence(data,spikes,params,L4,L23);
        
    
    
%%%%% LFPLFP
    meta.KO_exported.spectro.LFPLFP_L4_signal{file_idx}          = spectro.coherence_LFPLFP.input1;
    meta.KO_exported.spectro.LFPLFP_L23_signal{file_idx}         = spectro.coherence_LFPLFP.input2;
    meta.KO_exported.spectro.f_LFPLFP                            = spectro.coherence_LFPLFP.f;
    meta.KO_exported.spectro.LFPLFP_coherence(file_idx,:)        = spectro.coherence_LFPLFP.C_mean(:,L23);
    meta.KO_exported.spectro.LFPLFP_phase(file_idx,:)            = spectro.coherence_LFPLFP.phi_mean(:,L23);
    meta.KO_exported.spectro.LFP_spectrum_L4(file_idx,:)         = spectro.coherence_LFPLFP.S1_mean;
    meta.KO_exported.spectro.LFP_spectrum_L23(file_idx,:)        = spectro.coherence_LFPLFP.S2_mean;
%%%%% LFPspike
    % extract spike-LFP coherence/phase on strongest channel
    meta.KO_exported.spectro.f_LFPspike=spectro.coherence_LFPspike.f;
    for trial_id=data.sweep_sort.successful_sweeps
        temp_C_L4(:,trial_id)   = spectro.coherence_LFPspike.C_mean(:,L4);
        temp_phi_L4(:,trial_id) = spectro.coherence_LFPspike.phi_mean(:,L4);
        temp_C_L23(:,trial_id)   = spectro.coherence_LFPspike.C_mean(:,L23);
        temp_phi_L23(:,trial_id) = spectro.coherence_LFPspike.phi_mean(:,L23);
    end
    temp_C_L4(:,var(temp_C_L4)==0)=NaN; temp_phi_L4(:,var(temp_phi_L4)==0)=NaN;
    temp_C_L23(:,var(temp_C_L23)==0)=NaN; temp_phi_L23(:,var(temp_phi_L23)==0)=NaN;
    meta.KO_exported.spectro.LFPspike_coherence_L4(file_idx,:)      =  nanmean(temp_C_L4,2);
    meta.KO_exported.spectro.LFPspike_phase_L4(file_idx,:)          =  nanmean(temp_phi_L4,2);
    meta.KO_exported.spectro.LFPspike_coherence_L23(file_idx,:)     =  nanmean(temp_C_L23,2);
    meta.KO_exported.spectro.LFPspike_phase_L23(file_idx,:)         =  nanmean(temp_phi_L23,2);
   
    figure; hold on
    plot(meta.KO_exported.spectro.f_LFPLFP,meta.KO_exported.spectro.LFPLFP_coherence)
    set(gca,'XScale','log');axis([5 200 -Inf 1])


end
% clear data params burst spectro spikes CSD 
% clear active_file channel_index chosen_channel file_idx no_files progbar 
% clear functions
% clear hidden

% collate individual spectra
meta.KO_exported.spectro.LFP_spectrum_L4_mean    = nanmean(meta.KO_exported.spectro.LFP_spectrum_L4);
meta.KO_exported.spectro.LFP_spectrum_L4_SEM     = nansem(meta.KO_exported.spectro.LFP_spectrum_L4);
meta.KO_exported.spectro.LFP_spectrum_L23_mean   = nanmean(meta.KO_exported.spectro.LFP_spectrum_L23);
meta.KO_exported.spectro.LFP_spectrum_L23_SEM    = nansem(meta.KO_exported.spectro.LFP_spectrum_L23);

% collate LFP-LPF coherence 
meta.KO_exported.spectro.LFPLFP_coherence_mean   = nanmean(meta.KO_exported.spectro.LFPLFP_coherence);
meta.KO_exported.spectro.LFPLFP_coherence_SEM    = nansem(meta.KO_exported.spectro.LFPLFP_coherence);
meta.KO_exported.spectro.LFPLFP_phase_mean       = nanmean(meta.KO_exported.spectro.LFPLFP_phase);
meta.KO_exported.spectro.LFPLFP_phase_SEM        = nansem(meta.KO_exported.spectro.LFPLFP_phase);
% collate spike-field coherence 
meta.KO_exported.spectro.LFPspike_coherence_L4_mean = nanmean(meta.KO_exported.spectro.LFPspike_coherence_L4);
meta.KO_exported.spectro.LFPspike_coherence_L4_SEM  = nansem(meta.KO_exported.spectro.LFPspike_coherence_L4);
meta.KO_exported.spectro.LFPspike_phase_L4_mean     = nanmean(meta.KO_exported.spectro.LFPspike_phase_L4);
meta.KO_exported.spectro.LFPspike_phase_L4_SEM      = nansem(meta.KO_exported.spectro.LFPspike_phase_L4);

meta.KO_exported.spectro.LFPspike_coherence_L23_mean = nanmean(meta.KO_exported.spectro.LFPspike_coherence_L23);
meta.KO_exported.spectro.LFPspike_coherence_L23_SEM  = nansem(meta.KO_exported.spectro.LFPspike_coherence_L23);
meta.KO_exported.spectro.LFPspike_phase_L23_mean     = nanmean(meta.KO_exported.spectro.LFPspike_phase_L23);
meta.KO_exported.spectro.LFPspike_phase_L23_SEM      = nansem(meta.KO_exported.spectro.LFPspike_phase_L23);
%% moving t-tests
% frequency spectra
meta.stats.LFP_spectral_power_L4             =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFP_spectrum_L4)',smooth_hist(meta.KO_exported.spectro.LFP_spectrum_L4)');
meta.stats.LFP_spectral_power_L23            =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFP_spectrum_L23)',smooth_hist(meta.KO_exported.spectro.LFP_spectrum_L23)');
% LFP-LFP
meta.stats.LFPLFP_coherence_magnitude        =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFPLFP_coherence)',smooth_hist(meta.KO_exported.spectro.LFPLFP_coherence)');
meta.stats.LFPLFP_coherence_phase            =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFPLFP_phase)',smooth_hist(meta.KO_exported.spectro.LFPLFP_phase)');
% LFP-spike L4
meta.stats.LFPspike_coherence_magnitude_L4   =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFPspike_coherence_L4)',smooth_hist(meta.KO_exported.spectro.LFPspike_coherence_L4)');
meta.stats.LFPspike_coherence_phase_L4       =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFPspike_phase_L4)',smooth_hist(meta.KO_exported.spectro.LFPspike_phase_L4)');
% LFP-spike L23
meta.stats.LFPspike_coherence_magnitude_L23  =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFPspike_coherence_L23)',smooth_hist(meta.KO_exported.spectro.LFPspike_coherence_L23)');
meta.stats.LFPspike_coherence_phase_L23      =   ttestseries(smooth_hist(meta.WT_exported.spectro.LFPspike_phase_L23)',smooth_hist(meta.KO_exported.spectro.LFPspike_phase_L23)');

%% Plotting 
shading_colour=[0.7 0.7 0.7];
x_lims=[5 200];
% LFP spectrum - L23
figure('Name','Spectral power: L2/3'); 
% subplot(2,1,1);
hold on
bar(meta.WT_exported.spectro.f_LFPLFP,-100.*meta.stats.LFP_spectral_power_L23.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour,'BarWidth',1)        
plot(meta.WT_exported.spectro.f_LFPLFP,-1*meta.WT_exported.spectro.LFP_spectrum_L23_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPLFP,-1*meta.KO_exported.spectro.LFP_spectrum_L23_mean,'r','LineWidth',1.5)
ciplot(-1*meta.WT_exported.spectro.LFP_spectrum_L23_mean+meta.WT_exported.spectro.LFP_spectrum_L23_SEM,...
       -1*meta.WT_exported.spectro.LFP_spectrum_L23_mean-meta.WT_exported.spectro.LFP_spectrum_L23_SEM,...
       meta.WT_exported.spectro.f_LFPLFP,'b');
ciplot(-1*meta.KO_exported.spectro.LFP_spectrum_L23_mean+meta.KO_exported.spectro.LFP_spectrum_L23_SEM,...
       -1*meta.KO_exported.spectro.LFP_spectrum_L23_mean-meta.KO_exported.spectro.LFP_spectrum_L23_SEM,...
       meta.KO_exported.spectro.f_LFPLFP,'r')   
axis([x_lims(1) 100 -50 0])
xlabel('Frequency (Hz)')   
ylabel('Power (dB)') 

%% LFP spectrum - L4
figure;
% subplot(2,1,2);
hold on
bar(meta.WT_exported.spectro.f_LFPLFP,-100.*meta.stats.LFP_spectral_power_L4.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour,'BarWidth',1)          

plot(meta.WT_exported.spectro.f_LFPLFP,-1*meta.WT_exported.spectro.LFP_spectrum_L4_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPLFP,-1*meta.KO_exported.spectro.LFP_spectrum_L4_mean,'r','LineWidth',1.5)
ciplot(-1*meta.WT_exported.spectro.LFP_spectrum_L4_mean+meta.WT_exported.spectro.LFP_spectrum_L4_SEM,...
       -1*meta.WT_exported.spectro.LFP_spectrum_L4_mean-meta.WT_exported.spectro.LFP_spectrum_L4_SEM,...
       meta.WT_exported.spectro.f_LFPLFP,'b');
ciplot(-1*meta.KO_exported.spectro.LFP_spectrum_L4_mean+meta.KO_exported.spectro.LFP_spectrum_L4_SEM,...
       -1*meta.KO_exported.spectro.LFP_spectrum_L4_mean-meta.KO_exported.spectro.LFP_spectrum_L4_SEM,...
       meta.KO_exported.spectro.f_LFPLFP,'r')   
axis([x_lims(1) x_lims(2) -60 0])
xlabel('Frequency (Hz)')   
ylabel('L4 LFP spectral power (dB)') 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP-LFP coherence
figure; subplot(2,1,1); 
hold on
bar(meta.WT_exported.spectro.f_LFPLFP,meta.stats.LFPLFP_coherence_magnitude.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour,'BarWidth',1)       
% plot(meta.WT_exported.spectro.f_LFPLFP,meta.WT_exported.spectro.LFPLFP_coherence,'color',[0.5 0.5 0.9])
% plot(meta.KO_exported.spectro.f_LFPLFP,meta.KO_exported.spectro.LFPLFP_coherence,'color',[0.9 0.5 0.5])
plot(meta.WT_exported.spectro.f_LFPLFP,meta.WT_exported.spectro.LFPLFP_coherence_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPLFP,meta.KO_exported.spectro.LFPLFP_coherence_mean,'r','LineWidth',1.5)
ciplot(meta.WT_exported.spectro.LFPLFP_coherence_mean+meta.WT_exported.spectro.LFPLFP_coherence_SEM,...
       meta.WT_exported.spectro.LFPLFP_coherence_mean-meta.WT_exported.spectro.LFPLFP_coherence_SEM,...
       meta.WT_exported.spectro.f_LFPLFP,'b');
ciplot(meta.KO_exported.spectro.LFPLFP_coherence_mean+meta.KO_exported.spectro.LFPLFP_coherence_SEM,...
       meta.KO_exported.spectro.LFPLFP_coherence_mean-meta.KO_exported.spectro.LFPLFP_coherence_SEM,...
       meta.KO_exported.spectro.f_LFPLFP,'r')   

axis([x_lims(1) x_lims(2) 0.5 1])
xlabel('Frequency (Hz)')   
ylabel('Coherence')   
% set(gca,'XScale','log')
subplot(2,1,2); hold on

bar(meta.WT_exported.spectro.f_LFPLFP,meta.stats.LFPLFP_coherence_phase.each_comparison.h*90,'EdgeColor',shading_colour,'FaceColor',shading_colour,'BarWidth',1)   
plot(meta.WT_exported.spectro.f_LFPLFP,(180/pi)*meta.WT_exported.spectro.LFPLFP_phase_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPLFP,(180/pi)*meta.KO_exported.spectro.LFPLFP_phase_mean,'r','LineWidth',1.5)
ciplot((180/pi)*(meta.WT_exported.spectro.LFPLFP_phase_mean+meta.WT_exported.spectro.LFPLFP_phase_SEM),...
       (180/pi)*(meta.WT_exported.spectro.LFPLFP_phase_mean-meta.WT_exported.spectro.LFPLFP_phase_SEM),...
       meta.WT_exported.spectro.f_LFPLFP,'b');
ciplot((180/pi)*(meta.KO_exported.spectro.LFPLFP_phase_mean+meta.KO_exported.spectro.LFPLFP_phase_SEM),...
       (180/pi)*(meta.KO_exported.spectro.LFPLFP_phase_mean-meta.KO_exported.spectro.LFPLFP_phase_SEM),...
       meta.KO_exported.spectro.f_LFPLFP,'r')   
axis([x_lims(1) x_lims(2) 0 90])
xlabel('Frequency (Hz)')   
ylabel('Phase lag (deg.)')   
% set(gca,'XScale','log')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP-spike coherence
figure('Name','spike-field coherence'); subplot(4,1,1); hold on
% L23 LFP-spike coherence
bar(meta.WT_exported.spectro.f_LFPspike,meta.stats.LFPspike_coherence_magnitude_L23.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(meta.WT_exported.spectro.f_LFPspike,meta.WT_exported.spectro.LFPspike_coherence_L23_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPspike,meta.KO_exported.spectro.LFPspike_coherence_L23_mean,'r','LineWidth',1.5)
ciplot(meta.WT_exported.spectro.LFPspike_coherence_L23_mean+meta.WT_exported.spectro.LFPspike_coherence_L23_SEM,...
       meta.WT_exported.spectro.LFPspike_coherence_L23_mean-meta.WT_exported.spectro.LFPspike_coherence_L23_SEM,...
       meta.WT_exported.spectro.f_LFPspike,'b');
ciplot(meta.KO_exported.spectro.LFPspike_coherence_L23_mean+meta.KO_exported.spectro.LFPspike_coherence_L23_SEM,...
       meta.KO_exported.spectro.LFPspike_coherence_L23_mean-meta.KO_exported.spectro.LFPspike_coherence_L23_SEM,...
       meta.KO_exported.spectro.f_LFPspike,'r')   
xlabel('Frequency (Hz)')   
ylabel('Coherence')   
axis([x_lims(1) x_lims(2) 0.5 1])
subplot(4,1,2); hold on
bar(meta.WT_exported.spectro.f_LFPspike,meta.stats.LFPspike_coherence_phase_L23.each_comparison.h*180,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(meta.WT_exported.spectro.f_LFPspike,(180/pi)*meta.WT_exported.spectro.LFPspike_phase_L23_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPspike,(180/pi)*meta.KO_exported.spectro.LFPspike_phase_L23_mean,'r','LineWidth',1.5)
ciplot((180/pi)*(meta.WT_exported.spectro.LFPspike_phase_L23_mean+meta.WT_exported.spectro.LFPspike_phase_L23_SEM),...
       (180/pi)*(meta.WT_exported.spectro.LFPspike_phase_L23_mean-meta.WT_exported.spectro.LFPspike_phase_L23_SEM),...
       meta.WT_exported.spectro.f_LFPspike,'b');
ciplot((180/pi)*(meta.KO_exported.spectro.LFPspike_phase_L23_mean+meta.KO_exported.spectro.LFPspike_phase_L23_SEM),...
       (180/pi)*(meta.KO_exported.spectro.LFPspike_phase_L23_mean-meta.KO_exported.spectro.LFPspike_phase_L23_SEM),...
       meta.KO_exported.spectro.f_LFPspike,'r')   
axis([x_lims(1) x_lims(2) 0 180])
xlabel('Frequency (Hz)')   
ylabel('Phase (deg.)')   % L23 LFP-spike coherence
subplot(4,1,3); hold on
bar(meta.WT_exported.spectro.f_LFPspike,meta.stats.LFPspike_coherence_magnitude_L4.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(meta.WT_exported.spectro.f_LFPspike,meta.WT_exported.spectro.LFPspike_coherence_L4_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPspike,meta.KO_exported.spectro.LFPspike_coherence_L4_mean,'r','LineWidth',1.5)
ciplot(meta.WT_exported.spectro.LFPspike_coherence_L4_mean+meta.WT_exported.spectro.LFPspike_coherence_L4_SEM,...
       meta.WT_exported.spectro.LFPspike_coherence_L4_mean-meta.WT_exported.spectro.LFPspike_coherence_L4_SEM,...
       meta.WT_exported.spectro.f_LFPspike,'b');
ciplot(meta.KO_exported.spectro.LFPspike_coherence_L4_mean+meta.KO_exported.spectro.LFPspike_coherence_L4_SEM,...
       meta.KO_exported.spectro.LFPspike_coherence_L4_mean-meta.KO_exported.spectro.LFPspike_coherence_L4_SEM,...
       meta.KO_exported.spectro.f_LFPspike,'r')   
axis([x_lims(1) x_lims(2) 0.5 1])
xlabel('Frequency (Hz)')   
ylabel('Coherence')     
subplot(4,1,4); hold on
bar(meta.WT_exported.spectro.f_LFPspike,meta.stats.LFPspike_coherence_phase_L4.each_comparison.h*180,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(meta.WT_exported.spectro.f_LFPspike,(180/pi)*meta.WT_exported.spectro.LFPspike_phase_L4_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.f_LFPspike,(180/pi)*meta.KO_exported.spectro.LFPspike_phase_L4_mean,'r','LineWidth',1.5)
ciplot((180/pi)*(meta.WT_exported.spectro.LFPspike_phase_L4_mean+meta.WT_exported.spectro.LFPspike_phase_L4_SEM),...
       (180/pi)*(meta.WT_exported.spectro.LFPspike_phase_L4_mean-meta.WT_exported.spectro.LFPspike_phase_L4_SEM),...
       meta.WT_exported.spectro.f_LFPspike,'b');
ciplot((180/pi)*(meta.KO_exported.spectro.LFPspike_phase_L4_mean+meta.KO_exported.spectro.LFPspike_phase_L4_SEM),...
       (180/pi)*(meta.KO_exported.spectro.LFPspike_phase_L4_mean-meta.KO_exported.spectro.LFPspike_phase_L4_SEM),...
       meta.KO_exported.spectro.f_LFPspike,'r')   
axis([x_lims(1) x_lims(2) 0 180])
xlabel('Frequency (Hz)')   
ylabel('Phase (deg.)')   
%% plot raw LFP-LFP data
for slice_id= 1:size(meta.WT_exported.spectro.LFPLFP_L4_signal,2)
    tempL4=meta.WT_exported.spectro.LFPLFP_L4_signal{slice_id};
tempL23=meta.WT_exported.spectro.LFPLFP_L23_signal{slice_id};

tempL4(:,var(tempL4)==0)=[];
tempL23(:,var(tempL23)==0)=[];

if size(tempL4,2)>10
    tempL4(:,11:size(tempL4,2))=[];
    tempL23(:,11:size(tempL23,2))=[];
end
figure('Name',meta.archive.WT_filenames{slice_id}); hold on
plot(staggerplot(tempL4,0,3),'b','LineWidth',1.2)
plot(staggerplot(tempL23,0,3),'color',[0.3 0.3 0.3],'LineWidth',1.2)
axis ([-10 201 -Inf Inf])
box off; axis off
end
for slice_id= 1:size(meta.KO_exported.spectro.LFPLFP_L4_signal,2)
    tempL4=meta.KO_exported.spectro.LFPLFP_L4_signal{slice_id};
tempL23=meta.KO_exported.spectro.LFPLFP_L23_signal{slice_id};

tempL4(:,var(tempL4)==0)=[];
tempL23(:,var(tempL23)==0)=[];

if size(tempL4,2)>10
    tempL4(:,11:size(tempL4,2))=[];
    tempL23(:,11:size(tempL23,2))=[];
end
figure('Name',meta.archive.KO_filenames{slice_id}); hold on
plot(staggerplot(tempL4,0,3),'r','LineWidth',1.2)
plot(staggerplot(tempL23,0,3),'color',[0.3 0.3 0.3],'LineWidth',1.2)
axis ([-10 201 -Inf Inf])
box off; axis off
end