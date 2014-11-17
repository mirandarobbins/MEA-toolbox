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
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    [spectro]=AD_MP(data,params);    
    toc

    meta.WT_exported.spectro.MP.average_Energy_Diff(:,:,file_idx)              =    spectro.MP.rEnergy_Diff_mean_over_trials{1};
    meta.WT_exported.spectro.MP.average_window_Energy_Diff(:,file_idx)         =    spectro.MP.rEnergy_window_mean_over_trials{1};
    meta.WT_exported.spectro.MP.average_rEnergy_cutout(:,:,file_idx)           =    spectro.MP.rEnergy_cutout_mean_over_trials{1};
    meta.WT_exported.spectro.MP.data_cutout{file_idx}                          =    spectro.MP.data_cutout{1};
    meta.WT_exported.spectro.MP.data_cutout_mean{file_idx}                     =    spectro.MP.data_cutout_mean_over_trials{1};
    meta.WT_exported.spectro.MP.f = spectro.MP.f;
    meta.WT_exported.spectro.MP.t = spectro.MP.t;
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
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    [spectro]=AD_MP(data,params);    
    toc

    meta.KO_exported.spectro.MP.average_Energy_Diff(:,:,file_idx)              =    spectro.MP.rEnergy_Diff_mean_over_trials{1};
    meta.KO_exported.spectro.MP.average_window_Energy_Diff(:,file_idx)         =    spectro.MP.rEnergy_window_mean_over_trials{1};
    meta.KO_exported.spectro.MP.average_rEnergy_cutout(:,:,file_idx)           =    spectro.MP.rEnergy_cutout_mean_over_trials{1};
    meta.KO_exported.spectro.MP.data_cutout{file_idx}                          =    spectro.MP.data_cutout{1};
    meta.KO_exported.spectro.MP.data_cutout_mean{file_idx}                     =    spectro.MP.data_cutout_mean_over_trials{1};    
    meta.KO_exported.spectro.MP.f = spectro.MP.f;
    meta.KO_exported.spectro.MP.t = spectro.MP.t;
end

clear data params burst spectro spikes CSD 
clear active_file file_idx no_files channel_index chosen_channel sweep_id temp_LFP temp_PDF this_sweep
clear functions
clear hidden
%% calculate population means... 
meta.WT_exported.spectro.MP.average_Energy_Diff(isinf(meta.WT_exported.spectro.MP.average_Energy_Diff))=NaN;
meta.WT_exported.spectro.MP.average_window_Energy_Diff(isinf(meta.WT_exported.spectro.MP.average_window_Energy_Diff))=NaN;
meta.WT_exported.spectro.MP.average_Energy_Diff_grand_mean          = nanmean(meta.WT_exported.spectro.MP.average_Energy_Diff,3);
meta.WT_exported.spectro.MP.average_Energy_Diff_grand_SEM           = nansem (meta.WT_exported.spectro.MP.average_Energy_Diff,3);
meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_mean   = nanmean(meta.WT_exported.spectro.MP.average_window_Energy_Diff,2);
meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_SEM    = nansem (meta.WT_exported.spectro.MP.average_window_Energy_Diff,2);
meta.WT_exported.spectro.MP.average_rEnergy_cutout(isinf(meta.WT_exported.spectro.MP.average_rEnergy_cutout))=NaN;
meta.WT_exported.spectro.MP.average_rEnergy_cutout_grand_mean       = nanmean(meta.WT_exported.spectro.MP.average_rEnergy_cutout,3);
meta.WT_exported.spectro.MP.average_rEnergy_cutout_grand_SEM        = nansem (meta.WT_exported.spectro.MP.average_rEnergy_cutout,3);

meta.KO_exported.spectro.MP.average_Energy_Diff(isinf(meta.KO_exported.spectro.MP.average_Energy_Diff))=NaN;
meta.KO_exported.spectro.MP.average_window_Energy_Diff(isinf(meta.KO_exported.spectro.MP.average_window_Energy_Diff))=NaN;
meta.KO_exported.spectro.MP.average_Energy_Diff_grand_mean          = nanmean(meta.KO_exported.spectro.MP.average_Energy_Diff,3);
meta.KO_exported.spectro.MP.average_Energy_Diff_grand_SEM           = nansem (meta.KO_exported.spectro.MP.average_Energy_Diff,3);
meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_mean   = nanmean(meta.KO_exported.spectro.MP.average_window_Energy_Diff,2);
meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_SEM    = nansem (meta.KO_exported.spectro.MP.average_window_Energy_Diff,2);
meta.KO_exported.spectro.MP.average_rEnergy_cutout(isinf(meta.KO_exported.spectro.MP.average_rEnergy_cutout))=NaN;
meta.KO_exported.spectro.MP.average_rEnergy_cutout_grand_mean       = nanmean(meta.KO_exported.spectro.MP.average_rEnergy_cutout,3);
meta.KO_exported.spectro.MP.average_rEnergy_cutout_grand_SEM        = nansem (meta.KO_exported.spectro.MP.average_rEnergy_cutout,3);

%% Plot window mean spectra
figure; hold on
plot(meta.WT_exported.spectro.MP.f,meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.MP.f,meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_mean,'r','LineWidth',1.5)
% plot(meta.WT_exported.spectro.MP.f,meta.WT_exported.spectro.MP.average_window_Energy_Diff,'b','LineWidth',1)
% plot(meta.KO_exported.spectro.MP.f,meta.KO_exported.spectro.MP.average_window_Energy_Diff,'r','LineWidth',1)

ciplot(meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_mean - meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_SEM,...
       meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_mean + meta.WT_exported.spectro.MP.average_window_Energy_Diff_grand_SEM,...
       meta.WT_exported.spectro.MP.f,'b')

ciplot(meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_mean - meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_SEM,...
       meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_mean + meta.KO_exported.spectro.MP.average_window_Energy_Diff_grand_SEM,...
       meta.KO_exported.spectro.MP.f,'r')
axis([0 200 -10 50])
%% Plot peak spectra
% cutout peak time spectra
ROI=95:120;
meta.WT_exported.spectro.MP.energy_at_peak          =   squeeze(nanmean(meta.WT_exported.spectro.MP.average_rEnergy_cutout(:,ROI,:),2));
meta.WT_exported.spectro.MP.energy_at_peak_mean     =   nanmean(meta.WT_exported.spectro.MP.energy_at_peak,2);
meta.WT_exported.spectro.MP.energy_at_peak_SEM      =   nansem(meta.WT_exported.spectro.MP.energy_at_peak,2);

meta.KO_exported.spectro.MP.energy_at_peak          =   squeeze(nanmean(meta.KO_exported.spectro.MP.average_rEnergy_cutout(:,ROI,:),2));
meta.KO_exported.spectro.MP.energy_at_peak_mean     =   nanmean(meta.KO_exported.spectro.MP.energy_at_peak,2);
meta.KO_exported.spectro.MP.energy_at_peak_SEM      =   nansem(meta.KO_exported.spectro.MP.energy_at_peak,2);


figure; hold on
plot(meta.WT_exported.spectro.MP.f,meta.WT_exported.spectro.MP.energy_at_peak_mean,'b','LineWidth',1.5)
plot(meta.KO_exported.spectro.MP.f,meta.KO_exported.spectro.MP.energy_at_peak_mean,'r','LineWidth',1.5)

ciplot(meta.WT_exported.spectro.MP.energy_at_peak_mean - meta.WT_exported.spectro.MP.energy_at_peak_SEM,...
       meta.WT_exported.spectro.MP.energy_at_peak_mean + meta.WT_exported.spectro.MP.energy_at_peak_SEM,...
       meta.WT_exported.spectro.MP.f,'b')

ciplot(meta.KO_exported.spectro.MP.energy_at_peak_mean - meta.KO_exported.spectro.MP.energy_at_peak_SEM,...
       meta.KO_exported.spectro.MP.energy_at_peak_mean + meta.KO_exported.spectro.MP.energy_at_peak_SEM,...
       meta.KO_exported.spectro.MP.f,'r')
axis([0 200 0 30])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')



%% mean spectrograms
clims=[0 20];
figure;
subplot(211)
imagesc(meta.WT_exported.spectro.MP.t,meta.WT_exported.spectro.MP.f,meta.WT_exported.spectro.MP.average_Energy_Diff_grand_mean); shading interp;
        axis([0 1 0 200]); 
caxis(clims)
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis xy
subplot(212)
imagesc(meta.KO_exported.spectro.MP.t,meta.KO_exported.spectro.MP.f,meta.KO_exported.spectro.MP.average_Energy_Diff_grand_mean); shading interp;
        axis([0 1 0 200]); 
caxis(clims)
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis xy
%% mean window data
figure;
% subplot(211);
hold on
% temp=cell2mat(cellfun(@transpose,meta.WT_exported.spectro.MP.data_cutout_mean,'UniformOutput',0))';
temp=cell2mat(meta.WT_exported.spectro.MP.data_cutout); temp(temp==0)=NaN;
temp_mean=nanmean(temp,2);temp_sem=nansem(temp,2);
ciplot(temp_mean-temp_sem,temp_mean+temp_sem,1:400,'b')
% plot(temp,'color',[0.6 0.6 0.6])
plot(temp_mean,'b')
% axis([0 300 -0.5 0.2])

% subplot(212); hold on
% temp=cell2mat(cellfun(@transpose,meta.KO_exported.spectro.MP.data_cutout_mean,'UniformOutput',0))';
temp=cell2mat(meta.KO_exported.spectro.MP.data_cutout); temp(temp==0)=NaN;
temp_mean=nanmean(temp,2);temp_sem=nansem(temp,2);
ciplot(temp_mean-temp_sem,temp_mean+temp_sem,1:400,'r')
% plot(temp,'color',[0.6 0.6 0.6])
plot(temp_mean,'r')
% axis([0 300 -0.5 0.2])
line([100 100],[-0.25,0],'color','k','LineStyle',':','LineWidth',1.5)
axis([0 400 -0.3 0.1])
%%
figure; hold on
temp=cell2mat(meta.KO_exported.spectro.MP.data_cutout); temp(temp==0)=NaN;
temp_mean=nanmean(temp,2);temp_sem=nansem(temp,2);
ciplot(temp_mean-temp_sem,temp_mean+temp_sem,1:300,'r')
plot(temp_mean,'r')
% plot(temp,'r')
temp=cell2mat(meta.WT_exported.spectro.MP.data_cutout); temp(temp==0)=NaN;
temp_mean=nanmean(temp,2);temp_sem=nansem(temp,2);
ciplot(temp_mean-temp_sem,temp_mean+temp_sem,1:300,'b')
plot(temp_mean,'b')
% plot(temp,'b')
% line([100 100],[-0.25,0],'color','k','LineStyle','-','LineWidth',1.5)
rectangle('Position',[95 -0.23 10 0.17],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
axis([50 200 -0.25 0.05])
axis off;


% axis([0 300 -0.5 0.2])

%% plot LFP onto window spectra 
figure;
subplot(121);
% temp=cell2mat(cellfun(@transpose,meta.WT_exported.spectro.MP.data_cutout_mean,'UniformOutput',0))';
temp=cell2mat(meta.WT_exported.spectro.MP.data_cutout); temp(temp==0)=NaN;
temp_mean=nanmean(temp,2);temp_sem=nansem(temp,2);
    clims=[0 20];
ax1 = gca;
set(gca,'LineWidth',0.00001,...
        'box','off',...
        'Color','none')
imagesc(1:400,meta.WT_exported.spectro.MP.f,...
        smooth2a(meta.WT_exported.spectro.MP.average_rEnergy_cutout_grand_mean,1,1)...
        ); shading interp;
        axis([0 400 0 200]); 
        caxis(clims)
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        axis xy
        colormap(hot)
        ax2 = axes('Position',get(ax1,'Position'),...
                   'XAxisLocation','bottom',...
                   'YAxisLocation','right',...
                   'Visible','off',...
                   'XTick',[],...
                   'YTick',[],...
                   'box','off',...
                   'Color','none',...
                   'XColor','k','YColor','k');
        linkaxes([ax1 ax2],'x');
        hold on

        
        ciplot(temp_mean-temp_sem,temp_mean+temp_sem,1:400,[1 1 1 ])
        plot(1:400,temp_mean,'Parent',ax2,'color','w','LineWidth',1.5);
        line([100 100],[-1,1],'color','k','LineStyle',':','LineWidth',1.5)
        axis([0 400 -0.3 0.2])
%         axis off


subplot(122)
% temp=cell2mat(cellfun(@transpose,meta.KO_exported.spectro.MP.data_cutout_mean,'UniformOutput',0))';
temp=cell2mat(meta.KO_exported.spectro.MP.data_cutout); temp(temp==0)=NaN;
temp_mean=nanmean(temp,2);temp_sem=nansem(temp,2);
ax3 = gca;
set(gca,'LineWidth',0.00001,...
        'box','off',...
        'Color','none')
imagesc(1:400,meta.KO_exported.spectro.MP.f,...
        smooth2a(meta.KO_exported.spectro.MP.average_rEnergy_cutout_grand_mean,1,1)...
        ); shading interp;axis([0 400 0 200]); 
        caxis(clims)
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        axis xy
%         colorbar
        ax4 = axes('Position',get(ax3,'Position'),...
                   'XAxisLocation','bottom',...
                   'YAxisLocation','right',...
                   'Visible','off',...
                   'XTick',[],...
                   'YTick',[],...
                   'box','off',...
                   'Color','none',...
                   'XColor','k','YColor','k');
%         linkaxes([    ax3 ax4],'x');
        hold on
%         axis off

        plot(1:400,temp_mean,'Parent',ax4,'color','w','LineWidth',1.5);
        ciplot(temp_mean-temp_sem,temp_mean+temp_sem,1:400,[1 1 1 ])
        line([100 100],[-1,1],'color','k','LineStyle',':','LineWidth',1.5)
        line([250 250],[-0.2,-0.1],'color','w','LineStyle','-','LineWidth',2.5)
        axis([0 400 -0.3 0.2])
%% SEM around spectrograms
figure;
subplot(221)
imagesc(meta.WT_exported.spectro.MP.t,meta.WT_exported.spectro.MP.f,meta.WT_exported.spectro.MP.average_Energy_Diff_grand_SEM); shading interp;
        axis([0 1 0 100]); caxis([ 0 50])
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis xy
subplot(222)
imagesc(meta.KO_exported.spectro.MP.t,meta.KO_exported.spectro.MP.f,meta.KO_exported.spectro.MP.average_Energy_Diff_grand_SEM); shading interp;
        axis([0 1 0 100]); caxis([ 0 50])
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis xy

subplot(223)
imagesc(1:400,meta.WT_exported.spectro.MP.f,meta.WT_exported.spectro.MP.average_rEnergy_cutout_grand_SEM); shading interp;
        axis([0 200 0 100]); 
        caxis([0 5])
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis xy

subplot(224)
imagesc(1:400,meta.KO_exported.spectro.MP.f,meta.KO_exported.spectro.MP.average_rEnergy_cutout_grand_SEM); shading interp;
        axis([0 200 0 100]); 
        caxis([0 5])
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis xy
%%
figure; hold on
% subplot(211)

plot(meta.WT_exported.spectro.MP.average_rEnergy_cutout_grand_mean,'b')
axis([0 200 0 100]); 

% subplot(212)
plot(meta.KO_exported.spectro.MP.average_rEnergy_cutout_grand_mean,'r')
        axis([0 200 0 100]); 
        
%% stats on mean power
meta.WT_exported.spectro.MP.power_downsampled_WT=[];
meta.KO_exported.spectro.MP.power_downsampled_KO=[];
%downsample mean pwoer to 5Hz bins
bin_size=10;
    power_band=bin_size:bin_size:500;
for idx = 1:numel(power_band)
    this_power_band=power_band(idx); %upper power band limit

    freq_idx = ge(meta.WT_exported.spectro.MP.f,this_power_band-10) & lt(meta.WT_exported.spectro.MP.f,this_power_band);
                   
    
    meta.WT_exported.spectro.MP.power_downsampled_WT(idx,:) = mean(meta.WT_exported.spectro.MP.energy_at_peak(freq_idx,:),1); 
    meta.KO_exported.spectro.MP.power_downsampled_KO(idx,:) = mean(meta.KO_exported.spectro.MP.energy_at_peak(freq_idx,:),1);

end

figure; 
subplot(211);hold on
plot(power_band,meta.WT_exported.spectro.MP.power_downsampled_WT,'b')
axis([0 500 0 50]); 
% set(gca,'YScale','log')

plot(power_band,meta.KO_exported.spectro.MP.power_downsampled_KO,'r')

MP_stats=ttestseries(meta.WT_exported.spectro.MP.power_downsampled_WT,meta.KO_exported.spectro.MP.power_downsampled_KO);
subplot(212)
plot(power_band,MP_stats.each_comparison.p)
axis([0 500 0 inf]); 


