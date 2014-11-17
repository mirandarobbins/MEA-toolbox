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
    
    % update progressbar & load file to workspace
    waitbar(file_idx/no_files,progbar,strcat(['Loading file ',num2str(file_idx),'/',num2str(no_files)] ))
    load(active_file,'-mat');
    spikes = meta_spike_extract(data,params,4);
    spikes = meta_MUA64plot(data,params,spikes);
    toc
    close all
    L4_channel_temp=meta.archive.WT_channels_chosen(file_idx,1)+8;
    L23_channel_temp=meta.archive.WT_channels_chosen(file_idx,2);
    % LFP
    meta.exported.WT_L4_LFP{file_idx}          = squeeze(data.filtered_lfp(:,L4_channel_temp,data.sweep_sort.successful_sweeps));
    meta.exported.WT_L4_LFP_mean(:,file_idx)   = data.mean_channels(:,L4_channel_temp);
    meta.exported.WT_L4_LFP_std(:,file_idx)    = data.std_channels(:,L4_channel_temp);
    
    meta.exported.WT_L23_LFP{file_idx}         = squeeze(data.filtered_lfp(:,L23_channel_temp,data.sweep_sort.successful_sweeps));
    meta.exported.WT_L23_LFP_mean(:,file_idx)  = data.mean_channels(:,L23_channel_temp);
    meta.exported.WT_L23_LFP_std(:,file_idx)   = data.std_channels(:,L23_channel_temp);
    % monosynaptic LFP
    meta.exported.WT_L4_max_amp_mono(file_idx) = data.max_amp_mono(L4_channel_temp);
    meta.exported.WT_L4_latency_mono(file_idx) = data.latency_mono(L4_channel_temp);
    
    %trial-by-trial 
    temp_L4_latency=[];temp_L23_latency=[];temp_L4_amp=[];temp_L23_amp=[];
    
    for idx=data.sweep_sort.successful_sweeps
        temp_L4_latency(idx)  = data.burst_timing.latency{idx}(L4_channel_temp);
        temp_L23_latency(idx) = data.burst_timing.latency{idx}(L23_channel_temp);
        temp_L4_amp(idx)  = data.burst_timing.amp{idx}(L4_channel_temp);
        temp_L23_amp(idx) = data.burst_timing.amp{idx}(L23_channel_temp);
    end; clear idx
    
    meta.exported.WT_L4_latency_poly{file_idx}  = temp_L4_latency(temp_L4_latency>0);
    meta.exported.WT_L23_latency_poly{file_idx} = temp_L23_latency(temp_L4_latency>0);
    meta.exported.WT_L4_amp_poly{file_idx}      = temp_L4_amp(temp_L4_amp>0);
    meta.exported.WT_L23_amp_poly{file_idx}     = temp_L23_amp(temp_L23_amp>0); 
   
    clear temp_L4_latency temp_L23_latency temp_L4_amp temp_L23_amp
    
    %calculate L4-L2/3 transfer latency
    meta.exported.WT_transfer_latency{file_idx} = meta.exported.WT_L23_latency_poly{file_idx}-meta.exported.WT_L4_latency_poly{file_idx};
    meta.exported.WT_transfer_latency{file_idx}(meta.exported.WT_transfer_latency{file_idx}<0) = [];
    meta.exported.WT_transfer_latency_mean(file_idx) =nanmean(meta.exported.WT_transfer_latency{file_idx});
    
    meta.exported.WT_L4_MUA{file_idx}          = squeeze(data.filtered_spikes(:,L4_channel_temp,data.sweep_sort.successful_sweeps));
    meta.exported.WT_L23_MUA{file_idx}         = squeeze(data.filtered_spikes(:,L23_channel_temp,data.sweep_sort.successful_sweeps));
    
    meta.exported.WT_L4_MUA_PDF{file_idx}      = spikes.PDF.PDF_trimmed{L4_channel_temp};
    meta.exported.WT_L23_MUA_PDF{file_idx}     = spikes.PDF.PDF_trimmed{L23_channel_temp};
    meta.exported.WT_L4_MUA_PDF_mean(:,file_idx)  = spikes.PDF.mean_PDF(L4_channel_temp,:);
    meta.exported.WT_L23_MUA_PDF_mean(:,file_idx) = spikes.PDF.mean_PDF(L23_channel_temp,:);
end

clear data params burst spectro spikes CSD temp_C idx L23_channel_temp L4_channel_temp
clear active_file channel_index chosen_channel file_idx no_files progbar 
clear functions
clear hidden
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
    spikes = meta_spike_extract(data,params,4);
    spikes = meta_MUA64plot(data,params,spikes);
    toc
    close all
    L4_channel_temp=meta.archive.KO_channels_chosen(file_idx,1);
    L23_channel_temp=meta.archive.KO_channels_chosen(file_idx,2);
    % LFP
    meta.exported.KO_L4_LFP{file_idx}          = squeeze(data.filtered_lfp(:,L4_channel_temp,data.sweep_sort.successful_sweeps));
    meta.exported.KO_L4_LFP_mean(:,file_idx)   = data.mean_channels(:,L4_channel_temp);
    meta.exported.KO_L4_LFP_std(:,file_idx)    = data.std_channels(:,L4_channel_temp);
    
    meta.exported.KO_L23_LFP{file_idx}         = squeeze(data.filtered_lfp(:,L23_channel_temp,data.sweep_sort.successful_sweeps));
    meta.exported.KO_L23_LFP_mean(:,file_idx)  = data.mean_channels(:,L23_channel_temp);
    meta.exported.KO_L23_LFP_std(:,file_idx)   = data.std_channels(:,L23_channel_temp);
    % monosynaptic LFP
    meta.exported.KO_L4_max_amp_mono(file_idx) = data.max_amp_mono(L4_channel_temp);
    meta.exported.KO_L4_latency_mono(file_idx) = data.latency_mono(L4_channel_temp);
    
    %trial-by-trial 
    temp_L4_latency=[];temp_L23_latency=[];temp_L4_amp=[];temp_L23_amp=[];
    
    for idx=data.sweep_sort.successful_sweeps
        temp_L4_latency(idx)  = data.burst_timing.latency{idx}(L4_channel_temp);
        temp_L23_latency(idx) = data.burst_timing.latency{idx}(L23_channel_temp);
        temp_L4_amp(idx)  = data.burst_timing.amp{idx}(L4_channel_temp);
        temp_L23_amp(idx) = data.burst_timing.amp{idx}(L23_channel_temp);
    end; clear idx
    
    meta.exported.KO_L4_latency_poly{file_idx}  = temp_L4_latency(temp_L4_latency>0);
    meta.exported.KO_L23_latency_poly{file_idx} = temp_L23_latency(temp_L4_latency>0);
    meta.exported.KO_L4_amp_poly{file_idx}      = temp_L4_amp(temp_L4_amp>0);
    meta.exported.KO_L23_amp_poly{file_idx}     = temp_L23_amp(temp_L23_amp>0); 
   
    clear temp_L4_latency temp_L23_latency temp_L4_amp temp_L23_amp
    
    %calculate L4-L2/3 transfer latency
    meta.exported.KO_transfer_latency{file_idx} = meta.exported.KO_L23_latency_poly{file_idx}-meta.exported.KO_L4_latency_poly{file_idx};
    meta.exported.KO_transfer_latency{file_idx}(meta.exported.KO_transfer_latency{file_idx}<0) = [];
    meta.exported.KO_transfer_latency_mean(file_idx) =nanmean(meta.exported.KO_transfer_latency{file_idx});
    
    meta.exported.KO_L4_MUA{file_idx}          = squeeze(data.filtered_spikes(:,L4_channel_temp,data.sweep_sort.successful_sweeps));
    meta.exported.KO_L23_MUA{file_idx}         = squeeze(data.filtered_spikes(:,L23_channel_temp,data.sweep_sort.successful_sweeps));
    
    meta.exported.KO_L4_MUA_PDF{file_idx}      = spikes.PDF.PDF_trimmed{L4_channel_temp};
    meta.exported.KO_L23_MUA_PDF{file_idx}     = spikes.PDF.PDF_trimmed{L23_channel_temp};
    meta.exported.KO_L4_MUA_PDF_mean(:,file_idx)  = spikes.PDF.mean_PDF(L4_channel_temp,:);
    meta.exported.KO_L23_MUA_PDF_mean(:,file_idx) = spikes.PDF.mean_PDF(L23_channel_temp,:);
end

clear data params burst spectro spikes CSD temp_C idx L23_channel_temp L4_channel_temp
clear active_file channel_index chosen_channel file_idx no_files progbar 
clear functions
clear hidden

%% averages
% meta.exported.KO_L4_LFP_mean(:,2)=[];
% meta.exported.KO_L23_LFP_mean(:,2)=[];
% meta.exported.KO_L4_MUA_PDF_mean(:,2)=[];
% meta.exported.KO_L23_MUA_PDF_mean(:,2)=[];
% LFP
meta.exported.averages.WT_L4_LFP_mean  = nanmean(meta.exported.WT_L4_LFP_mean,2);
meta.exported.averages.WT_L4_LFP_SEM   = nansem(meta.exported.WT_L4_LFP_mean,2);
meta.exported.averages.WT_L23_LFP_mean = nanmean(meta.exported.WT_L23_LFP_mean,2);
meta.exported.averages.WT_L23_LFP_SEM  = nansem(meta.exported.WT_L23_LFP_mean,2);

meta.exported.averages.KO_L4_LFP_mean  = nanmean(meta.exported.KO_L4_LFP_mean,2);
meta.exported.averages.KO_L4_LFP_SEM   = nansem(meta.exported.KO_L4_LFP_mean,2);
meta.exported.averages.KO_L23_LFP_mean = nanmean(meta.exported.KO_L23_LFP_mean,2);
meta.exported.averages.KO_L23_LFP_SEM  = nansem(meta.exported.KO_L23_LFP_mean,2);

% MUA
meta.exported.averages.WT_L4_MUA_PDF_mean  = nanmean(meta.exported.WT_L4_MUA_PDF_mean,2);
meta.exported.averages.WT_L4_MUA_PDF_SEM   = nansem(meta.exported.WT_L4_MUA_PDF_mean,2);
meta.exported.averages.WT_L23_MUA_PDF_mean = nanmean(meta.exported.WT_L23_MUA_PDF_mean,2);
meta.exported.averages.WT_L23_MUA_PDF_SEM  = nansem(meta.exported.WT_L23_MUA_PDF_mean,2);

meta.exported.averages.KO_L4_MUA_PDF_mean  = nanmean(meta.exported.KO_L4_MUA_PDF_mean,2);
meta.exported.averages.KO_L4_MUA_PDF_SEM   = nansem(meta.exported.KO_L4_MUA_PDF_mean,2);
meta.exported.averages.KO_L23_MUA_PDF_mean = nanmean(meta.exported.KO_L23_MUA_PDF_mean,2);
meta.exported.averages.KO_L23_MUA_PDF_SEM  = nansem(meta.exported.KO_L23_MUA_PDF_mean,2);

%% moving t-tests
% frequency spectra
meta.exported.stats.MUA_L4   =   ttestseries(smooth_hist(meta.exported.WT_L4_MUA_PDF_mean),smooth_hist(meta.exported.KO_L4_MUA_PDF_mean));
meta.exported.stats.MUA_L23  =   ttestseries(smooth_hist(meta.exported.WT_L23_MUA_PDF_mean),smooth_hist(meta.exported.KO_L23_MUA_PDF_mean));

meta.exported.stats.MUA_4_vs_L23_KO  =   ttestseries(smooth_hist(meta.exported.KO_L4_MUA_PDF_mean),smooth_hist(meta.exported.KO_L23_MUA_PDF_mean));
%% plot
% Plot LFP
tb=1:20000;tb=tb/20;
figure; 
subplot(1,2,1); hold on
plot(tb,meta.exported.WT_L23_LFP_mean,'b')
plot(tb,meta.exported.WT_L4_LFP_mean-0.4,'b')
axis off; axis([0 1000 -Inf Inf])

% Plot LFP
subplot(1,2,2); hold on
plot(tb,meta.exported.KO_L23_LFP_mean,'r')
plot(tb,meta.exported.KO_L4_LFP_mean-0.4,'r')
axis off; axis([0 1000 -Inf Inf])
figure; 
subplot(2,1,1); hold on
plot(tb,meta.exported.averages.WT_L4_LFP_mean,'b')
plot(tb,meta.exported.averages.KO_L4_LFP_mean,'r')

ciplot(meta.exported.averages.WT_L4_LFP_mean+meta.exported.averages.WT_L4_LFP_SEM,...
       meta.exported.averages.WT_L4_LFP_mean-meta.exported.averages.WT_L4_LFP_SEM,...
       tb,'b')

ciplot(meta.exported.averages.KO_L4_LFP_mean+meta.exported.averages.KO_L4_LFP_SEM,...
       meta.exported.averages.KO_L4_LFP_mean-meta.exported.averages.KO_L4_LFP_SEM,...
       tb,'r')

subplot(2,1,2); hold on
plot(tb,meta.exported.averages.WT_L23_LFP_mean,'b')
plot(tb,meta.exported.averages.KO_L23_LFP_mean,'r')

ciplot(meta.exported.averages.WT_L23_LFP_mean+meta.exported.averages.WT_L23_LFP_SEM,...
       meta.exported.averages.WT_L23_LFP_mean-meta.exported.averages.WT_L23_LFP_SEM,...
       tb,'b')

ciplot(meta.exported.averages.KO_L23_LFP_mean+meta.exported.averages.KO_L23_LFP_SEM,...
       meta.exported.averages.KO_L23_LFP_mean-meta.exported.averages.KO_L23_LFP_SEM,...
       tb,'r')   
axis([0 500 -Inf Inf])
%% plot MUA (1) each area overlaid -half wave rectify
tb=1:20000;tb=tb/20;
alpha=1;

figure; 
subplot(2,1,1); hold on
for file_idx=1:size(meta.exported.WT_L4_MUA)
        temp=abs(meta.exported.WT_L4_MUA{file_idx});
        temp([1980:2060,... 
                                  2380:2460,... 
                                  2780:2860,...
                                  3180:3260,...
                                  3580:3660],:)=NaN;
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,temp,'edgecolor','b','linewidth',1,'edgealpha',1);

%     plot(abs(meta.exported.WT_L4_MUA{file_idx}),'b')
%     line([0, 20000],...
%          [8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000))))])
end;    
for file_idx=1:size(meta.exported.WT_L4_MUA)
        temp=abs(meta.exported.WT_L23_MUA{file_idx});
                temp([1980:2060,... 
                                  2380:2460,... 
                                  2780:2860,...
                                  3180:3260,...
                                  3580:3660],:)=NaN;
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,temp,'edgecolor',[0.6 0.6 0.6],'linewidth',1,'edgealpha',alpha);
%     plot(abs(meta.exported.WT_L23_MUA{file_idx}),'color',[0.5 0.5 0.5])

end;    
axis([0 550 0 .08])% axis([0 500 0 0.08])

axis off
subplot(2,1,2); hold on
for file_idx=1:size(meta.exported.KO_L4_MUA)
        temp=abs(meta.exported.KO_L4_MUA{file_idx});
                        temp([1980:2060,... 
                                  2380:2460,... 
                                  2780:2860,...
                                  3180:3260,...
                                  3580:3660],:)=NaN;
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,temp,'edgecolor','r','linewidth',1,'edgealpha',1);
%         line([0, 20000],[8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000))))],'color','r')


%     plot(abs(meta.exported.WT_L4_MUA{file_idx}),'b')
%     line([0, 20000],...
%          [8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000))))])
end;    
for file_idx=1:size(meta.exported.KO_L4_MUA)
        temp=abs(meta.exported.KO_L23_MUA{file_idx});
                        temp([1980:2060,... 
                                  2380:2460,... 
                                  2780:2860,...
                                  3180:3260,...
                                  3580:3660],:)=NaN;
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,temp,'edgecolor',[0.6 0.6 0.6],'linewidth',1,'edgealpha',alpha);
%     plot(abs(meta.exported.WT_L23_MUA{file_idx}),'color',[0.5 0.5 0.5])

end;    

line([450 500],[0.04 0.04],'LineWidth',1.5,'color','k')
line([450 450],[0.04 0.06],'LineWidth',1.5,'color','k')
axis([0 550 0 .08])
% axis([0 500 0 0.08])

axis off
%% plot MUA (1) each genotype overlaid -half wave rectify
tb=1:20000;
alpha=0.8;
figure; 
subplot(2,1,1); hold on


for file_idx=1:size(meta.exported.WT_L4_MUA)
%     plot(abs(meta.exported.WT_L4_MUA{file_idx}),'b')
    temp=abs(meta.exported.WT_L4_MUA{file_idx});
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,temp,'edgecolor','b','linewidth',1,'edgealpha',1);
%     line([0, 20000],[8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000))))])
end;
for file_idx=1:size(meta.exported.KO_L4_MUA)
%     plot(abs(meta.exported.KO_L4_MUA{file_idx}),'b')
    temp=abs(meta.exported.KO_L4_MUA{file_idx});
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,-1*temp,'edgecolor','r','linewidth',1,'edgealpha',alpha);
%     line([0, 20000],[8*mean(std(abs(meta.exported.KO_L4_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.KO_L4_MUA{file_idx}(1:2000))))])
end;
axis([0 Inf -0.1 .1])    
subplot(2,1,2); hold on
for file_idx=1:size(meta.exported.WT_L23_MUA)
%     plot(abs(meta.exported.WT_L23_MUA{file_idx}),'b')
    temp=abs(meta.exported.WT_L23_MUA{file_idx});
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,temp,'edgecolor','b','linewidth',1,'edgealpha',1);
%     line([0, 20000],[8*mean(std(abs(meta.exported.WT_L23_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.WT_L23_MUA{file_idx}(1:2000))))])
end;
for file_idx=1:size(meta.exported.KO_L23_MUA)
%     plot(abs(meta.exported.KO_L23_MUA{file_idx}),'b')
    temp=abs(meta.exported.KO_L23_MUA{file_idx});
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
    patchline(tb_temp,-1*temp,'edgecolor','r','linewidth',1,'edgealpha',alpha);
%     line([0, 20000],[8*mean(std(abs(meta.exported.KO_L23_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.KO_L23_MUA{file_idx}(1:2000))))])
end;

axis([0 Inf -.1 .1])    
%% plot MUA (2) each genotype overlaid
figure; 
subplot(2,1,1); hold on

for file_idx=1:size(meta.exported.KO_L4_MUA)
    plot(meta.exported.KO_L4_MUA{file_idx}-.05,'r')
end;   
for file_idx=1:size(meta.exported.WT_L4_MUA)
    plot(meta.exported.WT_L4_MUA{file_idx}+.05,'b')
end;axis([0 Inf -.15 .15])    

subplot(2,1,2); hold on
for file_idx=1:size(meta.exported.KO_L4_MUA)
    plot(meta.exported.KO_L23_MUA{file_idx}-.05,'r')
end;   
for file_idx=1:size(meta.exported.WT_L4_MUA)
    plot(meta.exported.WT_L23_MUA{file_idx}+.05,'b')
end;    axis([0 Inf -.1 .1])

%% Plot MUA_PDF
meta.exported.averages.WT_L4_MUA_PDF_mean(isnan(meta.exported.averages.WT_L4_MUA_PDF_mean))=0;
meta.exported.averages.WT_L4_MUA_PDF_SEM(isnan(meta.exported.averages.WT_L4_MUA_PDF_SEM))=0;
meta.exported.averages.WT_L23_MUA_PDF_mean(isnan(meta.exported.averages.WT_L23_MUA_PDF_mean))=0;
meta.exported.averages.WT_L23_MUA_PDF_SEM(isnan(meta.exported.averages.WT_L23_MUA_PDF_SEM))=0;

meta.exported.averages.KO_L4_MUA_PDF_mean(isnan(meta.exported.averages.KO_L4_MUA_PDF_mean))=0;
meta.exported.averages.KO_L4_MUA_PDF_SEM(isnan(meta.exported.averages.KO_L4_MUA_PDF_SEM))=0;
meta.exported.averages.KO_L23_MUA_PDF_mean(isnan(meta.exported.averages.KO_L23_MUA_PDF_mean))=0;
meta.exported.averages.KO_L23_MUA_PDF_SEM(isnan(meta.exported.averages.KO_L23_MUA_PDF_SEM))=0;

shading_colour=[0.7 0.7 0.7];


tb=1:1000;
figure; 
subplot(2,1,1); hold on
bar(1:1000,meta.exported.stats.MUA_L4.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   

plot(tb,meta.exported.averages.WT_L4_MUA_PDF_mean,'b')
plot(tb,meta.exported.averages.KO_L4_MUA_PDF_mean,'r')

ciplot(meta.exported.averages.WT_L4_MUA_PDF_mean+meta.exported.averages.WT_L4_MUA_PDF_SEM,...
       meta.exported.averages.WT_L4_MUA_PDF_mean-meta.exported.averages.WT_L4_MUA_PDF_SEM,...
       tb,'b')

ciplot(meta.exported.averages.KO_L4_MUA_PDF_mean+meta.exported.averages.KO_L4_MUA_PDF_SEM,...
       meta.exported.averages.KO_L4_MUA_PDF_mean-meta.exported.averages.KO_L4_MUA_PDF_SEM,...
       tb,'r')
axis([0 500 0 1.8])

subplot(2,1,2); hold on
bar(1:1000,meta.exported.stats.MUA_L23.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(tb,meta.exported.averages.WT_L23_MUA_PDF_mean,'b')
plot(tb,meta.exported.averages.KO_L23_MUA_PDF_mean,'r')

ciplot(meta.exported.averages.WT_L23_MUA_PDF_mean+meta.exported.averages.WT_L23_MUA_PDF_SEM,...
       meta.exported.averages.WT_L23_MUA_PDF_mean-meta.exported.averages.WT_L23_MUA_PDF_SEM,...
       tb,'b')

ciplot(meta.exported.averages.KO_L23_MUA_PDF_mean+meta.exported.averages.KO_L23_MUA_PDF_SEM,...
       meta.exported.averages.KO_L23_MUA_PDF_mean-meta.exported.averages.KO_L23_MUA_PDF_SEM,...
       tb,'r')      
axis([0 500 0 1.8])

%% Plot MUA_PDF
meta.exported.averages.WT_L4_MUA_PDF_mean(isnan(meta.exported.averages.WT_L4_MUA_PDF_mean))=0;
meta.exported.averages.WT_L4_MUA_PDF_SEM(isnan(meta.exported.averages.WT_L4_MUA_PDF_SEM))=0;
meta.exported.averages.WT_L23_MUA_PDF_mean(isnan(meta.exported.averages.WT_L23_MUA_PDF_mean))=0;
meta.exported.averages.WT_L23_MUA_PDF_SEM(isnan(meta.exported.averages.WT_L23_MUA_PDF_SEM))=0;

meta.exported.averages.KO_L4_MUA_PDF_mean(isnan(meta.exported.averages.KO_L4_MUA_PDF_mean))=0;
meta.exported.averages.KO_L4_MUA_PDF_SEM(isnan(meta.exported.averages.KO_L4_MUA_PDF_SEM))=0;
meta.exported.averages.KO_L23_MUA_PDF_mean(isnan(meta.exported.averages.KO_L23_MUA_PDF_mean))=0;
meta.exported.averages.KO_L23_MUA_PDF_SEM(isnan(meta.exported.averages.KO_L23_MUA_PDF_SEM))=0;

%% plot by genotype
shading_colour=[0.9 0.9 0.9];
tb=(1:1000)-100;
figure; 
% subplot(2,1,1); hold on
% plot(tb,meta.exported.averages.WT_L4_MUA_PDF_mean,'b','LineWidth',1.5)
% plot(tb,meta.exported.averages.WT_L23_MUA_PDF_mean,'color',[0.3 0.3 0.3],'LineWidth',1.5)
% 
% ciplot(meta.exported.averages.WT_L23_MUA_PDF_mean+meta.exported.averages.WT_L23_MUA_PDF_SEM,...
%        meta.exported.averages.WT_L23_MUA_PDF_mean-meta.exported.averages.WT_L23_MUA_PDF_SEM,...
%        tb,[0.3 0.3 0.3])      
%    
% ciplot(meta.exported.averages.WT_L4_MUA_PDF_mean+meta.exported.averages.WT_L4_MUA_PDF_SEM,...
%        meta.exported.averages.WT_L4_MUA_PDF_mean-meta.exported.averages.WT_L4_MUA_PDF_SEM,...
%        tb,'b')
% axis([100 300 0 0.5]); axis off
% 
% subplot(2,1,2); 
hold on
bar(tb,2*meta.exported.stats.MUA_4_vs_L23_KO.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   


plot(tb,meta.exported.averages.KO_L4_MUA_PDF_mean,'r','LineWidth',1.5)
plot(tb,meta.exported.averages.KO_L23_MUA_PDF_mean,'color',[0.3 0.3 0.3],'LineWidth',1.5)

ciplot(meta.exported.averages.KO_L23_MUA_PDF_mean+meta.exported.averages.KO_L23_MUA_PDF_SEM,...
       meta.exported.averages.KO_L23_MUA_PDF_mean-meta.exported.averages.KO_L23_MUA_PDF_SEM,...
       tb,[0.3 0.3 0.3])      
   ciplot(meta.exported.averages.KO_L4_MUA_PDF_mean+meta.exported.averages.KO_L4_MUA_PDF_SEM,...
       meta.exported.averages.KO_L4_MUA_PDF_mean-meta.exported.averages.KO_L4_MUA_PDF_SEM,...
       tb,'r')
axis([100 200 0 1]); axis off
%% plot by layer
shading_colour=[0.9 0.9 0.9];

tb=(1:1000)-100;
figure; 
subplot(2,1,1); hold on
bar(tb,2*meta.exported.stats.MUA_L23.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(tb,meta.exported.averages.WT_L23_MUA_PDF_mean,'color','b','LineWidth',1.5)
plot(tb,meta.exported.averages.KO_L23_MUA_PDF_mean,'color','r','LineWidth',1.5)
ciplot(meta.exported.averages.WT_L23_MUA_PDF_mean+meta.exported.averages.WT_L23_MUA_PDF_SEM,...
       meta.exported.averages.WT_L23_MUA_PDF_mean-meta.exported.averages.WT_L23_MUA_PDF_SEM,...
       tb,'b')     
ciplot(meta.exported.averages.KO_L23_MUA_PDF_mean+meta.exported.averages.KO_L23_MUA_PDF_SEM,...
       meta.exported.averages.KO_L23_MUA_PDF_mean-meta.exported.averages.KO_L23_MUA_PDF_SEM,...
       tb,'r')      
ylabel('p(spike/5ms)')
xlabel('Post-stimulus time (ms)')
axis([-10 300 0 1.2]); axis off

subplot(2,1,2); hold on
bar(tb,2*meta.exported.stats.MUA_L4.each_comparison.h,'EdgeColor',shading_colour,'FaceColor',shading_colour)   
plot(tb,meta.exported.averages.WT_L4_MUA_PDF_mean,'b','LineWidth',1.5)
plot(tb,meta.exported.averages.KO_L4_MUA_PDF_mean,'r','LineWidth',1.5)
   
ciplot(meta.exported.averages.WT_L4_MUA_PDF_mean+meta.exported.averages.WT_L4_MUA_PDF_SEM,...
       meta.exported.averages.WT_L4_MUA_PDF_mean-meta.exported.averages.WT_L4_MUA_PDF_SEM,...
       tb,'b')
ciplot(meta.exported.averages.KO_L4_MUA_PDF_mean+meta.exported.averages.KO_L4_MUA_PDF_SEM,...
       meta.exported.averages.KO_L4_MUA_PDF_mean-meta.exported.averages.KO_L4_MUA_PDF_SEM,...
       tb,'r')
ylabel('p(spike/5ms)')
xlabel('Post-stimulus time (ms)')
axis([-10 300 0 1.2]); %axis off
%% plot MUA (1) each area overlaid - half wave rectify
tb=1:20000;tb=tb/20;
alpha=1;

figure; 
% subplot(2,1,1); 
hold on
for file_idx=1:size(meta.exported.KO_L4_MUA,2)
        temp=abs(meta.exported.KO_L4_MUA{file_idx});
                        temp([1980:2060,... 
                                  2380:2460,... 
                                  2780:2860,...
                                  3180:3260,...
                                  3580:3660],:)=NaN;
    temp=nanmean(temp,2);              
    KO_temp(file_idx,:)=temp;
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
%     plot(tb_temp,temp,'r')
%     patchline(tb_temp,temp,'edgecolor','r','linewidth',1,'edgealpha',1);

end;    
for file_idx=1:size(meta.exported.WT_L4_MUA,2)
        temp=abs(meta.exported.WT_L4_MUA{file_idx});
        temp([1980:2060,... 
                                  2380:2460,... 
                                  2780:2860,...
                                  3180:3260,...
                                  3580:3660],:)=NaN;
    temp=nanmean(temp,2);   
    WT_temp(file_idx,:)=temp;
    tb_temp=repmat(tb,size(temp,2),1);tb_temp=tb_temp';
%     plot(tb_temp,temp,'b')
%     patchline(tb_temp,temp,'edgecolor','b','linewidth',1,'edgealpha',1);

%     plot(abs(meta.exported.WT_L4_MUA{file_idx}),'b')
%     line([0, 20000],...
%          [8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000)))), 8*mean(std(abs(meta.exported.WT_L4_MUA{file_idx}(1:2000))))])
end;    

    plot(tb_temp,nanmean(WT_temp),'b')
    plot(tb_temp,nanmean(KO_temp),'r')

line([450 500],[0.04 0.04],'LineWidth',1.5,'color','k')
line([450 450],[0.04 0.06],'LineWidth',1.5,'color','k')
% axis([0 550 -0.08 .08])
% axis([0 500 0 0.08])
axis off