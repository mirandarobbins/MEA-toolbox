% Cell_1=data.WC;
% Cell_2=data.WC;
%%
channel_to_plot=Cell_1.closestMEA;
figure; hold on
plot(staggerplot(Cell_1.chan_aligned(1:8000,:),8100,0),'color', [1 0.5 0.5])
plot(staggerplot(Cell_2.chan_aligned(1:8000,:),8100,0),'color', [0.5 0.5 1])
plot(100*staggerplot(squeeze(downsample(data.filtered_lfp(1:16000,channel_to_plot,:),2)),8100,0)-50,'color', [0.5 1 0.5])

%% bandpass filter to separate LFP and spikes

params.bp.gamma=[30/params.Nyquist 80/params.Nyquist];
    
[params.bp.B_gamma,  params.bp.A_gamma]   = butter(1,params.bp.gamma,'bandpass');

%    figure; fvtool(params.bp.B_gamma,params.bp.B_gamma,'Fs',params.Fs)  

% denoise 50Hz with Chronux
params.denoise.tapers=[3 1];
params.denoise.fpass=[30 80];
params.denoise.Fs=20000;

for trial_id=1:params.last_sweep
    for channel_id=1:64
        clear temp; 
        temp=data.raw_data(:,channel_id,trial_id);
        data.filtered_gamma(:,channel_id,trial_id)   = filter(params.bp.B_gamma,params.bp.A_gamma,temp);
        data.filtered_gamma(:,channel_id,trial_id)   = rmlinesc(data.filtered_gamma(:,channel_id,trial_id),params.denoise,0.001/params.denoise.Fs,'n',50);
        data.filtered_lfp(:,channel_id,trial_id)     = rmlinesc(data.filtered_lfp(:,channel_id,trial_id),params.denoise,0.001/params.denoise.Fs,'n',50);

    end
end
          
          
clear temp trial_id
% figure; plot(squeeze(data.filtered_gamma(:,Cell_1.closestMEA,:)));

%% plot IC
figure;
no_plots=size(Cell_1.chan_aligned,2);
for trial_id=1:no_plots
    subaxis(no_plots,1,trial_id, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0.01);hold on
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off; %axis off
    
    plot(params.tb_WC*1000,Cell_1.chan_aligned(:,trial_id),'color', [1 0.5 0.5])
    plot(params.tb_WC*1000,Cell_2.chan_aligned(:,trial_id),'color', [0.5 0.5 1])
    plot(data.tb,data.filtered_gamma(:,Cell_1.closestMEA,trial_id)*1000,'color', [0.5 0.8 0.5])
    
    axis([ 0 800 -50 50])
end

%% plot VC
figure;
no_plots=size(Cell_1.chan_aligned,2);
for trial_id=1:no_plots
    subaxis(no_plots,1,trial_id, 'Spacing', 0.0, 'Padding', 0, 'Margin', 0.01);hold on
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off; %axis off
    
    plot(params.tb_WC*1000,sgolayfilt(Cell_1.chan_aligned(:,trial_id),0,15),'color', [1 0.5 0.5])
    plot(params.tb_WC*1000,sgolayfilt(Cell_2.chan_aligned(:,trial_id),0,15),'color', [0.5 0.5 1])
    plot(data.tb,data.filtered_lfp(:,Cell_1.closestMEA,trial_id)*500,'color', [0.5 0.8 0.5])
    
    axis([ 0 800 -100 100])
end
