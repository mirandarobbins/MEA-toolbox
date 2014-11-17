%% bandpass filter to separate LFP and spikes

% denoise 50Hz with Chronux
params.denoise.tapers=[3 5];
params.denoise.fpass=[0 10000];
params.denoise.Fs=20000;
params.denoise.pad=-0;
p_val=0.0001/params.denoise.Fs;
for trial_id=1:params.last_sweep
    for channel_id=1:64
        clear temp; 
        temp=data.raw_data(:,channel_id,trial_id);
        
        data.filtered_lfp(:,channel_id,trial_id)   = rmlinesc(temp,params.denoise,p_val,'n',50);

    end
end
          
          
clear temp trial_id
figure;
subplot(2,1,1)
plot(squeeze(data.raw_data(:,9,:))); axis([0 20000 -0.1 0.1])

subplot(2,1,2)
plot(squeeze(data.filtered_lfp(:,9,:))); axis([0 20000 -0.1 0.1])