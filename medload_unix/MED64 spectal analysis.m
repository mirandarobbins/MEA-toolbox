%% OPTION 1 - single trace example

figure;
freqs=0:10:100;

t1=tic;
progbar = waitbar(0,'Initializing...',...
            'name','spectrogram progress')%,...
%             'position',[640 600 275 50]);

for channel_id=1:params.no_channels    
    
    % update progress bar
    waitbar(channel_id/params.no_channels,progbar,...
    strcat(['Analysing channel ' num2str(channel_id)  '/' num2str(params.no_channels) '...'] ))

        temp=data.filtered_lfp(:,channel_id,params.selected_rep);
        
        subaxis(8,8,channel_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.1);

        spectrogram(temp(:),hamming(1000),900,freqs,params.Fs,'yaxis');

set(gca,'xtick',[],...
        'ytick',[],...
        'clim',[-100 -20])
colormap(hot)
clear temp
%     drawnow ;pause (0.0000001)

end
% finalise progress bar
progbar=waitbar(channel_id/params.no_channels,  ['Finished: Processed ',num2str(params.no_channels),' channels in ', num2str(toc(t1)), ' seconds.']);
clear t1
%% OPTION 2 - loop through all channels
    no_channels = size(all_channels,2);
    this_trial=1;

freqs=0:1:100;
Fs=2E4;

t1=tic;
progbar = waitbar(0,'Initializing...',...
            'name','spectrogram progress')%,...
%             'position',[640 600 275 50]);


for  plot_id=1:no_channels
    
    % update progress bar
    waitbar(plot_id/no_channels,progbar,...
    strcat(['Analysing channel ' num2str(plot_id)  '/' num2str(no_channels) '...'] ))

        temp=all_channels(:,plot_id,this_trial);

    % Run the spectrogram
    [spectro.S(:,:,plot_id),...
     spectro.F(:,:,plot_id),...
     spectro.T(:,:,plot_id),...
     spectro.P(:,:,plot_id)]   =...
spectrogram(temp(:),hamming(1000),999,freqs,Fs,'yaxis');
end
   
% finalise progress bar
waitbar(plot_id/no_channels,  ['Finished: Processed ',num2str(no_channels),' channels in ', num2str(toc(t1)), ' seconds.']);
clear t1
%      spectro.mean=mean(10*log10(abs(spectro.P)),3); 
%      spectro.SD=std(10*log10(abs(spectro.P)),0,3);
%%
spectro.mean_smooth=smooth2a(spectro.mean,10, 100); %2D smoothing
spectro.SD_smooth=smooth2a(spectro.SD,10, 100); %2D smoothing
spectro.tb=(0:size(data,1)-1)./Fs;
x_temp=repmat(spectro.tb,no_sweeps,1);x_temp=x_temp';
spectro.freqs=freqs;
% spectro.data=data;
 figure;
subplot(311)
    plot(x_temp,data); 
    axis([0 1 -Inf Inf ])
    xlabel('time (ms)')
    ylabel('membrane potential (mV)')
subplot(312)

hndl = surf(spectro.T(:,:,1),...
            spectro.F(:,:,1),...
            spectro.mean,...
                                'EdgeColor','none');

            axis xy; axis tight;
            colormap(hot);
            title('Power Spectral Density estimate (dBuV/Hz)')
            ylabel('Frequency (Hz)');
            xlabel('Time (s)');
            view(0,90); % AZ = 0, EL = 90 is directly overhead and the default 2-D view.
subplot(313)
hndl = surf(spectro.T(:,:,1),...
            spectro.F(:,:,1),...
            spectro.SD_smooth,...
                                'EdgeColor','none');

            axis xy; axis tight;
            colormap(hot);
            title('Standard Deviation per frequency bin')
            ylabel('Frequency (Hz)');
            xlabel('Time (s)');
            view(0,90); % AZ = 0, EL = 90 is directly overhead and the default 2-D view.

    clear x_temp hndl sweep_id freqs no_sweeps data Fs progbar