function JointHist= JPSTH64(spikes)
% REMOVE Chronux from path first!


JointHist.coincidence_window_width =5;
t1=tic;
progbar = waitbar(0,'Initializing...',...
            'name','JPSTH analysis progress')%,...
%             'position',[640 600 275 50]);


run_id=1;
for idx=1:numel(spikes.PDF.timestamps)
    for idy=1:numel(spikes.PDF.timestamps)
        
        % update progressbar
        waitbar(run_id/(numel(spikes.PDF.timestamps)^2),progbar,...
        strcat(['Calculating JPSTH for pairing channel ' num2str(idx) ' vs. ' num2str(idy) , ' (', num2str(run_id)  '/' num2str((numel(spikes.PDF.timestamps)^2)) ') ...'] ))

        JointHist.jpsth_array{idx,idy}=jpsth(spikes.PDF.timestamps{idx},...
                                             spikes.PDF.timestamps{idy},...
                                             JointHist.coincidence_window_width);
        
        run_id=run_id+1;
    end
end

waitbar(run_id/(numel(spikes.PDF.timestamps)^2),progbar,  ['Finished: Processed ',num2str((numel(spikes.PDF.timestamps)^2)),' calcultiona in ', num2str(toc(t1)), ' seconds.']);


%% plot
chan_1=18;
chan_2=19;
temp=jpsth(spikes.PDF.timestamps{chan_1},spikes.PDF.timestamps{chan_2},50, chan_1, chan_2);
temp=jpsth(smooth_hist(spikes.PDF.timestamps{chan_1}),smooth_hist(spikes.PDF.timestamps{chan_2}),10, chan_1, chan_2);
plotJPSTH(temp,chan_1,chan_2,1:1000,spikes.PDF.binwidth);
%%
figure; hold on
tb= -50:50
h4=bar(tb, temp.xcorr_hist, 'histc');
h5=stairs(tb, temp.sig_low,'r')
h6=stairs(tb, temp.sig_high,'r')
