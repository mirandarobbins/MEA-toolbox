channel_to_plot=spectro.channeltoanalyse
for trial_id=data.sweep_sort.successful_sweeps
temp_x(trial_id)=-1*data.burst_timing.amp{trial_id}(channel_to_plot);
temp_y20(trial_id)=100+spectro.average_window_power(12,trial_id);
temp_y50(trial_id)=100+spectro.average_window_power(22,trial_id);
temp_y100(trial_id)=100+spectro.average_window_power(43,trial_id);
temp_y200(trial_id)=100+spectro.average_window_power(60,trial_id);
end
hold on
scatter(temp_x, temp_y20,'og')
scatter(temp_x, temp_y50,'ob')
scatter(temp_x, temp_y100,'or')
scatter(temp_x, temp_y200,'ok')