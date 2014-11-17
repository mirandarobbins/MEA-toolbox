classdef MEA_spikes
    % Multi-Electrode Array data set class properties
    properties
        coeffthreshold
        time_conversion
        f
        detection_thresh
        detection_rearm_time
        blanking_times
        padding_time
        spike_locs
        clusters
        waveforms_array
        spiketimes
    end
end