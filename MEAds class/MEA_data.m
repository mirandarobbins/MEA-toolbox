classdef MEA_data
    % Multi-Electrode Array data set class properties
    properties
        this_file       % file name
        tb              % timebase (in ms)
        raw_data        % archived as imported
        sweep_sort
        trimmed         % as raw_data but with response failures detected and removed
        filtered_lfp    % Low pass filtered version of trimmed
        filtered_spikes % High passed filterd version of trimmed
        mean_channels   % Experiment average (mean) LFP of detected successes
        std_channels    % Experiment average std LFP of detected successes
        max_amp         % Peak amplitude of LFP in detection window
        latency         % post-stimulus latency to peak LFP in detection window
    end
end