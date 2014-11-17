classdef MEA_params
    % Multi-Electrode Array data set class properties
    properties
        Fs              % sampling freq (Hz)
        Nyquist         % Nyq. frequency in Hz
        selected_rep    % trial of interest
        last_sweep      % maximum no of trials (can be modified by trimming failures)
        dead_channels   % user-specified channels to ignore
        frame           % for Sgolay filtering
        degree          % for Sgolay filtering
        window_retained % vector of sample no.s to keep
        baseline_win    % [min max] for normalisation and noise analysis (in ms)
        search_win      % [min max] for response analysis (in ms)
        baseline_win_samples    % [min max] for normalisation and noise analysis (sample numbers)
        search_win_samples      % [min max] for response analysis (sample numbers)
        detection_threshold     % criterion for detection response successes
        first_stim      % time stim was applied (in ms)
        channel_index   % symetrical array (8x8) of channel locations
        flags           % flags for import scripts (struct)
        no_points       % nuber of samples per trial, per channel
        no_channels     % number of channels for analysis
        denoise         % params for noise removal (struct)
        bp              % filter design variables (struct      
    end
end