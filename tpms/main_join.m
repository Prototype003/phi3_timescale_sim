%% Description

%{

Joins split results

%}

%% Settings

tpm_type = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_threshSplit_binAverage_100perState';

%% Setup

% Load parameters
params = load([tpm_type '/params.mat']);

% Source directory
source_dir = [tpm_type '/'];

%% Setup

dims = size(params.thresh_values);
dims = circshift(dims, -1); % shift threshold dimension to the end
dims = dims(2:end); % Get rid of channel dimension

% Full dimensions
channel_sets = params.networks;
nChannels = size(channel_sets, 2);
nStates = 2^nChannels;
dims_state_ind = [size(channel_sets, 1) dims];
dims_state_dep = [nStates dims_state_ind];

% State-by-state TPMs
tpms = zeros([nStates nStates dims_state_ind]);
tpms(:) = nan; % NaN for missing TPMs
state_counters = zeros([nStates dims_state_ind]);
state_counters(:) = nan;

%% Join split results

for fly = 1 : size(params.thresh_values, 4)
    tic;
    for condition = 1 : size(params.thresh_values, 5)
        for trial = 1 : size(params.thresh_values, 3)
            for network_c = 1 : size(params.networks, 1)
                for thresh_c = 1 : length(params.thresh_ps)
                    
                    % File name
                    source_file = [...
                        'fly' num2str(fly)...
                        '_condition' num2str(condition)...
                        '_trial' num2str(trial)...
                        '_network' num2str(network_c)...
                        'tau' num2str(1)...
                        '_thresh' num2str(thresh_c)...
                        ];
                    
                    try
                        % Load file
                        tmp = load([source_dir source_file '.mat']);
                    catch ME
                        disp(ME.message);
                        disp(['Failed file: ' source_file]);
                        continue; % Skip to the next file, leave the data entry structure for this entry empty
                    end
                    
                    % Place into large data structure
                    tpms(:, :, network_c, trial, fly, condition, thresh_c) = tmp.tpm;
                    state_counters(:, network_c, trial, fly, condition, thresh_c) = tmp.state_counters;
                end
            end
        end
    end
    toc
end

%% Save

save([source_dir 'tpms.mat'], 'tpms', 'state_counters', '-v7.3');
