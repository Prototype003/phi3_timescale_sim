%% Description

%{
Builds TPMs
%}

%% Settings

nChannels = 2;

% Number of samples to use per source state (at time t)
samplesPerState = 100;

% Binarise method
binariseMethod = 'threshSplit'; % 'percSplit'

% Binarise parameters
thresh_ps = (30:5:70);
binariseParams = struct();

% Timescales
taus = (1);
tau = 1;

% Downsampling method
downMethod = 'binAverage';

%% Setup

prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

% Source data
source_dir = '../fly_phi/bin/workspace_results/';
source_file = [prefix '.mat'];
load([source_dir source_file]);

data = fly_data;

% Output location
tpm_string = [...
    '_' binariseMethod...
    '_' downMethod...
    '_' num2str(samplesPerState) 'perState'];
mkdir('tpms', [prefix tpm_string]);
out_dir = ['tpms/' prefix tpm_string '/'];

% Convert percentile thresholds to actual values
% Overwrites original threshs vectors
if strcmp(binariseMethod, 'threshSplit')
    
    data = permute(data, [1 3 2 4 5]); % samples x trials x channels x flies x conditions
    data_dims = size(data);
    data = reshape(data, [data_dims(1)*data_dims(2) data_dims(3) data_dims(4) data_dims(5)]); % Combine samples and trials
    
    dims = size(data);
    threshs = zeros([length(thresh_ps) dims(2:end)]);
    
    for thresh_c = 1 : length(thresh_ps)
        tic;
        thresh_p = thresh_ps(thresh_c);
        disp(thresh_p);
        thresh = prctile(data, thresh_p, 1);
        threshs(thresh_c, :, :, :) = thresh;
        toc
    end
    
    % Add trials dimension
    threshs = permute(threshs, [1 2 5 3 4]);
    thresh_values = repmat(threshs, [1 1 size(fly_data, 3) 1 1]);
    
    data = reshape(data, data_dims);
    data = permute(data, [1 3 2 4 5]); % samples x channels x trials x flies x conditions
    
else
    % Thresh values determined when building TPMs
end

% Get channel sets
networks = nchoosek((1:size(data, 2)), nChannels);

%% Build TPMs with constant number of samples per source state

for fly = 1 : size(fly_data, 4)
    for condition = 1 : size(fly_data, 5)
        tic;
        for trial = 1 : size(fly_data, 3)
            disp(['f' num2str(fly) 'c' num2str(condition) 't' num2str(trial)]);
            for network_c = 1 : size(networks, 1)
                network = networks(network_c, :);
                net_data = data(:, network, trial, fly, condition);
                
                out_prefix = [...
                    'fly' num2str(fly)...
                    '_condition' num2str(condition)...
                    '_trial' num2str(trial)...
                    '_network' num2str(network_c)];
                
                for thresh_c = 1 : size(thresh_values, 1)
                    
                    if strcmp(downMethod, 'binAverage')
                        step = 1;
                    else
                        step = tau;
                    end
                    
                    % Downsample data
                    if tau > 1
                        if strcmp(downMethod, 'binAverage')
                            cum = cumsum(net_data, 1);
                            data_resampled = cum(tau:tau:end, :) / tau;
                            data_resampled = cat(1,...
                                data_resampled(1, :),...
                                data_resampled(2:end, :) - data_resampled(1:end-1, :));
                        end
                    else % No need to downsample
                        data_resampled = net_data;
                    end
                    
                    % Binarise data
                    if strcmp(binariseMethod, 'medianSplit')
                        data_binarised = binarise_median(data_resampled);
                        nValues = 2;
                    elseif strcmp(binariseMethod, 'threshSplit')
                        data_binarised = zeros(size(net_data));
                        ch_threshs = thresh_values(thresh_c, network, trial, fly, condition);
                        for ch = 1 : size(net_data, 2)
                            binariseParams.thresh = ch_threshs(1, ch);
                            data_binarised(:, ch) = binarise_threshValue(data_resampled(:, ch), binariseParams);
                        end
                        nValues = 2;
                    end
                    
                    % Output filename
                    out_suffix = [...
                        'tau' num2str(tau)...
                        '_thresh' num2str(thresh_c)];
                    
                    % Build TPM
                    try
                        [tpm, state_counters] = build_tpm_equalStates(data_binarised, step, nValues, samplesPerState);
                    catch ME
                        disp(ME.message);
                        disp(['Failed: ' out_prefix out_suffix]);
                        continue; % Skip to the next
                    end
                    
                    % Save TPM
                    save([out_dir out_prefix out_suffix],...
                        'tpm',...
                        'state_counters',...
                        'nValues',...
                        'fly',...
                        'condition',...
                        'trial',...
                        'network',...
                        'ch_threshs');
                    
                end
                
            end
            
        end
        toc
    end
end

%% Save general parameters

save([out_dir 'params.mat'],...
    'samplesPerState',...
    'binariseMethod',...
    'binariseParams',...
    'threshs',...
    'thresh_ps',...
    'taus',...
    'downMethod',...
    'thresh_values',...
    'networks');
