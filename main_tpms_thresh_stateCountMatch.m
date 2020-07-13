%% Description

%{
Builds TPMs
%}

%% Settings

% Number of samples to use per source state (at time t)
samplesPerState = 200;

% Binarise method
binariseMethod = 'threshSplit';

% Binarise parameters
threshs = (10:2:90);
threshs = (10:2:90);
binariseParams = struct();

% Timescales
taus = (1);
tau = 1;

% Downsampling method
downMethod = 'binAverage';

%% Setup

prefix = 'NLbidirNoInstOrder1Thresh0-6SpikeReset_nSamples200000_nRuns10';

% Source data
source_dir = 'sim_data/';
source_file = [prefix '.mat'];
load([source_dir source_file]);

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
    thresh_min = prctile(data(:), min(threshs));
    thresh_max = prctile(data(:), max(threshs));
    
    % Round to nearest 0.05
    % Ceiling for min
    thresh_min = 5*ceil((thresh_min*100)/5)/100;
    % Floor for max
    thresh_max = 5*floor((thresh_max*100)/5)/100;
    
    % Create new thresholds
    threshs = (thresh_min:0.05:thresh_max);
    thresh_values = repmat(threshs', [1 size(data, 2), size(data, 4), length(taus)]);
else
    % Thresh values determined when building TPMs
    thresh_values = zeros(length(threshs), size(data, 2), size(data, 4), length(taus));
end

%% Build TPMs with constant number of samples per source state

broken = 0;
for thresh_c = 1 : length(threshs)
    thresh = threshs(thresh_c);
    disp(thresh);
    binariseParams.thresh = thresh;
    tic;
    
    if strcmp(downMethod, 'binAverage')
        step = 1;
    else
        step = tau;
    end
    
    for run = 1 : size(data, 4)
        
        % Downsample data
        if tau > 1
            if strcmp(downMethod, 'binAverage')
                cum = cumsum(data(:, :, 1, run), 1);
                data_resampled = cum(tau:tau:end, :) / tau;
                data_resampled = cat(1,...
                    data_resampled(1, :),...
                    data_resampled(2:end, :) - data_resampled(1:end-1, :));
            end
        else % No need to downsample
            data_resampled = data(:, :, 1, run);
        end
        
        % Binarise data
        if strcmp(binariseMethod, 'medianSplit')
            data_binarised = binarise_median(data_resampled);
            nValues = 2;
        elseif strcmp(binariseMethod, 'percSplit')
            [data_binarised, thresh_values(thresh_c, :, run ,tau)] = binarise_perc(data_resampled, binariseParams);
            nValues = 2;
        elseif strcmp(binariseMethod, 'threshSplit')
            data_binarised = binarise_threshValue(data_resampled, binariseParams);
            nValues = 2;
        end
        
        % Build TPM
        try
            [tpm, state_counters] = build_tpm_equalStates(data_binarised, step, nValues, samplesPerState);
        catch
            thresh_c = thresh_c - 1;
            threshs = threshs(1:thresh_c);
            broken = 1;
            break
        end
        
        % Save TPM
        out_suffix = [...
            'tau' num2str(tau)...
            '_thresh' num2str(thresh)...
            '_run' num2str(run)];
        save([out_dir out_suffix], 'tpm', 'state_counters', 'nValues');
        
    end
    
    if broken == 1
        break;
    end
    
    toc
end

%% Save general parameters

nRuns = size(data, 4);
save([out_dir 'params.mat'],...
    'samplesPerState',...
    'binariseMethod',...
    'binariseParams',...
    'threshs',...
    'taus',...
    'downMethod',...
    'nRuns',...
    'thresh_values');