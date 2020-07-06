%% Description

%{
Builds TPMs
%}

%% Settings

% Number of samples to use per source state (at time t)
samplesPerState = 500;

% Binarise method
binariseMethod = 'percSplit';

% Binarise parameters
threshs = (10:90);
binariseParams = struct();

% Timescales
taus = (1);
tau = 1;

% Downsampling method
downMethod = 'binAverage';

%% Setup

prefix = 'NLbidirNoInstOrder1Thresh1_nSamples200000_nRuns10';

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

%% Build TPMs with constant number of samples per source state

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
            data_binarised = binarise_perc(data_resampled, binariseParams);
            nValues = 2;
        end
        
        % Build TPM
        [tpm, state_counters] = build_tpm_equalStates(data_binarised, step, nValues, samplesPerState);
        
        % Save TPM
        out_suffix = [...
            'tau' num2str(tau)...
            '_thresh' num2str(thresh)...
            '_run' num2str(run)];
        save([out_dir out_suffix], 'tpm', 'state_counters', 'nValues');
        
    end
    
    toc
end

%% Save general parameters

nRuns = size(data, 4);
save([out_dir 'params.mat'], 'samplesPerState', 'binariseMethod', 'binariseParams', 'threshs', 'taus', 'downMethod', 'nRuns');