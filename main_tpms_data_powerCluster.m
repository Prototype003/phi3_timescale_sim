%% Description

%{
Builds TPMs

Binarises by computing power and clustering
%}

%% Settings

nChannels = 2;

% Binarise method
binariseMethod = 'powerCluster'; %

% Timescales
taus = (4);
tau = 4;

% Downsampling method
downMethod = 'tauStep'; % 'binAverage' = average across timesamples (tau samples per bin)

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
    '_tau' num2str(tau)...
    '_' downMethod];
mkdir('tpms', [prefix tpm_string]);
out_dir = ['tpms/' prefix tpm_string '/'];

% Get channel sets
networks = nchoosek((1:size(data, 2)), nChannels);

%% Concatenate trials

data = permute(fly_data, [1 3 2 4 5]);
dims = size(data);
data = reshape(data, [dims(1)*dims(2) dims(3:end)]);

% turn trailing dimensions into a single dimension
dims = size(data);
data = reshape(data, [dims(1) prod(dims(2:end))]);

%% Compute power spectrum at each time window
% Assumes 1kHz sampling rate

Fs = 1000;
winSize = 1000; % size (in samples) of window
winStep = 1; % step (in number of samples) of window

nWins = size(data, 1) - winSize + 1; % total number of windows

params = struct();
params.Fs = Fs;
params.tapers = [2 3];
params.fpass = [1 100];

% get size of spectrum
[S, f] = mtspectrumc(data(1:winSize, :), params);

Ss = NaN([size(S) nWins]);

parfor win_c = 1 : nWins
    disp(win_c);
    tic;
    
    window = (win_c : win_c + winSize - 1);
    [S, f] = mtspectrumc(data(window, :), params);
    
    Ss(:, :, win_c) = S;
    
    toc
end

%% PCA on power spectra

dims = size(Ss);
[coeff, score, latent, tsquared, explained] = pca(log(reshape(permute(Ss, [3 2 1]), [dims(3)*dims(2) dims(1)])));


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
                        else
                            data_resampled = net_data;
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
                        [tpm, state_counters] = build_tpm(data_binarised, step, nValues);
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
    'binariseMethod',...
    'binariseParams',...
    'threshs',...
    'thresh_ps',...
    'taus',...
    'downMethod',...
    'thresh_values',...
    'networks');
