%% Description

%{
Builds TPMs
%}

%% Settings

nChannels = 2;

% timescales
taus = 2.^(0:10);

% Binarise method
binariseMethod = 'medianSplit';

% Binarise parameters
thresh_ps = (50);
binariseParams = struct();

% Downsampling method
%   'tauStep' = step tau timesamples
%   'binAverage' = average across timesamples (tau samples per bin)
downMethod = 'binAverage';

%% Setup

prefix = 'split13500_bPlrRerefTyp1_lineNoiseRemoved_postPuffPreStim';

% Source data
source_dir = '../fly_phi/bin/workspace_results/';
source_dir = '../flies/PHDCOHEND-Q1326/postPuffPreStim/';
source_file = [prefix '.mat'];
load([source_dir source_file]);

data = fly_data;

% Output location
tpm_string = [...
    '_' binariseMethod...
    '_tauSearch'...
    '_' downMethod];
mkdir('tpms', [prefix tpm_string]);
out_dir = ['tpms/' prefix tpm_string '/'];

% Get channel sets
networks = nchoosek((1:size(data, 2)), nChannels);

%% Build TPMs with constant number of samples per source state

for fly = 1 : size(fly_data, 4)
    for condition = 1 : size(fly_data, 5)
        for trial = 1 : size(fly_data, 3)
            disp(['f' num2str(fly) 'c' num2str(condition) 't' num2str(trial)]);
            tic;
            for network_c = 1 : size(networks, 1)
                network = networks(network_c, :);
                net_data = data(:, network, trial, fly, condition);
                
                out_prefix = [...
                    'fly' num2str(fly)...
                    '_condition' num2str(condition)...
                    '_trial' num2str(trial)...
                    '_network' num2str(network_c)];
                
                for tau_c = 1 : length(taus)
                    tau = taus(tau_c);
                    
                    % Output filename
                    out_suffix = [...
                        'tau' num2str(tau_c)];
                    
                    if strcmp(downMethod, 'tauStep')
                        % Binarise data
                        if strcmp(binariseMethod, 'medianSplit')
                            data_binarised = binarise_median(net_data);
                            nValues = 2;
                        elseif strcmp(binariseMethod, 'threshSplit')
                            data_binarised = zeros(size(net_data));
                            ch_threshs = thresh_values(thresh_c, network, trial, fly, condition);
                            for ch = 1 : size(net_data, 2)
                                binariseParams.thresh = ch_threshs(1, ch);
                                data_binarised(:, ch) = binarise_threshValue(net_data(:, ch), binariseParams);
                            end
                            nValues = 2;
                        end
                        
                        % Build TPM
                        try
                            [tpm, state_counters] = build_tpm(data_binarised, tau, nValues);
                        catch ME
                            disp(ME.message);
                            disp(['Failed: ' out_prefix out_suffix]);
                            continue; % Skip to the next
                        end
                        
                    else % strcmp(downMethod, 'binAverage')
                        % Build TPM
                        if strcmp(binariseMethod, 'medianSplit')
                            binarise_string = 'median';
                        end
                        nValues = 2;
                        try
                            [tpm, state_counters] = build_tpm_binOffsets(net_data, tau, nValues, binarise_string);
                        catch ME
                            disp(ME.message);
                            disp(['Failed: ' out_prefix out_suffix]);
                            continue; % Skip to the next
                        end
                    end
                    
                    % Save TPM
                    save([out_dir out_prefix out_suffix],...
                        'tpm',...
                        'state_counters',...
                        'nValues',...
                        'fly',...
                        'condition',...
                        'trial',...
                        'network');
                    
                end
                
            end
            toc
        end
    end
end

%% Save general parameters

save([out_dir 'params.mat'],...
    'binariseMethod',...
    'taus',...
    'downMethod',...
    'networks');
