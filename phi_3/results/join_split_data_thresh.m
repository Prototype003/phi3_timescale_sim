%% Description

%{

Joins split results

%}

%% Settings

tpm_type = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_threshSplit_tau4_tauStep';
tau = 4;

%% Setup

% Load parameters
params = load(['../../tpms/' tpm_type '/params.mat']);

% Source directory
source_dir = ['split/' tpm_type '/'];

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

% List of possible concepts
nConcepts = 0; % power-set of nChannels
concept_list_full = cell(0); concept_list_counter = 1;
for concept_order = 1 : nChannels
    subsets = nchoosek((1:nChannels), concept_order);
    nConcepts = nConcepts + size(subsets, 1);
    for concept = 1 : size(subsets, 1)
        concept_list_full{concept_list_counter} = subsets(concept, :) - 1;
        concept_list_counter = concept_list_counter + 1;
    end
end

%% Join split results

% NOTE: currently the actual parameters are used for indexing, so they need
%   from 1 and increment by 1
phis = cell(1);
phis{1} = struct();
phis{1}.nChannels = int8(nChannels);
phis{1}.channel_sets = int8(channel_sets);
phis{1}.params = params;

phis{1}.phis = single(zeros(dims_state_ind));
phis{1}.big_mips = cell(dims_state_dep);
phis{1}.big_mips = single(zeros([nStates 2 nConcepts dims_state_ind]));
phis{1}.state_counters = int16(zeros(dims_state_dep));
phis{1}.state_phis = single(zeros(dims_state_dep));

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
                        'tau' num2str(tau)...
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
                    phis{1}.phis(network_c, trial, fly, condition, thresh_c) = single(tmp.phi.phi);
                    phis{1}.state_counters(:, network_c, trial, fly, condition, thresh_c) = int16(tmp.phi.state_counters);
                    phis{1}.big_mips(:, :, :, network_c, trial, fly, condition, thresh_c) = constellation_parse(tmp.phi.big_mips, concept_list_full);
                    phis{1}.state_phis(:, network_c, trial, fly, condition, thresh_c) = single(tmp.phi.state_phis);
                end
            end
        end
    end
    toc
end

%% Save

save([source_dir 'joined.mat'], 'phis', '-v7.3');

%% Parse constellation

function stripped = constellation_parse(big_mips, concept_list)
% Goes through big_mip and gets only the vital stuff
% Assumes that concepts are sorted, with the same order of 'concept_list'!
% This is important as details of the mechanism are discarded - only phi
% values for each concept are kept and returned
%
% Inputs:
%   big_mips: cell vector holding big_mip structs (output from pyphi.compute.big_mip)
%   concept_list: cell vector, each cell holds a concept-mechanism
%       concepts in 'big_mips' should be sorted to have the same order as this
%
% Outputs:
%   stripped: matrix (nMips x 2 x concepts)
%       nMips: length of 'big_mips' (for 4ch, 16)
%       2: unpartitioned (1), and partitioned (2) constellations
%       concepts: number of possible concepts (for 4ch, 15)

stripped = single(zeros(length(big_mips), 2, length(concept_list)));


for mip_counter = 1 : length(big_mips)
    
    unpart_counter = 1;
    part_counter = 1;
    for concept =  1 : length(concept_list)
        
        % unpartitioned_constellation
        if unpart_counter <= length(big_mips{mip_counter}.unpartitioned_constellation) % Ensure that there are still concepts to index into
            unpart_concept = big_mips{mip_counter}.unpartitioned_constellation{unpart_counter}.mechanism;
            if isequal(unpart_concept, concept_list{concept})
                stripped(mip_counter, 1, concept) = single(big_mips{mip_counter}.unpartitioned_constellation{unpart_counter}.phi);
                unpart_counter = unpart_counter + 1;
            else
                stripped(mip_counter, 1, concept) = 0;
            end
        else
            stripped(mip_counter, 1, concept) = 0;
        end
        
        % partitioned_constellation
        if part_counter <= length(big_mips{mip_counter}.partitioned_constellation) % Ensure that there are still concepts to index into
            part_concept = big_mips{mip_counter}.partitioned_constellation{part_counter}.mechanism;
            if isequal(part_concept, concept_list{concept})
                stripped(mip_counter, 2, concept) = single(big_mips{mip_counter}.partitioned_constellation{part_counter}.phi);
                part_counter = part_counter + 1;
            else
                stripped(mip_counter, 2, concept) = 0;
            end
        else
            stripped(mip_counter, 2, concept) = 0;
        end
        
    end
    
end

end