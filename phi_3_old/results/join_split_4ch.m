%% Description

%{

Joins 4-channel results (across all parameters)

%}

%% Setup

nComponents = 2;

nChannels = 2;
global_tpm = 1;
flies = (1:10);
conditions = (1);
taus = (1:1000); %[1 4 100 200 300 400 500 600 700 800 900 1000]; % [1 2 3 4 5 10 20 30 40 50 75 100 125 150 175 200 225 250];
tau_type = 'bin'; % 'step' or 'bin'
bin_offsets = 1; % Currently this script only works for a single tau_offset
trials = (1);

% _nChannels4_globalTPM1_f01c2tauBin4500tauOffset21s1036t1

source_dir = '../results_split_nRuns10/';

prefix = 'unidir_no_inst_nSamples18000_nRuns10';
infix = '';
suffix = [...
    '_nChannels' num2str(nChannels)...
    '_globalTPM' num2str(global_tpm)...
    ];

source_prefix = [prefix infix suffix];

if strcmp(tau_type, 'step')
    tau_string = 'tau';
    binOffset_string = '';
    sbs_tpm = 1;
else % strcmp(tau_type, 'bin')
    tau_string = 'tauBin';
    binOffset_string = 'binOffset1';
    sbs_tpm = 0; % Currently phi-code stores SBN TPM when binning
end

channel_sets = nchoosek((1:nComponents), nChannels);
nStates = 2^nChannels;
dims_state_ind = [size(channel_sets, 1) length(trials), length(flies), length(conditions), length(taus)];
dims_state_dep = [nStates dims_state_ind];
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

output_file = [prefix infix '_phithree' suffix '_' binOffset_string '.mat'];

addpath('../../'); % For TPM conversion

%% Join split results

% NOTE: currently the actual parameters are used for indexing, so they need
%   from 1 and increment by 1
phis = cell(1);
phis{1} = struct();
phis{1}.nChannels = int8(nChannels);
phis{1}.channel_sets = int8(channel_sets);
phis{1}.taus = taus;

phis{1}.phis = single(zeros(dims_state_ind));
phis{1}.big_mips = cell(dims_state_dep);
phis{1}.big_mips = single(zeros([nStates 2 nConcepts dims_state_ind]));
phis{1}.state_counters = int16(zeros(dims_state_dep));
phis{1}.state_phis = single(zeros(dims_state_dep));
phis{1}.tpms = single(zeros([nStates dims_state_dep])); % Additional dimensions because TPM is state-by-state
for fly = flies
    disp(['fly' num2str(fly)]);
    for condition = conditions
        for tau_counter = 1 : length(taus)
            tau = taus(tau_counter);
            for set_counter = 1 : size(channel_sets, 1)
                disp(['set' num2str(set_counter)]);
                for trial = trials
                    
                    % File name
                    source_file = [source_prefix...
                        '_f' sprintf('%02d', fly)...
                        'c' num2str(condition)...
                        tau_string num2str(tau)...
                        binOffset_string...
                        's' sprintf('%04d', set_counter)...
                        't' num2str(trial)...
                        '.mat'...
                        ];
                    
                    try
                        % Load file
                        tmp = load([source_dir source_file]);
                    catch ME
                        disp(ME.message);
                        disp(['Failed file: ' source_file]);
                        continue; % Skip to the next file, leave the data entry structure for this entry empty
                    end
                    
                    % Place into large data structure
                    phis{1}.phis(set_counter, trial, fly, condition, tau_counter) = single(tmp.phi.phi);
                    phis{1}.state_counters(:, set_counter, trial, fly, condition, tau_counter) = int16(tmp.phi.state_counters);
                    phis{1}.big_mips(:, :, :, set_counter, trial, fly, condition, tau_counter) = constellation_parse(tmp.phi.big_mips, concept_list_full);
                    phis{1}.state_phis(:, set_counter, trial, fly, condition, tau_counter) = single(tmp.phi.state_phis);
                    if sbs_tpm == 1
                        phis{1}.tpms(:, :, set_counter, trial, fly, condition, tau_counter) = single(tmp.phi.tpm);
                    else % sbs_tpm == 0
                        phis{1}.tpms(:, :, set_counter, trial, fly, condition, tau_counter) = tpm_sbn2sbs(single(tmp.phi.tpm));
                    end
                    
                end
            end
        end
    end
end

%% Save

save(output_file, 'phis', '-v7.3');

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