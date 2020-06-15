function [ tpm, transition_counter ] = build_tpm_binOffsets(fly_data, tau, n_values, binarise_method)
% NOTE: Requires implementation of pyphi.convert.state2loli_index() from pyphi
%
% Builds a tpm for one fly and one condition
% Based off of build_tpm in fly_phi.py
%
% Inputs:
%   fly_data = matix (of discretised data) with dimensions (samples x channels)
%       Holds data for one fly, one condition
%   tau = integer - the lag between current and future states
%   n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
%   binarise_method = string; 'median' for binarising based on median or
%       'diff' for binarising based on gradient
%
% Outputs:
%   tpm = matrix with dimensions (n_values^channels x n_values^channels)
%       Each row holds the probabilities of a past state transitioning into future states (columns)

% Determine number of possible system states
n_states = n_values ^ size(fly_data, 2);

% Declare TPM
tpm = zeros(n_states, n_states);

transition_counter = zeros(n_states, 1); % transition counter for each state

for offset = 1 : tau
    fly_data_offset = fly_data(offset:end, :, :);
    
    % Downsample by averaging
    cum = cumsum(fly_data_offset, 1);
    fly_data_resampled = cum(tau:tau:end, :) / tau;
    fly_data_resampled = fly_data_resampled(2:end, :) - fly_data_resampled(1:end-1, :);
    fly_data_offset = fly_data_resampled;
    
    % Binarise
    if strcmp(binarise_method, 'median')
        fly_data_offset = binarise_median(fly_data_offset);
    else % strcmp(binarise_method, 'diff')
        error('binarise method not implemented here yet');
    end
    
    for sample = 1 : size(fly_data_offset, 1) - 1 % the last one doesn't transition
        sample_current = fly_data_offset(sample, :);
        sample_future = fly_data_offset(sample+1, :);
        
        % Identify current state
        state_current = state2loli_index(sample_current);
        
        % Identify future state
        state_future = state2loli_index(sample_future);
        
        % Increment TPM transition by 1
        tpm(state_current, state_future) = tpm(state_current, state_future) + 1;
        
        % Increment transition counter
        transition_counter(state_current) = transition_counter(state_current) + 1;
    end
    
end

% Divide elements in TPM by transition counter
% If counter is 0, then transition never occurred - to avoid dividing by 0
for counter = 1 : length(transition_counter)
    if transition_counter(counter) == 0
        transition_counter(counter) = 1;
    end
end
for future_state = 1 : size(tpm, 2)
    tpm(:, future_state) = tpm(:, future_state) ./ transition_counter;
end


    function [index] = state2loli_index(state)
        % Function to mimic pyphi.convert.state2loli_index()
        % Converts LOLI bit-index (fortran order) into decimal index
        % Inputs:
        %   state = vector of 1s and 0s
        % Outputs:
        %   index = corresponding loli index (1-indexed, not 0-indexed as in Python)
        
        % TODO: flip bit-order if using C order (as opposed to Fortran
        % order)
        
        % LOLI indexing is fortran order (as opposed to C order)
        index = 0;
        for bit = 1 : length(state)
            index = index + state(bit) * 2^(bit-1);
        end
        index = index + 1;
    end

end

