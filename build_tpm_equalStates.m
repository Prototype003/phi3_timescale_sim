function [ tpm, state_counter] = build_tpm_equalStates(fly_data, tau, n_values, n_samples)
% NOTE: Requires implementation of pyphi.convert.state2loli_index() from pyphi
%
% Builds a tpm for one fly and one condition
% Based off of build_tpm in fly_phi.py
%
% Builds TPMs such that the probability distributions of each state are
%   built using the same number of states
%
% Inputs:
%   fly_data = matix (of discretised data) with dimensions (samples x channels)
%       Holds data for one fly, one condition
%   tau = integer - the lag between current and future states
%   n_values = number of states each *node* can be in (e.g. 2 for ON and OFF)
%   n_samples = number of samples per source state (at time t)
%
% Outputs:
%   tpm = matrix with dimensions (n_values^channels x n_values^channels)
%       Each row holds the probabilities of a past state transitioning into future states (columns)

% Determine number of possible system states
n_states = n_values ^ size(fly_data, 2);

% Declare TPM
tpm = zeros(n_states, n_states);

state_counter = zeros(n_states, 1); % transition counter for each state

for state_current = 1 : n_states
    
    % Identify the actual state
    state_sample = loli_index2state(state_current, size(fly_data, 2));
    state_sample = logical(state_sample);
    
    % Filter data to get n_samples occurrences of state
    state_match = ismember(fly_data, state_sample, 'rows');
    state_positions = find(state_match);
    
    % Check that there are enough occurrences of the state
    if length(state_positions) < n_samples
        error(['number of state occurrences (' num2str(length(state_positions)) ') < n_samples']);
    end
    
    % Iterate through occurrences of state and track transitions
    for occ = 1 : n_samples
        
        % Identify position of occurrence
        position = state_positions(occ);
        
        % Identify future state
        sample_future = fly_data(position+tau, :);
        state_future = state2loli_index(sample_future);
        
        % Increment TPM transition by 1
        tpm(state_current, state_future) = tpm(state_current, state_future) + 1;
        
        % Increment transition counter
        state_counter(state_current) = state_counter(state_current) + 1;
        
    end
    
end

% Divide elements in TPM by transition counter
% If counter is 0, then transition never occurred - to avoid dividing by 0
for counter = 1 : length(state_counter)
    if state_counter(counter) == 0
        state_counter(counter) = 1;
    end
end
for future_state = 1 : size(tpm, 2)
    tpm(:, future_state) = tpm(:, future_state) ./ state_counter;
end

end

