function [sbs] = tpm_sbn2sbs(sbn)
% Converts a state-by-node form TPM (using LOLI/fortran index ordering)
% into a state-by-state form TPM (also using LOLI indexing)
%
% Assumes binary states
%
% Inputs:
%   sbn = TPM; matrix with dimensions (states x nodes)
%       Each row holds the probabilities of each node transitioning into a
%           particular state at t+1, given a state at t
%
% Outputs:
%   sbs = TPM; matrix with dimensions (states x states)

n_states = size(sbn, 1);

% Determine number of nodes
n_nodes = log2(n_states); % this can be made dynamic (change the log base)?

% State-by-state TPM
sbs = zeros(n_states, n_states);

% Build list of possible states
% Iterate through each state and convert decimal index into binary string
state_list = zeros(size(sbn));
for state = 1 : n_states
    state_list(state, :) = dec2bi(state-1, n_nodes);
end

% Iterate through each state, and compute probability
% Note that due to state-by-node condition, new probabilities are all
% multiplications of combinations of (1-p) and p (p is probability a given
% node = 1 at next time step) - combination is given by the future state
% Also, combinations of (1-p) and p are identical for all rows
for state = 1 : n_states
    % 1 -> use p (1 -> 0 -> -p -> p)
    % 0 -> use 1-p (0 -> 1 -> 1-p -> 1-p)
    probs = abs(repmat(~state_list(state, :), [n_states, 1]) - sbn);
    sbs(:, state) = prod(probs, 2);
end

    function [bit_string] = dec2bi(dec, nBits)
        % Converts decimal number into binary string with LOLI/fortran
        % order; lowest order bit is on the left
        %
        % Inputs:
        %   dec = decimal number; within range which can be represented by
        %       nBits bits
        %   nBits = number of bits to consider
        %
        % Outputs:
        %   bin = binary string of nBits bits
        
        bit_string = zeros(1, nBits);
        for bit = 1 : nBits
            fits = floor(dec / (2^(nBits-1))) > 0;
            if fits == 1
                bit_string(bit) = 1;
                dec = dec - 2^(nBits-1);
            end
            nBits = nBits - 1;
        end
        
        bit_string = fliplr(bit_string);
    end
end

