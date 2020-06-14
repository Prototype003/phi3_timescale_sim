function [ sbn ] = tpm_sbs2sbn(sbs)
% Converts a state-by-state form TPM (using LOLI/fortran index ordering)
% into a state-by-node form TPM (also using LOLI indexing)
%
% Assumes binary states
%
% Inputs:
%   sbs = TPM; matrix with dimensions (states x states)
%
% Outputs:
%   sbn = TPM; matrix with dimensions (states x nodes)
%       Each row holds the probabilities of each node transitioning into a
%           particular state at t+1, given a state at t

n_states = size(sbs, 1);

% Determine number of nodes
n_nodes = log2(n_states); % this can be made dynamic (change the log base)?

% State-by-node TPM
sbn = zeros(n_states, n_nodes);

% Build list of possible states
% Iterate through each state and convert decimal index into binary string
state_list = zeros(size(sbn));
for state = 1 : n_states
    state_list(state, :) = dec2bi(state-1, n_nodes);
end

% Iterate through each node, sum probabilities row-wise across columns
% where corresponding bit is '1'
for node = 1 : n_nodes
    on_columns = state_list(:, node) == 1;
    on_columns = repmat(on_columns', [n_states, 1]);
    tmp_sbs = sbs;
    tmp_sbs(~logical(on_columns)) = 0;
    sbn(:, node) = sum(tmp_sbs, 2);
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

