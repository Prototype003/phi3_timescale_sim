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