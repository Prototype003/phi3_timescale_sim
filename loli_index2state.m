function [state] = loli_index2state(index, bits)
% Function to reverse state2loli_index()
% Converts decimal index to LOLI big-index (fortran order)
% Inputs:
%   index = corresponding loli index (1-indexed, not 0-indexed as in Python)
%   bits = number of bits to use
% Outputs:
%   state = vector of 1s and 0s; (1 x bits)
%
% TODO: flip bit-order if using C order (as opposed to Fortran
% order)

% LOLI indexing is fortran order (as opposed to C order)
index = index - 1; % convert to 0-indexing
state = bitget(index, (1:bits));

end