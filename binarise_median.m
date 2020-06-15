function [binarised] = binarise_median(data)
% Binarised data based on median (median across first dimension
%
% Inputs:
%   fly_data = matrix (time-samples x ...)
%
% Outputs:
%   binarised = matrix (time-samples x ...)

% Get median across time
medians = repmat(median(data, 1), [size(data, 1) ones(1, length(size(data))-1)]);

% Binarised based on median
binarised = data > medians;

end