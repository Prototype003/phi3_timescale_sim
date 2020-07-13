function [binarised, thresh_values] = binarise_perc(data, params)
% Binarised data based on percentile threshold (in first dimension)
% Above thresh -> 1
% Equal/Below thresh -> 0
%
% Inputs:
%   data = matrix (time-samples x ...)
%   params =  struct;
%       .thresh = percentile value (0-100)
%
% Outputs:
%   binarised = matrix (time-samples x ...)

% Get threshold values
thresh_values = prctile(data, params.thresh, 1);
threshs = repmat(thresh_values, [size(data, 1) ones(1, length(size(data))-1)]);

% Binarised based on median
binarised = data >= threshs;

end

