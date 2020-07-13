function [binarised] = binarise_threshValue(data, params)
% Binarised data based on percentile threshold (in first dimension)
% Above thresh -> 1
% Equal/Below thresh -> 0
%
% Inputs:
%   data = matrix (time-samples x ...)
%   params =  struct;
%       .thresh = threshold value
%
% Outputs:
%   binarised = matrix (time-samples x ...)

% Binarised based on threshold value
binarised = data > params.thresh;

end

