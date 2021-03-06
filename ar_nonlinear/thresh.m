function [y_adj] = thresh(y,params)
%THRESH
%   Adjusts values:
%       If abovee or equal to threshold - preserve value
%       If below threshold - change to 0
%
% Inputs:
%   past = vector (channels x 1); values of each channel
%   params = struct;
%       .thresh = vector (channels x 1); thresholds for each channel
%
% Outputs:
%   y_adj = vector (channels x 1); adjusted values of each channel

greater = y >= params.thresh;

y_adj = y .* greater;

end

