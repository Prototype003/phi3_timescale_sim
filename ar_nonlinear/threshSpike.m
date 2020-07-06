function [y_adj] = threshSpike(y,params)
%THRESH
%   Adjusts values:
%       If above threshold - preserve value
%       If below threshold - change to 0
%
% Inputs:
%   past = vector (channels x 1); values of each channel
%   params = struct;
%       .thresh = vector (channels x 1); thresholds for each channel
%       .spike = vector (channels x 1); spiking value if threshold met
%
% Outputs:
%   y_adj = vector (channels x 1); adjusted values of each channel

greater = y > params.thresh;

y_adj = y;

y_adj(greater == 1) = params.spike;

end

