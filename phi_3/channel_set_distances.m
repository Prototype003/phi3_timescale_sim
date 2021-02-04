function [distances] = channel_set_distances(channel_sets)
% Finds all distances for all channel sets
%
% Inputs:
%   channel_sets: matrix (sets x channels); each row gives a channel set
%       Channel labels should be reflective of physical position
%
% Outputs:
%   distances: vector (sets x 1); each entry is the global distance (sum of distances)
%       between each channel pair in the corresponding channel set

distances = zeros(size(channel_sets, 1), 1);
for channel_set = 1 : size(channel_sets, 1)
    distances(channel_set) = channel_set_distance(channel_sets(channel_set, :));
end

end

