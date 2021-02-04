function [distance] = channel_set_distance(channel_set)
% Finds sum of distances between all channels in the channel set
%
% Inputs:
%   channel_set: vector of channels, each channel label is representative of its physical position
%
% Outputs:
%   distance: sum of differences between all channel pairs in the channel_set

% Sort channel set in ascending order (then distance is always dest-source)
channel_set = sort(channel_set);

distance = 0;
for source = 1 : length(channel_set)
    for dest = source+1 : length(channel_set)
        distance = distance + (channel_set(dest) - channel_set(source));
    end
end

end

