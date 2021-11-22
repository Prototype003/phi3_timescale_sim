function [downsampled] = downsample_mean(data, tau, offset)
% Downsamples a series by averaging consecutive samples
%
% Inputs:
%   data = matrix (samples x channels)
%   tau = integer; downsampling factor
%   offset = integer; number of samples by which to offset before
%       downsampling
%
% Outputs:
%   downsampled = matrix (samples x channels)

% Figure out downsampled matrix size
d = data(offset+1:end, 1);
[dFrames, dropped] = buffer(d, tau); % Bin into frames of size tau

downsampled = nan(size(dFrames, 2), size(data, 2));

for ch = 1 : size(data, 2)
    d = data(offset+1:end, ch);
    
    % Bin into frames of size tau
    [dFrames, dropped] = buffer(d, tau);
    
    % Average samples per bin
    downsampled(:, ch) = mean(dFrames, 1);
end

end

