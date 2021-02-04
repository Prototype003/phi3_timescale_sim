%% Description

%{

Plot simulated data

%}

%% Setup

data_source = 'NLbidirNoInstOrder1Thresh0-9SpikeReset_nSamples200000_nRuns10';

load([data_source '.mat']);

%% Plot

figure;
channel_colours = {'r', 'b'};
for channel = 1 : size(data, 2)
    
    % Plot all runs
    %plot(samples, squeeze(data(:, channel, :, :)), [channel_colours{channel} ':']);
    %hold on;
    
    % Plot average across runs
    %plot(samples, squeeze(mean(data(:, channel, :, :), 4)), [channel_colours{channel} '-']);
    %hold on;
    
    % Plot one run
    run = 1;
    plot(squeeze(data(:, channel, :, run)), [channel_colours{channel} '-']); % simulated value
    hold on;
    plot(squeeze(data_E(:, channel, :, run)), [channel_colours{channel} '-.']); % noise value
    
end

% Limit x-axis
xlim([1 100]);

title(data_source);
xlabel('t');