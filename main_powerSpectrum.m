%% Description

%{

Plots power spectrum of independent components

%}

%% Setup

prefix = 'NLbidirNoInstOrder1_nSamples200000_nRuns10';

% Source data
source_dir = 'sim_data/';
source_file = [prefix '.mat'];
load([source_dir source_file]);

% Sampling rate (Hz)
sample_rate = 1000;

%% Compute power spectra

L = size(data, 1); % Length of signal

Y = fft(data(:, :, :, :), [], 1);

% Two-sided spectrum
P2 = abs(Y/L);

% One-sided spectrum
P1 = P2(1:L/2+1, :, :, :);
P1(2:end-1, :, :, :) = 2*P1(2:end-1, :, :, :);

f = sample_rate*(0:(L/2))/L;

%% Plot power spectra

channel = 1;

cond_colours = {'r', 'b'};

figure;
for condition = 1 : size(P1, 5)
    plot_data = P1(:, channel, :, :);
    plot_data = mean(plot_data, 4); % average across runs
    plot_data = reshape(plot_data, [size(plot_data, 1) numel(plot_data)/size(plot_data, 1)]);
    plot(f, plot_data, cond_colours{condition}); hold on;
end
title('Single-Sided Amplitude Spectrum');
xlabel('f (Hz)');
ylabel('|P1(f)|');

%% Plot original signals

fly = 2;
trials = (8);

cond_colours = {'r:', 'b:'};

figure;
for condition = 1 : size(fly_data, 5)
    plot((1:size(fly_data, 1)), mean(fly_data(:, :, trials, fly, condition), 3), cond_colours{condition}); hold on;
end
xlim([1 size(fly_data, 1)]);
title('IC time-course');
xlabel('t');

%% Reformat fly data

fly_data = fly_data_original;

fly_data = permute(fly_data, [1 3 5 2 4]); % samples x trials x conditions x channels x flies

dims = size(fly_data);

fly_data = reshape(fly_data, [dims(1)*dims(2) dims(3) dims(4) dims(5)]);

fly_data = reshape(fly_data, [dims(1)*dims(2)*dims(3) dims(4) dims(5)]);

%% Autocorrelation

set(0, 'DefaultFigureVisible', 'on'); % acf function plots each time

nFlies = size(fly_data, 3); % Assuming samples x channels x flies
nChannels = size(fly_data, 2);
lags = 500;

timescales = zeros(lags+1, nChannels, nFlies);

figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    for channel = 1 : nChannels
        tic
        %timescales(:, fly) = acf(fly_data(:, channel, fly), lags);
        timescales(:, channel, fly) = autocorr(fly_data(:, channel, fly), 'NumLags', lags);
        autocorr(fly_data(:, channel, fly), 'NumLags', lags);
        hold on;
        toc;
    end
end

open_plots = findobj('type', 'figure'); % Use to get handles to all figures

% Plot
figure;
for fly = 1 : nFlies
    subplot(4, 4, fly)
    plot((1:lags+1), timescales(:, :, fly));
end
