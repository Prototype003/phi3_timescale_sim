%% Description

%{

Summary figures

%}

%% Load

tpm_type = 'NLbidirNoInstOrder1Thresh0-9Lag10_nSamples200000_nRuns10_medianSplit_tauSearch_binAverage';
%tpm_type = 'NLbidirNoInstOrder1Thresh0-9Lag10_nSamples200000_nRuns10_medianSplit_tauSearch_tauStep';

tpm_type = 'forwardNLbidirNoInstOrder1Thresh0-9Lag9-11_nSamples200000_nRuns10_medianSplit_tauSearch_binAverage';
%tpm_type = 'forwardNLbidirNoInstOrder1Thresh0-9Lag9-11_nSamples200000_nRuns10_medianSplit_tauSearch_tauStep';

tpm_type = '3chforwardNLbidirNoInstOrder1Thresh0-9Lag9-11_nSamples200000_nRuns10_medianSplit_tauSearch_binAverage_nCh2';

load(['results/split/' tpm_type '/joined.mat']);
params = load(['../tpms/' tpm_type '/params.mat']);

%%

values = phis{1}.phis(:, :, :, :, :);
values_ref = values;

% treat (sub)networks as runs
values = permute(values, [6 2 3 1 4 5]);
dims = size(values);
values = reshape(values, [dims(1:2) dims(3)*dims(4) dims(5:6)]);

%% Average phi across runs

set_values = mean(values, 2); % average across trials
set_mean = mean(set_values, 3); % average across runs
set_std = std(set_values, [], 3);

cond_colours = {'r', 'b'};
figure;
for cond = 1 : size(values, 4)
    
    x = params.taus;
    xstring = 'tau';
    
    y = squeeze(values(:, :, :, cond, :));
    ymean = squeeze(mean(values(:, :, :, cond, :), 3));
    ybar = squeeze(std(values(:, :, :, cond, :), [], 3));
    
    plot(x, y, [cond_colours{cond} ':']);
    hold on;
    
    errorbar(x, ymean, ybar, cond_colours{cond});
    
    set(gca, 'XScale', 'log');
    
end

xlabel('\tau');
ylabel('\Phi');
title('tau-bin downsampling');

%% Small-phis

values = phis{1}.big_mips;

% Get state weighted average
state_counts = permute(phis{1}.state_counters, [1 7 8 2 3 4 5 6]);
state_counts = repmat(state_counts, [1 size(values, 2) size(values, 3) 1 1 1 1 1]);
values = double(values) .* double(state_counts);
values = sum(values, 1) ./ sum(state_counts, 1);
values = permute(values, [2 3 4 5 6 7 8 1]); % get rid of state dimension

% 1 = unpartitioned; 2 = partitioned
part = 1;
values = permute(values(part, :, :, :, :, :, :), [2 3 4 5 6 7 1]);

% move channel-set dimension to the front
values = permute(values, [2 1 3 4 5 6 7]);

% Treat channel-sets as runs
values = permute(values, [2 3 4 1 5 6]);
dims = size(values);
values = reshape(values, [dims(1:2) dims(3)*dims(4) dims(5:6)]);
values = permute(values, [6 1 2 3 4 5]);

%% Plot average phi per channel as function of threshold value

% 2ch = 3 mechanisms
mech_colours = {'r', 'b', 'k'};
mech_lines = {':', ':', '-'};

% 3ch = 7 mechanisms
%mech_colours = {'r', 'r', 'r', 'b', 'b', 'b', 'k'};
%mech_lines = {':', ':', ':', '-', '-', '-' '-'};

figure;

for cond = 1 : size(values, 5)
    
    x = params.taus;
    xstring = '\tau';
    
    for mech = 1 : size(values, 2)
        
        % Plot individual runs
        y = squeeze(mean(values(1, mech, :, :, cond, :), 3));
        plot(x, y, [mech_colours{mech} mech_lines{mech}], 'LineWidth', 0.1);
        hold on;
        
        % Plot average acros runs
        y = squeeze(mean(mean(values(1, mech, :, :, cond, :), 3), 4));
        plot(x, y, [mech_colours{mech} mech_lines{mech}], 'LineWidth', 2);
        
        % Plot errorbars across runs
        ybar = squeeze(std(mean(values(1, mech, :, :, cond, :), 3), [], 4));
        errorbar(x, y, ybar, mech_colours{mech}, 'LineWidth', 0.1);
        
    end
    
end

xlabel(xstring);
ylabel('\phi');

axis tight;

set(gca, 'XScale', 'log');


%% Average phi per channel

dims = size(values);
ch_values = zeros([max(params.networks(:)) dims(2:end)]);
for channel = 1 : max(params.networks(:))
    
    % Get networks which include channel
    ch_networks = any(params.networks == channel, 2);
    ch_networks = find(ch_networks);
    relevant = values(ch_networks, :, :, :, :, :);
    
    % Exclude any 0 values (because they are due to not computing
    %   due to few samples for the TPM
    relevant(relevant==0) = nan;
    
    % Average phi values across networks
    ch_values(channel, :, :, :, :, :) = nanmean(relevant, 1);
    
end

%% Plot average phi per channel as function of threshold value

% 2ch = 3 mechanisms
mech_colours = {'r', 'b', 'k'};
mech_lines = {':', ':', '-'};

% 3ch = 7 mechanisms
%mech_colours = {'r', 'r', 'r', 'b', 'b', 'b', 'k'};
%mech_lines = {':', ':', ':', '-', '-', '-' '-'};

figure;

for cond = 1 : size(ch_values, 5)
    
    x = params.taus;
    xstring = '\tau';
    
    for mech = 1 : size(ch_values, 2)
        
        % Plot individual runs
        y = squeeze(mean(ch_values(channel, mech, :, :, cond, :), 3));
        plot(x, y, [mech_colours{mech} mech_lines{mech}], 'LineWidth', 0.1);
        hold on;
        
        % Plot average acros runs
        y = squeeze(mean(mean(ch_values(channel, mech, :, :, cond, :), 3), 4));
        plot(x, y, [mech_colours{mech} mech_lines{mech}], 'LineWidth', 2);
        
        % Plot errorbars across runs
        ybar = squeeze(std(mean(ch_values(channel, mech, :, :, cond, :), 3), [], 4));
        errorbar(x, y, ybar, mech_colours{mech}, 'LineWidth', 0.1);
        
    end
    
end

xlabel(xstring);
ylabel('\phi');

axis tight;

set(gca, 'XScale', 'log');

% consistent y-axis across channels
%range = mean(ch_values, 2);
%ylim([min(range(:)) max(range(:))]);
