%% Description

%{

Summary figures

%}

%% Load

tpm_type = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_threshSplit_binAverage_100perState';

load(['results/split/' tpm_type '/joined.mat']);
params = load(['../tpms/' tpm_type '/params.mat']);

%% Plot phi as a function of threshold value

flies = 1;
values = phis{1}.phis(:, :, flies, :, :);

% Average phi per channel
dims = size(values);
ch_values = zeros([max(params.networks(:)) dims(2:end)]);
for channel = 1 : max(params.networks(:))
    
    % Get networks which include channel
    ch_networks = any(params.networks == channel, 2);
    ch_networks = find(ch_networks);
    relevant = values(ch_networks, :, :, :, :);
    
    % Exclude any 0 values (because they are due to not computing
    %   due to few samples for the TPM
    relevant(relevant==0) = nan;
    
    % Average phi values across networks
    ch_values(channel, :, :, :, :) = nanmean(relevant, 1);
    
end

% Plot average phi as function of threshold value
cond_colours = {'r', 'b'};
figure;
for channel = 1 : max(params.networks(:))
    subplot(3, 5, channel);
    
    for condition = 1 : 1%size(ch_values, 3)
        
        x = params.thresh_values(:, channel, 1, flies, condition); % Assuming equal thresholds across trials
        
        % Plot individual trials
        y = squeeze(ch_values(channel, :, condition, :));
        plot(x, y, [cond_colours{condition} '-.']);
        hold on;
        
        % Plot average acros trials
        y = squeeze(mean(ch_values(channel, :, condition, :), 2));
        plot(x, y, [cond_colours{condition} '-'], 'LineWidth', 2);
        
    end
    
    title(['ch' num2str(channel)]);
    xlabel('thresh (V)');
    ylabel('\Phi');
    
    axis tight;
    
    % consistent y-axis across channels
    range = ch_values;
    ylim([min(range(:)) max(range(:))]);
    
end

%% Plot phi as function of actual threshold value

params = load(['../tpms/' tpm_type '/params.mat']);

% Assume consistent threshold values at each percentile, for both nodes
params.thresh_values = params.thresh_values(1:length(params.threshs), :, :);
threshs = mean(mean(params.thresh_values, 3), 2);

values = permute(phis{1}.phis, [3 5 1 2 4]);

figure;
% Plot individual runs
plot(threshs, values, ':');
hold on;
% Plot average across runs
plot(threshs, mean(values, 1), 'LineWidth', 2);

title(tpm_type);
ylabel('\Phi');
xlabel('thresh value');




