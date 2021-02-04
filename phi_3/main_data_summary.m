%% Description

%{

Summary figures

%}

%% Load

tpm_type = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_threshSplit_tau4_tauStep';

load(['results/split/' tpm_type '/joined.mat']);
params = load(['../tpms/' tpm_type '/params.mat']);

%% Plot phi as a function of threshold value

values = phis{1}.phis(:, :, :, :, :);

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
flies = (1:13);
cond_colours = {'r', 'b'};
figure;
for channel = 1 : max(params.networks(:))
    subplot(3, 5, channel);
    
    for condition = 1 : 2%size(ch_values, 3)
        
        x = mean(params.thresh_values(:, channel, 1, flies, condition), 4); % Assuming equal thresholds across trials
        xstring = 'thresh (V)';
        x = params.thresh_ps;
        xstring = '%tile thresh';
        
        % Plot individual trials
        y = squeeze(mean(ch_values(channel, :, flies, condition, :), 3));
        plot(x, y, [cond_colours{condition} '-.']);
        hold on;
        
        % Plot average acros trials
        y = squeeze(mean(mean(ch_values(channel, :, flies, condition, :), 2), 3));
        plot(x, y, [cond_colours{condition} '-'], 'LineWidth', 2);
        
    end
    
    title(['ch' num2str(channel)]);
    xlabel(xstring);
    ylabel('\Phi');
    
    axis tight;
    
    % consistent y-axis across channels
    %range = mean(ch_values, 2);
    %ylim([min(range(:)) max(range(:))]);
    
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

%% Plot phi as function of threshold, per fly

values = phis{1}.phis(:, :, :, :, :);

% Average across trials
values = mean(values, 2);
%values = values(:, :, :, 1, :) ./ values(:, :, :, 2, :); % wake/anest

cond_colours = {'r', 'b'};
figure;
for fly = 1 : size(values, 3)
    subplot(4, 4, fly);
    
    for condition = 1 : size(values, 4)
        
        x = params.thresh_ps;
        xstring = '%tile thresh';
        
        % Plot average across channel sets
        y = squeeze(mean(values(:, 1, fly, condition, :), 1));
        plot(x, y, [cond_colours{condition} '-'], 'LineWidth', 2);
        
        hold on
        
        % Plot errorpatch
        ybar = squeeze(std(values(:, 1, fly, condition, :), [], 1));
        patch(...
            [x fliplr(x)],...
            [y'+ybar' fliplr(y'-ybar')],...
            cond_colours{condition},...
            'FaceAlpha', 0.3,...
            'LineStyle', 'none');
        
    end
    
    title(['fly' num2str(fly)]);
    xlabel(xstring);
    ylabel('\Phi');
    
    axis tight;
    
end

%% Print

figure_name = 'figS_binThresh_raw';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG