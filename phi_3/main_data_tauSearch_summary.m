%% Description

%{

Summary figures

%}

%% Load

tpm_type = 'split13500_bPlrRerefTyp1_lineNoiseRemoved_postPuffPreStim_medianSplit_tauSearch_tauStep';
tpm_type = 'split13500_bPlrRerefTyp1_lineNoiseRemoved_postPuffPreStim_medianSplit_tauSearch_binAverage';

load(['results/split/' tpm_type '/joined.mat']);
params = load(['../tpms/' tpm_type '/params.mat']);

%%

values = phis{1}.phis(:, :, :, :, :);
s = 24414.0625 / 4; % sampling rate

%% Average phi across channel sets

flies = (1);

set_values = mean(values(:, :, flies, :, :), 2); % average across trials
set_mean = mean(set_values, 1); % average across sets
set_std = std(set_values, [], 1);

cond_colours = {'r', 'b'};
figure;
for cond = 1 : size(values, 4)
    
    x = params.taus .* (1/s); % steps * sampling rate
    xstring = 'tau';
    
    y = squeeze(set_mean(:, :, :, cond, :));
    ybar = squeeze(set_std(:, :, :, cond, :));
    
    plot(x, y, cond_colours{cond});
    hold on;
    
    errorbar(x, y, ybar, cond_colours{cond});
    
    set(gca, 'XScale', 'log');
    
end

%% Phi across channel sets (wake/anest)

flies = (1);

set_values = mean(values(:, :, flies, 1, :) - values(:, :, flies, 2, :), 2); % average across trials
set_mean = mean(set_values, 1); % average across sets
set_std = std(set_values, [], 1);

cond_colours = {'r', 'b'};
figure;

x = params.taus .* (1/s); % steps * sampling rate
xstring = 'tau';

y = squeeze(set_values(:, :, :, :, :));
ybar = squeeze(set_std(:, :, :, :, :));

plot(x, y, cond_colours{cond});
hold on;

%errorbar(x, y, ybar, cond_colours{cond});

set(gca, 'XScale', 'log');

%% Average phi per channel

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

%% Plot average phi per channel as function of threshold value

flies = (1);
cond_colours = {'r', 'b'};
figure;
for channel = 1 : max(params.networks(:))
    subplot(3, 5, channel);
    
    for cond = 1 : 2%size(ch_values, 3)
        
        x = params.taus .* (1/s); % steps * sampling rate
        xstring = 'tau';
        
        % Plot individual trials
        y = squeeze(mean(ch_values(channel, :, flies, cond, :), 3));
        plot(x, y, [cond_colours{cond} '-.']);
        hold on;
        
        % Plot average acros trials
        y = squeeze(mean(mean(ch_values(channel, :, flies, cond, :), 3), 2));
        plot(x, y, [cond_colours{cond} '-'], 'LineWidth', 2);
        
        % Plot errorbars across trials
        ybar = squeeze(std(mean(ch_values(channel, :, flies, cond, :), 3), [], 2));
        errorbar(x, y, ybar, cond_colours{cond});
        
    end
    
    title(['ch' num2str(channel)]);
    xlabel(xstring);
    ylabel('\Phi');
    
    axis tight;
    
    set(gca, 'XScale', 'log');
    
    % consistent y-axis across channels
    %range = mean(ch_values, 2);
    %ylim([min(range(:)) max(range(:))]);
    
end

%% Plot average phi per channel as function of threshold value (wake / anest)
flies = (1);
cond_colours = {'r', 'b'};
figure;
for channel = 1 : max(params.networks(:))
    subplot(3, 5, channel);
    
    x = params.taus .* (1/s); % steps * sampling rate
    xstring = 'tau';
    
    % Plot individual trials
    y = squeeze(mean(ch_values(channel, :, flies, 1, :) ./ ch_values(channel, :, flies, 2, :), 3));
    plot(x, y, [cond_colours{cond} '-.']);
    hold on;
    
    % Plot average acros trials
    y = squeeze(mean(mean(ch_values(channel, :, flies, 1, :) ./ ch_values(channel, :, flies, 2, :), 2), 3));
    plot(x, y, [cond_colours{cond} '-'], 'LineWidth', 2);
    hold on;
    
    % Plot errorbars
    y = squeeze(mean(mean(ch_values(channel, :, flies, 1, :) ./ ch_values(channel, :, flies, 2, :), 2), 3));
    e = squeeze(std(mean(ch_values(channel, :, flies, 1, :) ./ ch_values(channel, :, flies, 2, :), 3), [], 2));
    errorbar(x, y, e, 'r');
    
    title(['ch' num2str(channel)]);
    xlabel(xstring);
    ylabel('\Phi');
    
    axis tight;
    
    set(gca, 'XScale', 'log');
    
    % consistent y-axis across channels
    %range = mean(ch_values, 2);
    %ylim([min(range(:)) max(range(:))]);
    
end

%% Plot as function of distance between channels
flies = (1);
cond = 2;
cond_colours = {'r', 'b'};

taus = repmat(params.taus .* (1/s), [size(phis{1}.channel_sets, 1), 1]);
dists = repmat(channel_set_distances(phis{1}.channel_sets), [1, length(params.taus)]);
vals = permute(mean(values(:, :, flies, cond, :), 2), [1 5 2 3 4]);

% % Plot at each tau, as function of distance
% figure;
% for tau_c = 1 : length(params.taus)
%     subplot(3, 4, tau_c);
%     x = channel_set_distances(phis{1}.channel_sets);
%     y = mean(values(:, :, flies, cond, tau_c), 2);
%     scatter(x, y, '.');
%     title(['\tau = ' num2str(params.taus(tau_c) * (1/s)) 's']);
% end

% Plot at each distance, as function of tau
figure;
for dist = 1 : max(dists)
    subplot(4, 4, dist);
    for cond = 1 : size(values, 4)
        x = params.taus .* (1/s);
        ids = dists(:, 1) == dist;
        y = permute(mean(values(ids, :, flies, cond, :), 2), [1 5 2 3 4]);
        plot(x, y, cond_colours{cond});
        hold on;
    end
    set(gca, 'XScale', 'log');
    title(['dist = ' num2str(dist)]);
end

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

flies = (1);
cond_colours = {'r', 'b'};
mech_lines = {':', ':', '-'};
figure;
for channel = 1 : max(params.networks(:))
    subplot(3, 5, channel);
    
    for cond = 1 : size(ch_values, 5)
        
        x = params.taus .* (1/s); % steps * sampling rate
        xstring = 'tau';
        
        for mech = 1 : size(ch_values, 2)
            
            % Plot individual trials
            y = squeeze(mean(ch_values(channel, mech, :, flies, cond, :), 4));
            plot(x, y, [cond_colours{cond} mech_lines{mech}], 'LineWidth', 0.1);
            hold on;
            
            % Plot average acros trials
            y = squeeze(mean(mean(ch_values(channel, mech, :, flies, cond, :), 4), 3));
            plot(x, y, [cond_colours{cond} mech_lines{mech}], 'LineWidth', 2);
            
            % Plot errorbars across trials
            ybar = squeeze(std(mean(ch_values(channel, mech, :, flies, cond, :), 4), [], 3));
            errorbar(x, y, ybar, cond_colours{cond}, 'LineWidth', 0.1);
            
        end
        
    end
    
    title(['ch' num2str(channel)]);
    xlabel(xstring);
    ylabel('\phi');
    
    axis tight;
    
    set(gca, 'XScale', 'log');
    %set(gca, 'YScale', 'log');
    
    % consistent y-axis across channels
    %range = mean(ch_values, 2);
    %ylim([min(range(:)) max(range(:))]);
    
end

%% Plot average phi per channel as function of threshold value (wake/anest)

flies = (1);
mech_colours = {'b', 'm', 'k'};
mech_lines = {':', ':', '-'};
figure;
for channel = 1 : max(params.networks(:))
    subplot(3, 5, channel);
    
    for cond = 1 : 1 %size(ch_values, 5)
        
        x = params.taus .* (1/s); % steps * sampling rate
        xstring = 'tau';
        
        for mech = 1 : size(ch_values, 2)
            
            % Plot individual trials
            y = squeeze(mean(ch_values(channel, mech, :, flies, 1, :) ./ ch_values(channel, mech, :, flies, 2, :), 4));
            plot(x, y, [mech_colours{mech} mech_lines{mech}], 'LineWidth', 0.1);
            hold on;
            
            % Plot average acros trials
            y = squeeze(mean(mean(ch_values(channel, mech, :, flies, 1, :) ./ ch_values(channel, mech, :, flies, 2, :), 4), 3));
            plot(x, y, [mech_colours{mech} mech_lines{mech}], 'LineWidth', 2);
            
            % Plot errorbars across trials
            ybar = squeeze(std(mean(ch_values(channel, mech, :, flies, 1, :) ./ ch_values(channel, mech, :, flies, 2, :), 4), [], 3));
            errorbar(x, y, ybar, mech_colours{mech}, 'LineWidth', 0.1);
            
        end
        
    end
    
    title(['ch' num2str(channel)]);
    xlabel(xstring);
    ylabel('\phi');
    
    axis tight;
    
    set(gca, 'XScale', 'log');
    %set(gca, 'YScale', 'log');
    
    % consistent y-axis across channels
    %range = mean(ch_values, 2);
    %ylim([min(range(:)) max(range(:))]);
    
end

%% Plot as function of distance between channels
flies = (1);
cond = 2;
cond_colours = {'r', 'b'};
mech_lines = {'-.', ':', '-'};
mech_lineWidths = [0.5 0.5 2];

taus = repmat(params.taus .* (1/s), [size(phis{1}.channel_sets, 1), 1]);
dists = repmat(channel_set_distances(phis{1}.channel_sets), [1, length(params.taus)]);

% Plot at each distance, as function of tau
figure;
for dist = 1 : max(dists)
    subplot(4, 4, dist);
    for cond = 1 : size(values, 5)
        for mech = 1 : size(values, 2)
            x = params.taus .* (1/s);
            ids = dists(:, 1) == dist;
            y = permute(mean(values(ids, mech, :, flies, cond, :), 3), [1 6 2 3 4 5]);
            %plot(x, y, [cond_colours{cond} mech_lines{mech}], 'LineWidth', 0.1);
            hold on;
            errorbar(x, mean(y, 1), [], std(y, [], 1), [cond_colours{cond} mech_lines{mech}], 'LineWidth', mech_lineWidths(mech));
        end
    end
    set(gca, 'XScale', 'log');
    %set(gca, 'YScale', 'log');
    title(['dist = ' num2str(dist)]);
    ylabel('\phi');
    xlabel('\tau');
end

%% Plot as function of distance between channels (wake/anest)
flies = (1);
cond = 2;
mech_colours = {'b', 'm', 'k'};
mech_lines = {':', ':', '-'};
mech_lineWidths = [1 1 1];

taus = repmat(params.taus .* (1/s), [size(phis{1}.channel_sets, 1), 1]);
dists = repmat(channel_set_distances(phis{1}.channel_sets), [1, length(params.taus)]);

% Plot at each distance, as function of tau
figure;
for dist = 1 : max(dists)
    subplot(4, 4, dist);
    for cond = 1 : 1%size(values, 5)
        for mech = 1 : size(values, 2)
            x = params.taus .* (1/s);
            ids = dists(:, 1) == dist;
            y = permute(mean(values(ids, mech, :, flies, 1, :) ./ values(ids, mech, :, flies, 2, :), 3), [1 6 2 3 4 5]);
            %plot(x, y, [cond_colours{cond} mech_lines{mech}], 'LineWidth', 0.1);
            hold on;
            errorbar(x, mean(y, 1), [], std(y, [], 1), [mech_colours{mech} mech_lines{mech}], 'LineWidth', mech_lineWidths(mech));
        end
    end
    set(gca, 'XScale', 'log');
    %set(gca, 'YScale', 'log');
    title(['dist = ' num2str(dist)]);
    ylabel('\phi');
    xlabel('\tau');
end
