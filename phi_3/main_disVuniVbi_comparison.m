%% Description

%{

Compare disconnected vs unidirectionally connected vs bidirectionally
connected model phi

%}

%% Load

model_names = {...
    'NLnodirNoInstOrder1Thresh0-9_nSamples200000_nRuns10_percSplit_binAverage_200perState',...
    'NLunidirNoInstOrder1Thresh0-9_nSamples200000_nRuns10_percSplit_binAverage_200perState',...
    'NLbidirNoInstOrder1Thresh0-9_nSamples200000_nRuns10_percSplit_binAverage_200perState'};

models = cell(1);
for m = 1 : length(model_names)
    models{m} = load(['results/split/' model_names{m} '/joined.mat']);
end

%% Get specific threshold to look at

thresh = 50;

phis = zeros(size(models{1}.phis{1}.phis, 3), 3); % runs x models
for m = 1 : length(models)
    thresh_pos = find(models{m}.phis{1}.threshs == thresh);
    phis(:, m) = squeeze(models{m}.phis{1}.phis(1, 1, :, 1, thresh_pos));
end

%% Plot average big-phi per model

vals_mean = mean(phis, 1);
vals_std = std(phis, [], 1);

figure;
bar(vals_mean);
hold on;
errorbar((1:length(models)), vals_mean, vals_std, 'k.');

set(gca, 'XTick', (1:length(models)), 'XTickLabel', {'dis.', 'uni.' 'bi.'});
xlabel('model');
ylabel('\Phi');

%% Print

figure_name = 'figures/figS_disVuniVbi_raw';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG