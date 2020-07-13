%% Description

%{

Summary figures

%}

%% Load

tpm_type = 'unidir_no_inst_nSamples200000_nRuns10_medianSplit_binAverage';
tpm_type = 'unidirNoInstOrder1_nSamples200000_nRuns10_medianSplit_binAverage';
%tpm_type = 'disconnectedNoInstOrder1_nSamples200000_nRuns10_medianSplit_binAverage';
%tpm_type = 'bidirNoInstOrder1_nSamples200000_nRuns10_medianSplit_binAverage';
tpm_type = 'NLbidirNoInstOrder1Thresh-inf_nSamples200000_nRuns10_percSplit_binAverage';
tpm_type = 'NLbidirNoInstOrder1Thresh1_nSamples200000_nRuns10_threshSplit_binAverage_200perState';
%tpm_type = 'NLbidirNoInstOrder1Thresh0-7SpikeReset_nSamples200000_nRuns10_threshSplit_binAverage_200perState';

load(['results/split/' tpm_type '/joined.mat']);

%% Plot phi as function of tau

% figure;
% 
% values = permute(phis{1}.phis, [3 5 1 2 4]);
% 
% % Plot individual runs
% plot(phis{1}.taus, values, ':');
% hold on;
% % Plot average across runs
% plot(phis{1}.taus, mean(values, 1), 'LineWidth', 2);
% 
% title(tpm_type);
% ylabel('\Phi');
% xlabel('\tau (ms)');
% 
% ylim([0 0.15]); xlim([0 20]);

%% Plot phi as function of percentile thresh

% figure;
% 
% values = permute(phis{1}.phis, [3 5 1 2 4]);
% 
% % Plot individual runs
% plot(phis{1}.threshs, values, ':');
% hold on;
% % Plot average across runs
% plot(phis{1}.threshs, mean(values, 1), 'LineWidth', 2);
% 
% title(tpm_type);
% ylabel('\Phi');
% xlabel('thresh percentile');
% 
% %ylim([0 0.15]); xlim([0 20]);

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




