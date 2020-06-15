%% Description

%{

Summary figures

%}

%% Load

load('results/split/unidir_no_inst_nSamples200000_nRuns10_medianSplit_binAverage/joined.mat');

%% Plot phi as function of tau

figure;

values = permute(phis{1}.phis, [3 5 1 2 4]);

% Plot individual runs
plot(phis{1}.taus, values, ':');
hold on;
% Plot average across runs
plot(phis{1}.taus, mean(values, 1), 'LineWidth', 2);

title('tau-bin using all offsets (median-split)');
ylabel('\Phi');
xlabel('\tau (ms)');