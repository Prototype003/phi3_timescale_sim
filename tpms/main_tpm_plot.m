%% Description

%{

Plot TPMs

%}

%% Setup

tpm_type = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_threshSplit_binAverage_100perState';

addpath('../');

%% Load TPM

load([tpm_type '/' 'tpms.mat']);

params = load([tpm_type '/' 'params.mat']);

%% Plot TPM

network = 60;
fly = 1;
condition = 1;

figure;
vals = nanmean(tpms(:, :, network, :, fly, condition, :), 4); % Average across trials
clim = [min(vals(:)) max(vals(:))];
for thresh = 1 : size(tpms, 7)
    subplot(3, 3, thresh);
    imagesc(vals(:, :, 1, 1, 1, 1, thresh), clim); colorbar;
    %imagesc(tpm_sbs2sbn(vals(:, :, 1, 1, 1, 1, thresh)), clim); colorbar;
    title(['thresh at ' num2str(params.thresh_ps(thresh)) 'th %tile']);
end
