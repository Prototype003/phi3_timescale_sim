%% Description

%{

Plot TPMs

%}

%% Setup

tpm_type = 'NLbidirNoInstOrder1Thresh1_nSamples200000_nRuns10_percSplit_binAverage';
tpm_type = 'NLbidirNoInstOrder1Thresh1_nSamples200000_nRuns10_percSplit_binAverage_500perState';

tau = 1;
thresh = 50;
run = 1;

addpath('../');

%% Load TPM

tpm_source = [...
    'tau' num2str(tau)...
    '_thresh' num2str(thresh)...
    '_run' num2str(run)];

load([tpm_type '/' tpm_source '.mat']);

%% Plot TPM

% State-by-state TPM
figure;
imagesc(tpm); colorbar;

% State-by-node TPM
figure;
imagesc(tpm_sbs2sbn(tpm)); colorbar;

% Origin state distribution
figure; bar((1:size(tpm, 1)), state_counters);
xlabel('state'); ylabel('occurrences');
