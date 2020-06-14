%% Description

%{

Computes GC for data at varying timescales

%}

%% Setup

load('unidir_no_inst_nSamples18000_phithree_nChannels2_globalTPM1_binOffset1.mat');

taus = phis{1}.taus;

%% Load timeseries

load('sim_data/unidir_no_inst_nSamples18000');

%% Downsample and compute GC for each tau

for tau_c = 1 : length(taus)
    tau = taus(tau_c);
    
    
end