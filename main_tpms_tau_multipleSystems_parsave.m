function [] = main_tpms_tau_multipleSystems_parsave(filename, tpm, state_counters, nValues, fly, condition, trial, network, downMethod)
% Function for saving within parallel loop, check main script for required
% inputs

save(filename,...
    'tpm',...
    'state_counters',...
    'nValues',...
    'fly',...
    'condition',...
    'trial',...
    'network',...
    'downMethod');

end

