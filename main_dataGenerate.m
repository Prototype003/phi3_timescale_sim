%% Description

%{

Generates data from vector auto-regressive models

Models specified in Cohen et al. 2019
https://doi.org/10.1016/j.jneumeth.2019.108443

%}

%% Settings

% Specify which model to take from the Example_systems (from the paper)
model_type = 'unidir_no_inst';
model_type = 'unidirNoInstOrder1';
model_type = 'disconnected_with_inst';
model_type = 'disconnectedNoInstOrder1';
model_type = 'bidirNoInstOrder1';

% Time series length to simulate
nSamples = 200000;

% Number of time series to simulate
nRuns = 10;

%% Get model parameters
% From paper

% Each row gives params for a model in the paper
% A = matrix of regression coefficients (nVars x nVars x auto-regressive
%   order)
% SIG = covariance matrix of residuals (nVars x nVars)
model_params = Example_systems();

% Take desired model
desired = cellfun(@(x) strcmp(x, model_type), {model_params.nm});
model_params = model_params(desired);

%% Make model parameters
% If not using model from paper

if isempty(model_params)
    paper_models = Example_systems();
    
    % Coefficient matrix
    model_params(1).A =...
        paper_models(5).A(:, :, 1);
    
    % Innovation matrix
    model_params(1).SIG =...
        paper_models(4).SIG;
    
    % nm
    model_params(1).nm =...
        model_type;
    
end

%% Build model
% Assumes all orders up to the maximum order are used
% Uses 0 as intercepts

% Get regression coefficients into cell format
coeffs = cell(size(model_params.A, 3), 1);
for order = 1 : length(coeffs)
    coeffs{order} = model_params.A(:, :, order);
end

model = varm('AR', coeffs, 'Covariance', model_params.SIG, 'Constant', zeros(size(coeffs{1}, 1), 1));

%% Generate data

% Same format as fly_data (samples x channels x trials x flies x conditions)
data = zeros(nSamples, 2, 1, nRuns, 1);
data_E = zeros(size(data));

for run = 1 : nRuns
    tic;
    
    [Y, E] = simulate(model, nSamples);
    data(:, :, 1, run, 1) = Y;
    data_E(:, :, 1, run, 1) = E;
    
    toc
end

%% Save data

out_dir = 'sim_data/';
out_file = [model_type '_nSamples' num2str(nSamples) '_nRuns' num2str(nRuns)];

save([out_dir out_file], 'model_type', 'nSamples', 'model_params', 'model', 'data', 'data_E');