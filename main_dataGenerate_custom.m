%% Description

%{

Generate non-linear auto-regressive model and simulate

General form - y = Ax + Bx' + e
    x' is nonlinearity

%}

%%

addpath('ar_nonlinear/');

%% Simulation parameters

% Time series length to simulate
nSamples = 200000;

% Number of time series to simulate
nRuns = 10;

%% Model parameters
model_type = ''; % name of the model

model_params = struct();

% Constants/starting points (nodes x 1)
model_params.constant = [0; 0];

% Linear connection matrix (nodes x nodes x time-lag)
% Only linear effects from self, none from others
model_params.A =[...
    -0.3 0;...
    0 -0.3];

% Nonlinear connection matrix (nodes x nodes x time-lag)
model_params.B = [...
    0 0.35;...
    0.35, 0];

% Instantaneous noise covariance (nodes x nodes)
% To generate noise from covariance:
%   https://au.mathworks.com/matlabcentral/answers/358204-how-to-generate-gaussian-noise-with-certain-covariance-and-zero-mean
%   See also MATLAB function simulate(varm)
model_params.Covariance = [...
    0.5 0;...
    0 0.5];

% Nonlinear function name ('NL' = 'nonlinear')
model_params.NL_func = 'thresh';
model_params.NL_params = struct();
model_params.NL_params.thresh = [1; 1]; % (nodes x 1)

%% Simulate data
% See MATLAB functions simulate(varm) and filter(varm) for reference

lags = size(model_params.A, 3);

% Same format as fly_data (samples x channels x trials x flies x conditions)
data = zeros(nSamples, size(model_params.A, 1), 1, nRuns, 1);
data_E = zeros(size(data)); % Contributions from noise

for run = 1 : nRuns
    tic;
    
    % Generate noise (see filter(); using ?lower Cholesky factor)
    noise_independent = randn(size(model_params.A, 1), nSamples); % (channels x samples), see filter()
    L = chol(model_params.Covariance, 'lower'); % Lower Cholesky factor
    noise = (L * noise_independent);
    
    % Store noise
    data_E(:, :, 1, run, 1) = noise'; % (samples x channels)
    
    % Prepend (zeros) and add starting conditions
    % Remember, noise is (channels x samples)
    Y = cat(2, zeros(size(model_params.A, 1), lags), noise + model_params.constant);
    
    % Step through time to determine samples
    for t = lags + 1 : nSamples + lags
        for lag = 1 : size(model_params.A, 3)
            % Add lagged components to the noise
            Y(:, t) = Y(:, t) +...
                model_params.A(:, :, lag) * Y(:, t-lag) +...
                model_params.B(:, :, lag) * feval(model_params.NL_func, Y(:, t-lag), model_params.NL_params);
        end
    end
    
    % Store simulated data
    data(:, :, 1, run, 1) = Y(:, lags+1:end)'; % (samples x channels); excluded prepended points

    toc
end

%% Save

