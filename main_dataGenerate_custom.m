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
model_type = 'NLbidirNoInstOrder1Thresh0-6'; % name of the model
%model_type = 'NLbidirNoInstOrder1Thresh0-7SpikeReset';

model_params = struct();

% Constants/starting points (nodes x 1)
model_params.constant = [0; 0];

% Linear connection matrix (nodes x nodes x time-lag)
% Only linear effects from self, none from others
model_params.A =[...
    -0.1 0;...
    0 -0.1];

% Nonlinear connection matrix (nodes x nodes x time-lag)
model_params.B = [...
    0 0.9;...
    0.9, 0];

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
model_params.NL_params.thresh = [0.6; 0.6]; % (nodes x 1)

% If value exceeds thresohld, spike to given value
model_params.NL_params.spike = 0;
model_params.NL_params.spikeVals = model_params.NL_params.thresh; % (nodes x 1)

% Reset if thresholds exceeded
model_params.NL_params.reset = 0;
model_params.NL_params.resetVals = -model_params.NL_params.thresh;

% Check for stability (so system doesn't explode)
%   Not sure if check is correct (and what about imaginary eigenvalues?)
%   Check applies to linear systems
% Note for stability/stationarity:
%   https://kevinkotze.github.io/ts-7-var/#fn6 - see section 1.3
%   Lutkepohl, H. 2005. New Introduction to Multiple Time Series Analysis. Heidelberg: Springer-Verlang.
eigs = eig(det(model_params.A + model_params.B));
if any(eigs > 1)
    warning('VAR probably unstable');
end

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
            
            % Spike
            if model_params.NL_params.spike == 1
                for node = 1 : size(model_params.A, 1)
                    if Y(node, t) >= model_params.NL_params.thresh(node)
                        Y(node, t) = model_params.NL_params.spikeVals(node);
                    end
                end
            end
            
            % Reset if previous value met threshold
            % Resets will override all of the previous things
            %   TODO - move this part to the top?
            if model_params.NL_params.reset == 1
                for node = 1 : size(model_params.A, 1)
                    if Y(node, t-1) >= model_params.NL_params.thresh(node)
                        Y(node, t) = model_params.NL_params.resetVals(node);
                    end
                end
            end
            
        end
    end
    
    % Store simulated data
    data(:, :, 1, run, 1) = Y(:, lags+1:end)'; % (samples x channels); excluded prepended points

    toc
end

%% Save

out_dir = 'sim_data/';
out_file = [model_type '_nSamples' num2str(nSamples) '_nRuns' num2str(nRuns)];

save([out_dir out_file], 'model_type', 'nSamples', 'model_params', 'data', 'data_E');
