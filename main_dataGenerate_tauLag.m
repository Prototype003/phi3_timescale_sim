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
model_type = '3chforwardNLbidirNoInstOrder1Thresh0-9Lag9-11'; % name of the model

model_params = struct();

% Constants/starting points (nodes x 1)
model_params.constant = [0; 0];
model_params.constant = [0; 0; 0];

% Linear connection matrix (nodes x nodes x time-lag)
% Only linear effects from self, none from others
model_params.A =[...
    -0.1 0;...
    0 -0.1];

model_params.A =[...
    -0.1 0 0;...
    0 -0.1 0;...
    0 0 -0.1];

% Nonlinear connection matrix (nodes x nodes x time-lag)
model_params.B = [...
    0 0.9;...
    0.9, 0];

model_params.B = [...
    0 0.4 0.4;...
    0.4 0 0.4;...
    0.4 0.4 0];

% Instantaneous noise covariance (nodes x nodes)
% To generate noise from covariance:
%   https://au.mathworks.com/matlabcentral/answers/358204-how-to-generate-gaussian-noise-with-certain-covariance-and-zero-mean
%   See also MATLAB function simulate(varm)
model_params.Covariance = [...
    0.5 0;...
    0 0.5];

model_params.Covariance = [...
    0.5 0 0;...
    0 0.5 0;...
    0 0 0.5];

% Nonlinear function name ('NL' = 'nonlinear')
model_params.NL_func = 'thresh';
model_params.NL_params = struct();
model_params.NL_params.thresh = [0.9; 0.9]; % (nodes x 1)
model_params.NL_params.thresh = [0.9; 0.9; 0.9]; % (nodes x 1)

% If value exceeds thresohld, spike to given value
model_params.NL_params.spike = 0;
model_params.NL_params.spikeVals = model_params.NL_params.thresh; % (nodes x 1)

% Reset if thresholds exceeded
model_params.NL_params.reset = 0;
model_params.NL_params.resetVals = -model_params.NL_params.thresh;

% Set time delay for causal interactions
model_params.t_params.delay = [9 10 11]; % timesteps before 
model_params.t_params.delay_probs = [25 50 25];

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

lags = size(model_params.A, 3); % this is AR order

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
    prepend_length = max(lags, max(model_params.t_params.delay));
    Y = cat(2, zeros(size(model_params.A, 1), prepend_length), noise + model_params.constant);
    
    % Step through time to determine samples in the future
    %   At t, determine sample in the future at t+lag
    for t = 1 : nSamples
        
        % Check if spiking
        if model_params.NL_params.spike == 1
            for node = 1 : size(model_params.A, 1)
                if Y(node, t) >= model_params.NL_params.thresh(node)
                    Y(node, t) = model_params.NL_params.spikeVals(node);
                end
            end
        end
        
        % Determine delays per every pairing of nodes
        if length(model_params.t_params.delay) == 1
            delays = repmat(model_params.t_params.delay, [size(model_params.A, 1) size(model_params.A, 1)]);
        else
            delays = randsample(model_params.t_params.delay, size(model_params.A, 1).^2, true, model_params.t_params.delay_probs);
        end
        delays = reshape(delays, [size(model_params.A, 1) size(model_params.A, 1)]);
        
        % Iterate through each AR lag
        for lag = 1 : size(model_params.A, 3)
            
            % Node delays are independent
            for node = 1 : size(model_params.A, 1) % for each source nocde
                for link = 1 : size(delays, 2)
                    nl = feval(model_params.NL_func, Y(node, t), model_params.NL_params);
                    % Add current sample to noise in the future
                    Y(link, t+delays(link, node)) = Y(link, t+delays(link, node)) +...
                        model_params.A(link, node, lag) * Y(node, t) +...
                        model_params.B(link, node, lag) .* nl(node);
                end
            end
            
        end
        
    end
    
    % Check if spiking, for remaining t
    for t = nSamples+1 : size(Y, 2)
        if model_params.NL_params.spike == 1
            for node = 1 : size(model_params.A, 1)
                if Y(node, t) >= model_params.NL_params.thresh(node)
                    Y(node, t) = model_params.NL_params.spikeVals(node);
                end
            end
        end
    end
    
%     % Step through time to determine samples based on the past
%     %   At t, look backwards to determine sample at t
%     for t = lags + prepend_length + 1 : nSamples + lags
%         
%         % Determine delay
%         if length(model_params.t_params.delay) == 1
%             delay = model_params.t_params.delay;
%         else
%             delay = randsample(model_params.t_params.delay, 1, true, model_params.t_params.delay_probs);
%         end
%         
%         for lag = 1 : size(model_params.A, 3)
%             
%             % Add lagged components to the noise
%             Y(:, t) = Y(:, t) +...
%                 model_params.A(:, :, lag) * Y(:, t-(lag-1)-delay) +...
%                 model_params.B(:, :, lag) * feval(model_params.NL_func, Y(:, t-(lag-1)-delay), model_params.NL_params);
%             
%             % Spike
%             if model_params.NL_params.spike == 1
%                 for node = 1 : size(model_params.A, 1)
%                     if Y(node, t) >= model_params.NL_params.thresh(node)
%                         Y(node, t) = model_params.NL_params.spikeVals(node);
%                     end
%                 end
%             end
%             
%             % Reset if previous value met threshold
%             % Resets will override all of the previous things
%             %   TODO - move this part to the top?
%             if model_params.NL_params.reset == 1
%                 for node = 1 : size(model_params.A, 1)
%                     if Y(node, t-1) >= model_params.NL_params.thresh(node)
%                         Y(node, t) = model_params.NL_params.resetVals(node);
%                     end
%                 end
%             end
%             
%         end
%     end
    
    % Store simulated data
    data(:, :, 1, run, 1) = Y(:, prepend_length+1:end)'; % (samples x channels); excluded prepended points

    toc
end

%% Save

out_dir = 'sim_data/';
out_file = [model_type '_nSamples' num2str(nSamples) '_nRuns' num2str(nRuns)];

save([out_dir out_file], 'model_type', 'nSamples', 'model_params', 'data', 'data_E');

disp(['saved: ' out_file]);
