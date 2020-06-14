% compute several measures of integrated information in a multivariate autoregressive model, 
% X^t = A*X^{t-1} + E,
% where A is a connectivity matrix and E is gaussian noise.

clear all;

addpath(genpath('../PhiToolbox'))

%% parameters for computing phi
Z = [1 2]; % partition with which phi is computed
tau = 1; % time delay

%% important
params.tau = tau;

%% generate random gaussian time series X
N = 2; % the number of elements
A = [0.2 0.1; 0.5 0.2]; % connectivity matrix
Cov_E = eye(N,N); % covariance matrix of E
T = 10^5;
X = zeros(N,T);
X(:,1) = randn(N,1);
for t=2: T
    E = randn(N,1);
    X(:,t) = A*X(:,t-1) + E;
end


%% compute phi from time series 
options.type_of_dist = 'Gauss';
% type_of_dist:
%    'Gauss': Gaussian distribution
%    'discrete': discrete probability distribution

% available options of type_of_phi for Gaussian distributions:
%    'MI1': Multi (Mutual) information, e.g., I(X_1; X_2). (IIT1.0)
%    'SI': phi_H, stochastic interaction
%    'star': phi_star, based on mismatched decoding
%    'Geo': phi_G, information geometry version

% mutual information
options.type_of_phi = 'MI1';
MI = phi_comp(X, Z, params, options);

% stochastic interaction
options.type_of_phi = 'SI';
SI = phi_comp(X, Z, params, options);

% phi*
options.type_of_phi = 'star';
phi_star = phi_comp(X, Z, params, options);

% geometric phi
options.type_of_phi = 'Geo';
phi_G = phi_comp(X, Z, params, options);

%%
fprintf('MI=%f SI=%f phi*=%f phi_G=%f\n',MI,SI,phi_star,phi_G);

% %% compute phi from pre-computed covariance matrices
% % mutual information
% isjoint = 0;
% probs = Cov_comp(X,tau,isjoint);
% MI = MI1_Gauss(probs.Cov_X,Z);
% 
% % stochastic interaction
% isjoint = 1;
% probs = Cov_comp(X,tau,isjoint);
% SI = SI_Gauss(probs.Cov_X,probs.Cov_XY,probs.Cov_Y,Z);
% 
% % phi*
% phi_star = phi_star_Gauss(probs.Cov_X,probs.Cov_XY,probs.Cov_Y,Z);
% 
% % geometric phi
% type_of_phi = 'Geo';
% phi_G = phi_G_Gauss(probs.Cov_X,probs.Cov_XY,probs.Cov_Y,Z);
% 
% %%
% fprintf('MI=%f SI=%f phi*=%f phi_G=%f\n',MI,SI,phi_star,phi_G);
