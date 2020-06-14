clear all;

addpath(genpath('../PhiToolbox'))

N = 6; % number of units
T = 10^6; % number of iterations
tau = 1; % time delay
N_st = 2;  % number of states

%% parameters
params.tau = tau;
params.number_of_states = N_st;

%% generate time series from Boltzman machine
beta = 4; % inverse temperature

% connectivity matrix
W = zeros(N,N);
community(1:ceil(N/2)) = 1;
community(ceil(N/2)+1:N) = 2;

for i=1: N
    for j=1: N
        if i~=j
            if community(i) == community(j)
                W(i,j) = 0.3;
            else
                W(i,j) = 0;
            end
        else
            W(i,i) = 0.1;
        end
    end
end
X = generate_Boltzmann(beta,W,N,T); 

%% 

T_seg = 1000;
figure(1)
t_vec1 = 1: T_seg;
t_vec2 = 2*10^3: 2*10^3+T_seg;
t_vec3 = 10^4: 10^4+T_seg;
t_vec4 = 10^5: 10^5+T_seg;
t_vec5 = T-300: T;
subplot(3,2,1),imagesc(X(:,t_vec1));
subplot(3,2,2),imagesc(X(:,t_vec2));
subplot(3,2,3),imagesc(X(:,t_vec3));
subplot(3,2,4),imagesc(X(:,t_vec4));
subplot(3,2,5),imagesc(X(:,t_vec5));
drawnow;

%% compute correlation 
R = corrcoef(X');
disp('Correlation Matrix')
disp(R);

%% find the Minimum Information Partition (MIP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.type_of_dist = 'discrete';
options.type_of_phi = 'star';

%%%%%%%%%% without pre-computed probability distributions %%%%%%%%%%
disp('Find the MIP without pre-computed probability distributions')
%% Exhaustive search %%
options.type_of_MIPsearch = 'Exhaustive';

disp('Exhaustive Search...')
tic;
[Z_MIP_without, phi_MIP_without] = MIP_search(X, params, options);
t_Exhaustive_without = toc;
disp(['Exhaustive Search finished. CalcTime=', num2str(t_Exhaustive_without)] )
disp(['phi at the MIP: ', num2str(phi_MIP_without)])
disp(['the MIP: ', num2str(Z_MIP_without)])
disp(' ')

%% Queyranne's algorithm %%
options.type_of_MIPsearch = 'Queyranne';

disp('Queyranne''s algorithm...')
tic;
[Z_MIP_Q_without, phi_MIP_Q_without] = MIP_search( X, params, options);
t_Queyranne_without = toc;
disp(['Queyranne''s algorithm Search finished. CalcTime=', num2str(t_Queyranne_without)])
disp(['phi at the MIP: ', num2str(phi_MIP_Q_without)])
disp(['the MIP: ', num2str(Z_MIP_Q_without)])
disp(' ')


% %%%%%%%%%% with pre-computed probability distributions %%%%%%%%%% 
% disp('Find the MIP with pre-computed probability distributions')
% 
% %% estimate probability distributions
% disp('Estimating probability distributions...')
% probs = data_to_probs(X, params, options);
% 
% %% Exhaustive Search %%
% disp('Exhaustive Search...')
% tic;
% [Z_MIP_with, phi_MIP_with, Zs, phis] = MIP_Exhaustive( probs, options );
% t_Exhaustive_with = toc;
% disp( ['Exhaustive Search finished. CalcTime=', num2str(t_Exhaustive_with)])
% disp(['phi at the MIP: ', num2str(phi_MIP_with)])
% disp(['the MIP: ', num2str(Z_MIP_with)])
% disp(' ')
% 
% %% Queyeranne's algorithm %%%
% disp('Queyranne''s algorithm...')
% tic;
% [Z_MIP_Q_with, phi_MIP_Q_with] = MIP_Queyranne( probs, options);
% t_Queyranne_with = toc;
% disp(['Queyranne''s algorithm finished. CalcTime=', num2str(t_Queyranne_with)])
% disp(['phi at the estimated MIP: ', num2str(phi_MIP_Q_with)])
% disp(['the estimated MIP: ', num2str(Z_MIP_Q_with)])
% disp(' ')
% 
