function [samp_mean,samp_var] = ukf_pred(model_param,X_HAT,P,uprev,state_eqn)
% PERFORMS MODEL PREDICTION STEP IN UKF
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters, 
%   X_HAT {COL vector}: A Posteriori State Estimation [K-1|K-1],
%   P {matrix} = A Posteriori State Covariance Matrix [K-1|K-1],
%   uprev {vector} = Previous Input Vector [K-1]
%   state_eqn {func. handle}: Handle to the function defining state
%       equation. Inputs are model_param, previous state, current input.
%
% OUTPUTS: 
%   X_HAT {vector} A Priori State Estimation [K|K-1], 
%   P {matrix} A Priori State Covariance Matrix [K|K-1],

% ---------------------------------------------------------------

% Unpack Model Parameters:
Q = model_param.Q;
N = model_param.N;

%% Compute Cholesky Decomp, Form 2*N Sigma Points
if P == zeros(N,N)
    choles = zeros(N,N);
else
    choles = chol(N*P);
end

sigX_HAT = zeros(N,2*N); % Each column contains a state vector sample
for i = 1:1:N   % [k-1|k-1]
    % result = 2*N sigma point vector 
    sigX_HAT(:,i)= X_HAT + choles(i,:).';  % i-th row of Cholesky decomposition
    sigX_HAT(:,N+i) = X_HAT - choles(i,:).';
    
end

%% Propagate through STATE eqn

x_hatk = zeros(N,2*N);
for j = 1:1:2*N   % [k|k-1]
    % progagate each sig pt
    x_hatk(:,j) = state_eqn(sigX_HAT(:,j),uprev,model_param); 
    
end

%% Compute sample mean, covariance

samp_mean = mean(x_hatk,2);  % compute mean along dir 2 (mean of each row)
VARS = zeros(N,N); % matrix of variances
for i = 1:1:2*N  % [k|k-1], i-th page of 3D variance array
    % to compute samp var
    vars(:,:,i) = (x_hatk(:,i) - samp_mean) * (x_hatk(:,i) - samp_mean).';
    VARS = vars(:,:,i) + VARS;
end
samp_var = (VARS ./ (2*N)) + Q;  % N x N matrix
