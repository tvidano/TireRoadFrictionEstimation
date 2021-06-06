function [XKK,PKK] = ukf_upd(model_param,X_HATK,PK,uprev,yk,output_eqn)
% PERFORMS MEASUREMENT UPDATE STEP IN UKF
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   X_HATK {COL vector} A Priori State Estimation [K|K-1],
%   PK {matrix} A Priori State Covariance Matrix [K|K-1],
%   uprev {vector}: A Priori Input Vector [k-1],
%   yk {COL vector}: measurement 
%   output_eqn {func. handle}: Handle to the function defining output
%       equation. Inputs are model_param, previous state, current input.
%
% OUTPUTS: 
%   X_HAT {COL vector}: A Posteriori State Estimation [K|K],
%   P {N x N matrix}: A Posteriori State Covariance Matrix [K|K]

% ---------------------------------------------------------------

% Unpack Model Parameters:
R = model_param.R;      % Output Covariance Matrix
N = model_param.N;      % Number of states
M = model_param.M;      % Number of output measurements

%% CHOLESKY DECOMP, NEW SET OF SIGMA POINTS
try
    choles = chol(N*PK);
catch
    choles = chol(1e-8*eye(1) + N*PK);
end

sigX_HAT = zeros(N,2*N);
for i = 1:1:N   % x[k|k-1]
    % result = 2*N sigma point vector 
    sigX_HAT(:,i)= X_HATK + choles(i,:).';  % i-th row of Cholesky decomposition
    sigX_HAT(:,N+i) = X_HATK - choles(i,:).';
    
end

%% Propagate through OUTPUT eqn

y_hatk = zeros(M,2*N);
for j = 1:1:2*N   % y[k|k-1]
    % progagate each sig pt
    y_hatk(:,j) = output_eqn(model_param, sigX_HAT(:,j), uprev); 
    
end

%% Compute sample mean, variances

ymean = mean(y_hatk,2); % compute mean along dir 2 (mean of each row)
xyvar = zeros(N,M); % variance is N x M matrix
yvar = zeros(M,M);
for ii = 1:1:2*N  % ii-th page of the 3D variance matrix
    % compute those variances
    XYvar(:,:,ii) = (sigX_HAT(:,ii) - X_HATK) * (y_hatk(:,ii) - ymean).'; 
    Yvar(:,:,ii) = (y_hatk(:,ii) - ymean) * (y_hatk(:,ii) - ymean).';
    xyvar = XYvar(:,:,ii) + xyvar; 
    yvar = Yvar(:,:,ii) + yvar;
    
end

Pxy = xyvar ./ (2*N);
Py = (yvar ./ (2*N)) + R;

%% Posteriori Estimate

XKK = X_HATK + (Pxy/Py)*(yk - ymean);
PKK = PK - (Pxy/Py)*Pxy.';

