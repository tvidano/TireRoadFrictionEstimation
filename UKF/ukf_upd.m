function [XKK,PKK] = ukf_upd(model_param,X_HATK,PK,uprev,yk)
% PERFORMS MEASUREMENT UPDATE STEP IN UKF
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   X_HAT {vector} A Priori State Estimation [K|K-1],
%   P {matrix} A Priori State Covariance Matrix [K|K-1],
%   u {vector}: A Priori Input Vector [k-1],
%   yk {vector}: measurement 
%
% OUTPUTS: 
%   X_HAT {vector}: A Posteriori State Estimation [K|K],
%   P {matrix}: A Posteriori State Covariance Matrix [K|K]

% ---------------------------------------------------------------
%global R
% Unpack Model Parameters:
R = model_param.R;
N = model_param.N;

%% CHOLESKY DECOMP, NEW SET OF SIGMA POINTS
choles = chol(N*PK);

sigX_HAT = zeros(1,2*N);
for i = 1:1:N   % x[k|k-1]
    % result = 2*N sigma point vector 
    sigX_HAT(i)= X_HATK + choles(i,:).';  % i-th row of Cholesky decomposition
    sigX_HAT(N+i) = X_HATK - choles(i,:).';
    
end

%% Propagate through OUTPUT eqn

y_hatk = zeros(1,2*N);
for j = 1:1:2*N   % y[k|k-1]
    % progagate each sig pt
    y_hatk(j) = output_eqn(model_param, sigX_HAT(j), uprev); 
    
end

%% Compute sample mean, variances

ymean = sum(y_hatk)/(2*N);
xyvar = zeros(1,2*N);
yvar = xyvar;
for ii = 1:1:2*N
    % compute those variances
    xyvar(ii) = (sigX_HAT(ii) - X_HATK) * (y_hatk(ii) - ymean).'; 
    yvar(ii) = (y_hatk(ii) - ymean) * (y_hatk(ii) - ymean).';

end

Pxy = sum(xyvar)/(2*N);
Py = sum(yvar)/(2*N) + R;

%% Posteriori Estimate

XKK = X_HATK + (Pxy/Py)*(yk - ymean);
PKK = PK - (Pxy/Py)*Pxy.';

