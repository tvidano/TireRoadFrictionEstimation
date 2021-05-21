function [samp_mean,samp_var] = ukf_pred(N,X_HAT,P,uprev)
% PERFORMS MODEL PREDICTION STEP IN UKF

% INPUTS: N [dimension of inputs], X_HAT[K-1|K-1], P[K-1|K-1], u[k-1] 

% OUTPUTS: X_HAT[K|K-1], P[K|K-1]

% ---------------------------------------------------------------

global Q

%% Compute Cholesky Decomp, Form 2*N Sigma Points
if P == 0
    choles = 0;
else
    choles = chol(N*P);
end

for i = 1:1:N   % [k-1|k-1]
    % result = 2*N sigma point vector 
    sigX_HAT(i)= X_HAT + choles(i,:).';  % i-th row of Cholesky decomposition
    sigX_HAT(N+i) = X_HAT - choles(i,:).';
    
end

%% Propagate through STATE eqn

for j = 1:1:2*N   % [k|k-1]
    % progagate each sig pt
    x_hatk(j) = state_eqn(sigX_HAT(j),uprev); 
    
end

%% Compute sample mean, variance

samp_mean = sum(x_hatk)/(2*N);
for ii = 1:1:2*N  % [k|k-1]
    % to compute samp var
    var(ii) = (x_hatk(ii) - samp_mean) * (x_hatk(ii) - samp_mean).';
    
end
samp_var = sum(var)/(2*N) + Q;
