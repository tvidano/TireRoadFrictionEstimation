function [XKK,PKK] = ukf_upd(N,X_HATK,PK,uprev,yk)
% PERFORMS MEASUREMENT UPDATE STEP IN UKF

% INPUTS: N [dimension of inputs], X_HAT[K|K-1], P[K|K-1], u[k-1], yk measurement 

% OUTPUTS: X_HAT[K|K], P[K|K]

% ---------------------------------------------------------------
global R

%% CHOLESKY DECOMP, NEW SET OF SIGMA POINTS
choles = chol(N*PK);

for i = 1:1:N   % x[k|k-1]
    % result = 2*N sigma point vector 
    sigX_HAT(i)= X_HATK + choles(i,:).';  % i-th row of Cholesky decomposition
    sigX_HAT(N+i) = X_HATK - choles(i,:).';
    
end

%% Propagate through OUTPUT eqn

for j = 1:1:2*N   % y[k|k-1]
    % progagate each sig pt
    y_hatk(j) = output_eqn(sigX_HAT(j),uprev); 
    
end

%% Compute sample mean, variances

ymean = sum(y_hatk)/(2*N);
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

