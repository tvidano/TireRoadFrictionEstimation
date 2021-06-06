function [XKK, PKK] = ekf_upd(model_param, X_HATK,PK,uk,yk,C_pr,F_pr,output_eqn)
% PERFORMS MEASUREMENT UPDATE STEP IN EKF
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   X_HATK {COL vector} A Priori State Estimation [K|K-1],
%   PK {matrix} A Priori State Covariance Matrix [K|K-1],
%   uk {vector}: current time step Input Vector [k],
%   yk {COL vector}: measurement 
%   C_pr {M x N matrix}: constant
%   F_pr {M x M matrix}: constant
%   output_eqn {func. handle}: Handle to the function defining output
%       equation. Inputs are model_param, previous state, current input.
%
% OUTPUTS: 
%   XKK {COL vector}: A Posteriori State Estimation [K|K],
%   PKK {N x N matrix}: A Posteriori State Covariance Matrix [K|K]

% ---------------------------------------------------------------

% Unpack Model Parameters:
R = model_param.R;

% compute Kalman Gain
LK = PK*C_pr.'/(C_pr*PK*C_pr.' + F_pr*R*F_pr.');

% measurement update --> XKK
XKK = X_HATK + LK*(yk - output_eqn(X_HATK, uk, model_param));

% measure update --> PKK
PKK = PK - LK*C_pr*PK;

end

