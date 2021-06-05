function [XK,PK] = ekf_pred(model_param,X_HAT,P,uprev,Epr,state_eqn)
% PERFORMS MODEL PREDICTION STEP IN EKF
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters, 
%   X_HAT {COL vector}: A Posteriori State Estimation [K-1|K-1],
%   P {N x N matrix}: A Posteriori State Covariance Matrix [K-1|K-1],
%   uprev {vector}: Previous Input Vector [K-1]
%   E_pr {matrix}: re-linearized matrix at [K-1]
%   state_eqn {func. handle}: Handle to the function defining state
%       equation. Inputs are model_param, previous state, current input.
%
% OUTPUTS: 
%   XK {vector} A Priori State Estimation [K|K-1], 
%   PK {matrix} A Priori State Covariance Matrix [K|K-1],

% ---------------------------------------------------------------

Q = model_param.Q;
ts = model_param.tk - model_param.tj;  % to accomodate variable time step 

% Compute A'[K-1]
Apr = A_pr(model_param,X_HAT,ts);

% Propagate state estimate through nonlinear state equation
XK = state_eqn(model_param,X_HAT,uprev); 

% Propagate variance
PK = (Apr * P * Apr.') + (Epr * Q * Epr.');

end

