function [XK,PK] = ekf_pred(model_param,X_HAT,P,uprev,Epr,state_eqn,ts)
% PERFORMS MODEL PREDICTION STEP IN EKF
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters, 
%   X_HAT {vector}: A Posteriori State Estimation [K-1|K-1],
%   P {matrix}: A Posteriori State Covariance Matrix [K-1|K-1],
%   uprev {vector}: Previous Input Vector [K-1]
%   E_pr {matrix}: re-linearized matrix at [K-1]
%   state_eqn {func. handle}: Handle to the function defining state
%       equation. Inputs are model_param, previous state, current input.
%   ts: scalar value; sampling interval
%
% OUTPUTS: 
%   XK {vector} A Priori State Estimation [K|K-1], 
%   PK {matrix} A Priori State Covariance Matrix [K|K-1],

% ---------------------------------------------------------------

Q = model_param.Q;

% Compute A'[K-1]
Apr = A_pr(model_param,X_HAT,ts);

% Propagate state estimate through nonlinear state equation
XK = state_eqn();

% Propagate variance
PK = (Apr * P * Apr.') + (Epr * Q * Epr.');

end

