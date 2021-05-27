function [outputArg1,outputArg2] = wheel_state_eqn(model_param,xprev,uprev)
%WHEEL_STATE_EQN Calls state equation for wheel model
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   xprev {vector}: Previous A Posteriori State Vector [k-1|k-1],
%   u {vector}: Previous Input Vector [k-1]
%
% OUTPUTS: 
%   x_hat {vector}: Current A Priori State Vector [k|k-1]

%-------------------------------------------------------
% Unpack Model Parameters:
C = model_param.C;
B = model_param.B;
E = model_param.E;
r_e = model_param.r_e;
J = model_param.J;
m = model_param.m;
Fz = model_param.Fz;

% Define helper functions:
s = @(w,U) r_e*w/U - 1;
Fx = @(U,w,mu) mu*Fz*sin(C*atan(B*(1 - E)*s(w,u) + E*atan(B*s(w,u))));

% State equations:


end

