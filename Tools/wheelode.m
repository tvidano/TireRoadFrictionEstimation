function [dydt,ext] = wheelode(t,y,model_param,inputs)
%WHEELODE ODE of single wheel with pacejka tire model for use in ODE45.
%
% INPUTS:
%   t: time,
%   y: packed states, see unpack states section,
%   model_param {struct}: A struct containing model parameters,
%   inputs {struct}: A struct containing the inputs to the simulation.
%
% OUTPUTS:
%   dydt: derivatives of packed states, see pack outputs section,
%   ext: extraneous variables useful for plotting and debugging.
%
% ---------------------------------------------------------------

%% Unpack model parameters
C = model_param.C;
B = model_param.B;
E = model_param.E;
r_e = model_param.r_e;
J = model_param.J;
m = model_param.m;
Fz = model_param.Fz;
mu = model_param.mu;

%% Unpack states:
U = y(1);
w = y(2);
if w <= 0
    w = 0;
end
if U <= 0
    U = 0;
end

%% Unpack inputs:
tm = inputs.time;
T = inputs.torque;
torque = interp1(tm,T,t);

%% Define Helper Functions:
calc_slip = @(w,U) r_e*w/U - 1;
get_force = @(U,w,mu) mu*Fz*sin(C*atan(B*(1 - E)*calc_slip(w,U)...
                      + E*atan(B*calc_slip(w,U))));

%% Cont. Time Equations:
Fx = -get_force(U,w,mu);
dU = Fx/(m/4);
domega = (torque - r_e*Fx)/J;
% Prevent negative ang. vel.
if (w == 0) && domega < 0
    domega = 0;
end

%% Pack Outputs:
dydt(1) = dU;
dydt(2) = domega;
dydt = dydt';

ext(1) = torque;
end

