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
% dU = model_param.dU;

Tire.D_x = mu*Fz; Tire.D_y = 0;
Tire.B_x = B; Tire.B_y = 0;
Tire.C_x = C; Tire.C_y = 0;
Tire.E_x = E; Tire.E_y = 0;
Tire.B_xalpha = 0; Tire.C_xalpha = 0; Tire.E_xalpha = 0;
Tire.B_ykappa = 0; Tire.C_ykappa = 0; Tire.S_Hykappa = 0;

%% Unpack states:
U = y(1);
w = y(2);
% if w <= 0
%     w = 0;
% end
% 
% if U <= 0
%     U = 0;
% end

%% Unpack inputs:
tm = inputs.time;
T = inputs.torque;
if length(tm) > 2
    torque = interp1(tm,T,t);
else
    torque = (T(2) - T(1))/(tm(2) - tm(1))*(t - tm(1)) + T(1);
end
% torque = interp1(tm,T,t);
% torque = -2000;

%% Define Helper Functions:
% calc_slip = @(w,U) -(U - r_e*w)/U;
% get_force = @(U,w,mu) ;

%% Cont. Time Equations:
% Prevent Divide by Zero:
if abs(U) < 1e-10
    U = 1e-10;
%     if dU < 0
%         U = -U;
%     end
end
if abs(w) < 1e-10
    w = 1e-10;
end
if U < model_param.r_e*w
    dU = 1;
else
    dU = -1;
end
if dU < 0
    kappa = -(U - r_e*w)/U;
elseif dU >= 0
    kappa = (r_e*w - U)/(r_e*w);
else
    kappa = 0;
end

[Fx,~,~,~] = PacSimple(Tire,Fz,0,kappa);
Fx = -Fx;
% Fx = 1*mu*Fz*kappa;
if Fx > 5e4
    Fx = 5e4;
end
dU = (Fx)/(m/4);
domega = (torque - r_e*Fx)/J;

if U < 1
    dU = 0;
end
if w < 0
    domega = 0;
end

%% Pack Outputs:
dydt(1) = dU;
dydt(2) = domega;
dydt = dydt';

ext(1) = torque;
end

