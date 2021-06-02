function [dydt,ext] = wheelode(t,y,model_param,inputs)
%WHEELODE Ode of single wheel with pacejka tire model
%   Detailed explanation goes here

% Unpack model parameters
C = model_param.C;
B = model_param.B;
E = model_param.E;
r_e = model_param.r_e;
J = model_param.J;
m = model_param.m;
Fz = model_param.Fz;
mu = model_param.mu;

% Unpack states:
U = y(1);
w = y(2);

% Unpack inputs:
% if t > 14
%     Tm = -4000;
% else
%     Tm = 0;
% end
tm = inputs.time;
T = inputs.torque;

torque = interp1(tm,T,t);

% Cont. Time Equations:
slip = @(w,U) (r_e*w/U - 1);
Fx = -mu*Fz*sin(C*atan(B*slip(w,U) - E*(B*slip(w,U) - atan(B*slip(w,U)))));
dU = Fx/(m/4);
domega = (torque - r_e*Fx)/J;

% Pack Outputs:
dydt(1) = dU;
dydt(2) = domega;
dydt = dydt';

ext(1) = torque;
end

