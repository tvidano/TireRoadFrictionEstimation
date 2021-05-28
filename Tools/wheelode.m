function dydt = wheelode(t,y,model_param,inputs,mu)
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

% Unpack inputs
U = y(1);
w = y(2);
if t > 14
    Tm = -500;
else
    Tm = 0;
end
% [~,i] = min(abs(t - inputs.time));
% Tm = inputs.torque(i);
% Tm = interp1(inputs.time, inputs.torque, t);


slip = @(w,U) (r_e*w/U - 1);
% Fx = mu*Fz*calc_slip(w,U);
Fx = mu*Fz*sin(C*atan(B*slip(w,U) - E*(B*slip(w,U) - atan(B*slip(w,U)))));
dU = Fx/(m/4);
domega = (Tm - r_e*Fx)/J;

% Pack Outputs:
dydt(1) = dU;
dydt(2) = domega;
dydt = dydt';
end

