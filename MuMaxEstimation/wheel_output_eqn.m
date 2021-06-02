function [y_hat] = wheel_output_eqn(model_param,xprev,uprev)
%WHEEL_OUTPUT_EQN Calls output equation for wheel model estimation.
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   xprev {vector}: Previous A Posteriori State Vector [k-1|k-1],
%   uprev {vector}: Previous Input Vector [k-1]
%
% OUTPUTS: 
%   y_hat {vector}: Current A Priori Output Vector [k|k-1]

%-------------------------------------------------------
% Unpack Model Parameters:
C = model_param.C;
B = model_param.B;
E = model_param.E;
r_e = model_param.r_e;
J = model_param.J;
m = model_param.m;
Fz = model_param.Fz;
ts = model_param.ts;

% Unpack previous states and inputs:
U = xprev(1);
w = xprev(2);
mu = xprev(3);

torque = uprev;

% Define helper functions:
calc_slip = @(w,U) r_e*w/U - 1;
get_force = @(U,w,mu) mu*Fz*sin(C*atan(B*(1 - E)*calc_slip(w,U)...
                      + E*atan(B*calc_slip(w,U))));

% Euler Integration to Estimate Next Discrete State:
n = 10;
dt = ts/n;
for i = 2:n
    Fx = -get_force(U(i-1),w(i-1),mu);
    dU = Fx/(m/4);
    U(i) = U(i-1) + dU*(dt);
    
    domega = (torque - r_e*Fx)/J;
    w(i) = w(i-1) + domega*(dt);
    if w(i) <= 0 
        w(i) = 0;
    end
end

% Pack current states:
y_hat = [U(end),w(end)]';

end

