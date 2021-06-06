function [x_hat] = wheel_state_eqn(model_param,xprev,uprev)
%WHEEL_STATE_EQN Calls state equation for wheel model estimation.
%
% INPUTS: 
%   model_param {struct}: A struct containing model parameters,
%   xprev {vector}: Previous A Posteriori State Vector [k-1|k-1],
%   uprev {vector}: Previous Input Vector [k-1]
%
% OUTPUTS: 
%   x_hat {vector}: Current A Priori State Vector [k|k-1]

%-------------------------------------------------------
% Unpack Model Parameters:
% C = model_param.C;
% B = model_param.B;
% E = model_param.E;
% r_e = model_param.r_e;
% J = model_param.J;
% m = model_param.m;
% Fz = model_param.Fz;
% ts = model_param.ts;
tk = model_param.tk; % Current time
tj = model_param.tj; % Previous time

% Pack torque table for ODE inputs:
k = find(model_param.t == tk);
j = k - 1;
torques = [model_param.torque(j),model_param.torque(k)];
inputs = struct('time',[tj, tk],'torque',torques);

% Unpack previous states and inputs:
U = xprev(1);
w = xprev(2);
mu = xprev(3);
model_param.mu = mu; % Pass to ODE

% RK 4 Integration to Estimate Next Discrete State:
h = 2e-3;
tspan = tj:h:tk;
y0 = [U;w];
t = tspan;
y(1,:) = y0;
for i=1:(length(t)-1)
    K1 = wheelode(t(i),y(i,:),model_param,inputs)';
    K2 = wheelode(t(i)+0.5*h,y(i,:)+0.5*h*K1,model_param,inputs)';
    K3 = wheelode(t(i)+0.5*h,y(i,:)+0.5*h*K2,model_param,inputs)';
    K4 = wheelode(t(i)+h,y(i,:)+h*K3,model_param,inputs)';
    
    y(i+1,:) = y(i,:) + (1/6)*(K1+2*K2+2*K3+K4)*h;
end

% Unpack outputs:
U = y(:,1); 
w = y(:,2);

% Pack current states:
x_hat = [U(end),w(end),mu(end)]';

end

