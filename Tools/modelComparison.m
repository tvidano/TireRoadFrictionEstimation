% Comparison of Measured Data from Multi-Body Dynamic Simulation, single
% wheel model solved with ODE45, and single wheel model solved with euler
% integration.
clear;clc;close all;
%% Collect measurement data:
mu8 = matfile('mu0.80.mat');
tm = mu8.t;
Um = mu8.U;
sm = mu8.s;
Fxm = mu8.Fx;
Tbm = mu8.Tb;
wm = mu8.w;
Twm = mu8.Tw;

% Define Model Parameters:
C = 1.5833;         % Pac. Tire Hyperparam.
B = -15.0975;       % Pac. Tire Hyperparam.
E = 0.6099;         % Pac. Tire Hyperparam.
r_e = 0.4013;       % Effective Tire Radius [m]; 0.37338;
J = 2.5462;         % Wheel Rotational Inertia [kg-m^2]
m = 2714.3;         % Vehicle Mass [kg]
Fz = m*9.81/4;      % Tire Normal Force [N]

%% Simulate single wheel model with Euler Integration

% Define helper functions:
calc_slip = @(w,U) r_e*w/U - 1;
w = @(s,U) (s + 1)*U/r_e;
get_force = @(U,w,mu) mu*Fz*sin(C*atan(B*(1 - E)*calc_slip(w,U)...
                      + E*atan(B*calc_slip(w,U))));

% Initialize Variables
t = [tm(1):1e-4:tm(end)]';
n = length(t);

U = zeros(n,1);
omega = U;
s = omega;
U(1) = Um(1);
s(1) = sm(1);
omega(1) = wm(1);
mu = 0.8;

% Euler Integration
for i = 2:n
    Fx(i) = -get_force(U(i-1),omega(i-1),mu);
    dU = Fx(i)/(m/4);
    U(i) = U(i-1) + dU*(t(i) - t(i-1));
    
    Tw = interp1(tm,Twm,t(i));
    Tb = interp1(tm,Tbm,t(i));
    torque = Tw - Tb;
    domega = (torque - r_e*Fx(i))/J;
    omega(i) = omega(i-1) + domega*(t(i) - t(i-1));
    s(i) = calc_slip(omega(i),U(i));
end

% Plot trajectories
figure();subplot(2,1,1);
plot(t,U,tm,Um); title('Euler Long. Vel.'); legend('Euler','Measured');
subplot(2,1,2);
plot(t,omega,tm,wm); title('Euler Long. Slip'); legend('Euler','Measured');

% Plot Errors:
for i = 1:length(tm)
    U_eul = interp1(t,U,tm(i));
    err_U(i) = Um(i) - U_eul;
    w_eul = interp1(t,omega,tm(i));
    err_w(i) = wm(i) - w_eul;
end
figure();subplot(2,1,1);
histogram(err_U,'Normalization','pdf');
title('Longitudinal Velocity [m/s]');
subplot(2,1,2);
histogram(err_w,'Normalization','pdf');
title('Wheel Angular Velocity [m/s]');
sgtitle('Euler Model Errors (pdf)');

%% Simulate single wheel model with ODE45 Integration
model_param = struct('C',C,'B',B,'E',E,'r_e',r_e,...
                     'J',J,'m',m,'Fz',Fz,'mu',0.8);
torque = Twm - Tbm;
inputs = struct('time',tm,'torque',torque);

% Initialization:
tspan = [tm(1), tm(end)];
y0 = [Um(1), wm(1)];

% Simulation:
options = odeset('RelTol',1e-8);
[t,y] = ode45(@(t,y) wheelode(t,y,model_param,inputs), tspan, y0, options);

% Get extra outputs:
for i = 1:length(t)
    [dy(i,:), ext(i,:)] = wheelode(t(i),y(i,:),model_param,inputs);
end

% Unpack outputs:
U = y(:,1);
omega = y(:,2);
T = ext(:,1);

% Plot Errors:
for i = 1:length(tm)
    U_45 = interp1(t,U,tm(i));
    err_U(i) = Um(i) - U_45;
    w_45 = interp1(t,omega,tm(i));
    err_w(i) = wm(i) - w_45;
end
figure();subplot(2,1,1);
histogram(err_U,'Normalization','pdf');
title('Longitudinal Velocity [m/s]');
subplot(2,1,2);
histogram(err_w,'Normalization','pdf');
title('Wheel Angular Velocity [m/s]');
sgtitle('ODE45 Model Errors (pdf)');