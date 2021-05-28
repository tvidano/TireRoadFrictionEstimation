% Comparison of Measured Data from Multi-Body Dynamic Simulation, single
% wheel model solved with ODE45, and single wheel model solved with euler
% integration.
clear;clc;close all;
%% Collect measurement data:
mu8 = matfile('mu0.80.mat');
tm = mu8.t;
Um = mu8.U;
sm = mu8.s;
Tm = mu8.T;
n = length(tm);

%% Simulate single wheel model with Euler Integration
% Define Model Parameters:
C = 1.5833;         % Pac. Tire Hyperparam.
B = -15.0975;       % Pac. Tire Hyperparam.
E = 0.6099;         % Pac. Tire Hyperparam.
r_e = 0.4013;       % Effective Tire Radius [m]; 0.37338;
J = 2.5462;         % Wheel Rotational Inertia [kg-m^2]
m = 2714.3;         % Vehicle Mass [kg]
Fz = m*9.81/4;      % Tire Normal Force [N]

% Define helper functions:
calc_slip = @(w,U) r_e*w/U - 1;
w = @(s,U) (s + 1)*U/r_e;
get_force = @(U,w,mu) mu*Fz*sin(C*atan(B*(1 - E)*calc_slip(w,U)...
                      + E*atan(B*calc_slip(w,U))));

% Initialize Variables
t = linspace(tm(1),tm(end),1e6);
U = zeros(n,1);
omega = U;
s = omega;
U(1) = Um(1);
s(1) = sm(1);
omega(1) = w(sm(1),Um(1));
mu = 0.8;

% Euler Integration
for i = 2:n
    Fx = get_force(U(i-1),omega(i-1),mu);
    dU = Fx/(m/4);
    U(i) = U(i-1) + dU*(t(i) - t(i-1));
    domega = (Tm(i-1) - r_e*Fx)/J;
    omega(i) = omega(i-1) + domega*(t(i) - t(i-1));
    s(i) = calc_slip(omega(i),U(i));
end

%% Simulate single wheel model with ODE45 Integration
model_param = struct('C',C,'B',B,'E',E,'r_e',r_e,'J',J,'m',m,'Fz',Fz);
inputs = struct('time',tm,'torque',Tm);

% Initialization:
tspan = [tm(1), tm(end)];
y0 = [Um(1), w(s(1),U(1))];

% Simulation:
options = odeset('RelTol',1e-2);
[t,y] = ode15s(@(t,y) wheelode(t,y,model_param,inputs,mu), tspan, y0, options);

% Plot trajectories
% figure();subplot(2,1,1);
% plot(tm,U,tm,Um); title('Long. Vel.'); legend('Euler','Measured');
% subplot(2,1,2);
% plot(tm,s,tm,sm); title('Long. Slip'); legend('Euler','Measured');