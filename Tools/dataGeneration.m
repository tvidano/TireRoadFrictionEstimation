% Generate data from low fidelity model for testing UKF on mu max
% estimation.

clear;clc;close all;

%% Collect measurement data:
% -------------------------------------------------------------------------
% Select which mu to compare models with the measurements:
mu = 0.80;
% mu = 0.50;
% mu = 0.30;
% -------------------------------------------------------------------------

% Define Model Parameters:
C = 1.5833;         % Pac. Tire Hyperparam.
B = -15.0975;       % Pac. Tire Hyperparam.
E = 0.6099;         % Pac. Tire Hyperparam.
r_e = 0.4013;       % Effective Tire Radius [m]; 0.37338;
J = 2.5462;         % Wheel Rotational Inertia [kg-m^2]
m = 2714.3;         % Vehicle Mass [kg]
Fz = m*9.81/4;      % Tire Normal Force [N]

model_param = struct('C',C,'B',B,'E',E,'r_e',r_e,...
                     'J',J,'m',m,'Fz',Fz,'mu',mu);

% Define Torque Inputs:
t = 0:1e-2:10;
torque = zeros(length(t),1);
for i = 1:length(t)
    if t(i) < 3
        torque(i) = 0;
    else
        torque(i) = -4000;
    end
end
inputs = struct('time',t,'torque',torque);

%% Simulate single wheel model with ODE45 Integration
% Initialization:
tspan = [t(1), t(end)];
U0 = 27.0;
w0 = U0/r_e;
y0 = [U0;w0];

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

% Plot trajectories
figure();subplot(2,1,1);
plot(t,U); title('Euler Longitudinal Velocity');
xlabel('Time [s]'); ylabel('U [m/s]');
subplot(2,1,2);
plot(t,omega); title('Euler Angular Velocity'); 
xlabel('Time [s]'); ylabel('\omega [rad/s]');
sgtitle("R-K 45 Model Trajectories (pdf) mu = " + num2str(mu,'%.2f'));