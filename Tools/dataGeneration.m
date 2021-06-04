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
Fz = 1.5*Fz;

model_param = struct('C',C,'B',B,'E',E,'r_e',r_e,...
                     'J',J,'m',m,'Fz',Fz,'mu',mu);

% Define Torque Inputs:
t_torque = 0:2e-3:6;
torque = zeros(length(t_torque),1);
for i = 1:length(t_torque)
    if t_torque(i) < 1
        torque(i) = 0;
    else
        torque(i) = -4000;
    end
end
inputs = struct('time',t_torque,'torque',torque);
U0 = 27.0;
w0 = U0/r_e;

%% Simulate single wheel model with ODE45 Integration
% Initialization:
tspan = t_torque(1):2e-3:t_torque(end);
y0 = [U0;w0];

% Simulation:
options = odeset('RelTol',1e-12);
[t,y] = ode45(@(t,y) wheelode(t,y,model_param,inputs), tspan, y0, options);

% Get extra outputs:
for i = 1:length(t)
    [dy(i,:), ext(i,:)] = wheelode(t(i),y(i,:),model_param,inputs);
end

% Unpack outputs:
U = y(:,1); U(isnan(U))=0;
w = y(:,2); w(isnan(w))=0;
T = ext(:,1);
s = r_e*w./U - 1;

%% Simulate single wheel model with Euler Integration

% % Define helper functions:
% calc_slip = @(w,U) r_e*w/U - 1;
% get_force = @(U,w,mu) mu*Fz*sin(C*atan(B*(1 - E)*calc_slip(w,U)...
%                       + E*atan(B*calc_slip(w,U))));
% 
% % Initialize Variables
% t = 0:1e-6:10;
% n = length(t);
% 
% U = zeros(n,1);
% omega = U;
% s = omega;
% U(1) = U0;
% omega(1) = w0;
% 
% % Euler Integration
% for i = 2:n
%     % Prevent divide by 0 issues:
%     if U(i) ~= 0 
%         Fx(i) = -get_force(U(i-1),omega(i-1),mu);
%     else 
%         Fx(i) = 0;
%     end
%     
%     dU = Fx(i)/(m/4);
%     U(i) = U(i-1) + dU*(t(i) - t(i-1));
%     
%     % Torque Input
%     if t(i) < 3
%         T(i) = 0;
%     else
%         T(i) = -4000;
%     end
%     
%     domega = (T(i) - r_e*Fx(i))/J;
%     omega(i) = omega(i-1) + domega*(t(i) - t(i-1));
%     
%     % Prevent the vehicle from traveling backwards: 
%     if U(i) <= 0
%         U(i) = 0;
%     end
%     % Prevent wheel from rotating backwards:
%     if omega(i) <= 0 
%         omega(i) = 0;
%     end
% 
%     s(i) = calc_slip(omega(i),U(i));
% end

%% Plot trajectories
figure();subplot(3,1,1);
plot(t,U); title('R-K 4,5 Longitudinal Velocity');
xlabel('Time [s]'); ylabel('U [m/s]');
subplot(3,1,2);
plot(t,w); title('R-K 4,5 Angular Velocity'); 
xlabel('Time [s]'); ylabel('\omega [rad/s]');
sgtitle("R-K 45 Model Trajectories (pdf) mu = " + num2str(mu,'%.2f'));
subplot(3,1,3);
plot(t,s); title('R-K 4,5 Longitudinal Slip');
xlabel('Time [s]'); ylabel('S');

%% Save output:
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
relPath = fullfile('..','MuMaxEstimation','data');
dataPath = fullfile(pathstr,relPath);
if ~exist(dataPath, 'dir')
    mkdir(dataPath)
end
save(fullfile(dataPath, "LF_mu" + num2str(mu,'%.2f') + ".mat"),'t','U','s','w','T');