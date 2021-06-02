% MAE 298: Estimation Final Project
% Authors:
% Li, Yihui
% Vidano, Trevor
% Zhong, Anna
%
% Extended Kalman Filter and Unscented Kalman Filter Applied to estimate
% the maximum tire road friction coefficient.

clear; close all; clc; 

%% Setup
% Define model parameters:
model_param.C = 1.5833;         % Pac. Tire Hyperparam.
model_param.B = -15.0975;       % Pac. Tire Hyperparam.
model_param.E = 0.6099;         % Pac. Tire Hyperparam.
model_param.r_e = 0.4013;       % Effective Tire Radius [m]; 0.37338;
model_param.J = 2.5462;         % Wheel Rotational Inertia [kg-m^2]
model_param.m = 2714.3;         % Vehicle Mass [kg]
model_param.Fz = model_param.m*9.81/4; % Tire Normal Force [N]

model_param.Q = eye(3);
model_param.R = eye(2)*0.01;
model_param.N = 3;
model_param.M = 2;

% -------------------------------------------------------------------------
% Select which mu to load measurements:
mu = 0.80;
% mu = 0.50;
% mu = 0.30;
% -------------------------------------------------------------------------

% Collect measurement data:
muData = matfile("mu" + num2str(mu,'%.2f') + ".mat");
t = muData.t;      % time
U = muData.U;      % longitudinal speed
s = muData.s;      % long. tire slip
Tb = muData.Tb;    % brake torque
Tw = muData.Tw;    % wheel torque (accel.)
w = muData.w;      % wheel omega (ang. vel.)
torque = Tw - Tb;

% Determine sampling interval
ts = abs(t(2) - t(1)); 
model_param.ts = ts;

% FOR EKF: the following are constant for all k 
C_pr = [1 0 0; 0 1 0];
E_pr = eye(3);
F_pr = eye(2);

% INITIAL values
states_ukf(:,1) = [U(1),w(1),mu]';
var_ukf(:,:,1) = eye(3);
%% Implement UKF, EKF

for k = 2:1:length(t)
    
    % Get measurement yk
    yk = [U(k),w(k)]';
    
    % model prediction step, UKF
    j = k - 1;
    [xk_ukf,pk_ukf] = ukf_pred(model_param,states_ukf(:,j),var_ukf(:,:,j),...
                               torque(j),@wheel_state_eqn);
    
    % model prediction step, EKF
%     [xk_ekf,pk_ekf] = ekf_pred(model_param,states_ekf(j),var_ekf(j),I(j),E_pr,@wheel_state_eqn,ts);
    
    
    % measurement update step, UKF
    [states_ukf(:,k),var_ukf(:,:,k)] = ukf_upd(model_param, xk_ukf,... 
                                pk_ukf, torque(j), yk, @wheel_output_eqn);
                            
    % measurement update step, EKF
%     [states_ekf(k),var_ekf(k)] = ekf_upd(model_param, xk_ekf, pk_ekf, I(k), yk, C_pr, F_pr, @output_eqn);
    
end

% Unpack States:
U_ukf = states_ukf(1,:);
w_ukf = states_ukf(2,:);
mu_ukf = states_ukf(3,:);

%% Data Visualization 

figure();
plot(t,mu_ukf,t,mu*ones(length(t),1)); 
xlabel('Time [s]'); ylabel('\mu_{max}');
legend('UKF','Measurement');

figure();
plot(t,U_ukf,t,U);
xlabel('Time [s]'); ylabel('U [m/s]');
legend('UKF','Measurement');

