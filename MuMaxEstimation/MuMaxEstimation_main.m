% MAE 298: Estimation Final Project
% Authors:
% Li, Yihui
% Vidano, Trevor
% Zhong, Anna
%
% Extended Kalman Filter and Unscented Kalman Filter Applied to estimate
% the maximum tire road friction coefficient.

clear; close all; clc; 

%% Model Setup
% Define model parameters:
model_param.C = 1.5833;         % Pac. Tire Hyperparam.
model_param.B = -15.0975;       % Pac. Tire Hyperparam.
model_param.E = 0.6099;         % Pac. Tire Hyperparam.
model_param.r_e = 0.4013;       % Effective Tire Radius [m]; 0.37338;
model_param.J = 2.5462;         % Wheel Rotational Inertia [kg-m^2]
model_param.m = 2714.3;         % Vehicle Mass [kg]
Fz = model_param.m*9.81/4;      % Tire Normal Force [N]
model_param.Fz = 1.4*Fz;

model_param.Q = diag([1e-4,1e-4,1e-3]);%diag([3.1093,274.5482,0])
model_param.R = diag([1e-4,1e-1]);
model_param.N = 3;
model_param.M = 2;
model_param.ts = 2e-3;

% -------------------------------------------------------------------------
% Select which mu to load measurements:
mu = 0.80;
% mu = 0.50;
% mu = 0.30;
% -------------------------------------------------------------------------

% Collect measurement data:

% To use measurements from high fidelity model:
muData = matfile("mu" + num2str(mu,'%.2f') + ".mat");
Tb = muData.Tb;    % brake torque
Tw = muData.Tw;    % wheel torque (accel.)
torque = - Tb; 

% To use measurements from low fidelity model:
% muData = matfile("LF_mu" + num2str(mu,'%.2f') + ".mat");
% torque = muData.T;

% Measurement data from either model:
t = muData.t;      % time
U = muData.U;      % longitudinal speed
s = muData.s;      % long. tire slip
w = muData.w;      % wheel omega (ang. vel.)

% Store torque table:
model_param.t = t;
model_param.torque = torque;
model_param.options = odeset('RelTol',1e-3);

%% EKF/UKF Setup: 

% FOR EKF: the following are constant for all k 
C_pr = [1 0 0; 0 1 0];
E_pr = eye(3);
F_pr = eye(2);

% INITIAL values
states_ukf(:,1) = [U(1),w(1),.7]';  % COL VEC
var_ukf(:,:,1) = zeros(3,3);  % N x N MATRIX

states_ekf(:,1) = [U(1),w(1),mu]';
var_ekf(:,:,1) = zeros(3,3);

%% Implement UKF, EKF
for k = 2:1:length(t)
    % Determine sampling interval
    j = k - 1;
    
    % Get measurement yk
    yk = [U(k),w(k)]';
    
    % model prediction step, UKF
    model_param.tk = t(k);
    model_param.tj = t(j);
    [xk_ukf,pk_ukf] = ukf_pred(model_param,states_ukf(:,j),var_ukf(:,:,j),...
                               torque(j),@wheel_state_eqn);
    
    % model prediction step, EKF
    [xk_ekf,pk_ekf] = ekf_pred(model_param,states_ekf(:,j),var_ekf(:,:,j),...
                                torque(j),E_pr,@wheel_state_eqn);
    
    
    % measurement update step, UKF
    [states_ukf(:,k),var_ukf(:,:,k)] = ukf_upd(model_param, xk_ukf,... 
                                pk_ukf, torque(j), yk, @wheel_output_eqn);
                            
    % measurement update step, EKF
    [states_ekf(:,k),var_ekf(:,:,k)] = ekf_upd(model_param, xk_ekf, pk_ekf,...
                                torque(k), yk, C_pr, F_pr, @wheel_output_eqn);
    
end

% Unpack States: UKF
U_ukf = states_ukf(1,:);
w_ukf = states_ukf(2,:);
mu_ukf = states_ukf(3,:);

s_ukf = model_param.r_e*w_ukf./U_ukf - 1;

% Unpack states: EKF
U_ekf = states_ekf(1,:);
w_ekf = states_ekf(2,:);
mu_ekf = states_ekf(3,:);

s_ekf = model_param.r_e*w_ekf./U_ekf - 1;

%% Data Visualization
figure();subplot(3,1,1);
plot(t,mu_ukf,t,mu*ones(length(t),1));ylabel('mu');
legend('UKF','Measurement');
subplot(3,1,2);
plot(t,U_ukf,t,U);ylabel('U');
legend('UKF','Measurement');
subplot(3,1,3);
plot(t,torque);ylabel('torque');

figure();
% plot(t,mu_ukf); 
plot(t,mu_ukf,t,mu_ekf,t,mu*ones(length(t),1)); 
xlabel('Time [s]'); ylabel('\mu_{max}');
legend('UKF','EKF','Measurement');
grid on;

figure();
% plot(t,U_ukf);
plot(t,U_ukf,t,U_ekf,t,U);
xlabel('Time [s]'); ylabel('U [m/s]');
legend('UKF','EKF','Measurement');

% Estimation Error
mu_err_ukf = mu - mu_ukf;
mu_err_ekf = mu - mu_ekf;
figure
plot(t,mu_err_ukf,t,mu_err_ekf);
title('Estimation Error');
xlabel('Time');
ylabel('Error');
legend('UKF','EKF');

% Estimation Error
mu_err = mu - mu_ukf; 
figure
plot(t,mu_err);
title('Estimation Error');
xlabel('Time');
ylabel('Error');

% Plot distribution of errors
% PDF of Estimation Error
intv = 0.05;
xvals = -7:intv:7;
yvals = normpdf(xvals,0,sqrt(var_ukf(3,3,end)));
% bins = 2 * xvals(end) / intv;
% newDat = histBins(mu_err, bins, xvals(end));
% newLen = length(xvals) - 1;

figure
% plot(xvals(1:newLen),newDat);
histogram(mu_err_ukf,'Normalization','pdf','DisplayStyle','stairs');
hold on
histogram(mu_err_ekf,'Normalization','pdf','DisplayStyle','stairs');
% plot(xvals(1:newLen),yvals(1:newLen));
plot(xvals(1:end-1),yvals(1:end-1));
legend('UKF','EKF','Theoretical PDF');
xlabel('Range');
ylabel('Frequency');
ylim([0 20]);
title('KF Estimation Error PDF');
grid on;

figure();
% plot(t,s_ukf);
plot(t,s_ukf,t,s_ekf,t,s);
xlabel('Time [s]'); ylabel('\kappa');
ylim([-16,1.1]);
legend('UKF','EKF','Measurement');
grid on;
