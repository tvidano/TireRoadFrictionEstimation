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

<<<<<<< HEAD
Q = diag([1e-4,1e-3,1e-6]);%diag([3.1093,274.5482,0])
model_param.Q = Q;
R = diag([1e-6,1e-3]);
model_param.R = R;
=======
model_param.Q = diag([1e-6,3.2e-1,1e-8]);%diag([3.1093,274.5482,0])
model_param.R = diag([1e-6,2.6e-5]);
>>>>>>> stash
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

% To use measurements from high fidelity model 
% without noise:
% muData = matfile("mu" + num2str(mu,'%.2f') + ".mat");
% with noise: 
muData = matfile("noisy_mu" + num2str(mu,'%.2f') + ".mat");

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
mu0 = 0.7;
states_ukf(:,1) = [U(1),w(1),mu0]';  % COL VEC
var_ukf(:,:,1) = zeros(3,3);  % N x N MATRIX

states_ekf(:,1) = [U(1),w(1),mu0]';
var_ekf(:,:,1) = zeros(3);

%% Implement UKF, EKF
for k = 2:1:length(t)
    % Determine sampling interval
    j = k - 1;
    
    % Get measurement yk
    yk(:,k) = [U(k),w(k)]';
    
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
                                pk_ukf, torque(j), yk(:,k), @wheel_output_eqn);
                            
    % measurement update step, EKF
    [states_ekf(:,k),var_ekf(:,:,k)] = ekf_upd(model_param, xk_ekf, pk_ekf,...
                                torque(k), yk(:,k), C_pr, F_pr, @wheel_output_eqn);

    if abs(s(k)) < 1e-2
%         states_ukf(3,k) = mu;
%         var_ukf(3,3,k) = 1e-6;
%         states_ekf(3,k) = mu;
%         var_ekf(3,3,k) = 1e-6;
    end
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

%% MATLAB UKF validation

% initializations
states_mat(:,1) = states_ukf(:,1);
var_mat(:,:,1) = var_ukf(:,:,1);

% create object UKF
init_guess = states_ukf(:,1);
% Construct the filter
ukf = unscentedKalmanFilter(...
    @wheel_state_eqn,... % State transition function
    @wheel_output_eqn,... % Measurement function
    init_guess,...
    'HasAdditiveMeasurementNoise',true);
ukf.MeasurementNoise = R;
ukf.ProcessNoise = Q;

for ii = 1:1:length(t)
     
    % Let ii denote the current time.
    %
    % Residuals (or innovations): Measured output - Predicted output
    res_matlab(:,ii) = yk(:,ii) - wheel_output_eqn(ukf.State,torque(ii),model_param); % ukf.State is x[k|k-1] at this point
    % Incorporate the measurements at time k into the state estimates by
    % using the "correct" command. This updates the State and StateCovariance
    % properties of the filter to contain x[k|k] and P[k|k]. These values
    % are also produced as the output of the "correct" command.
    [states_mat(:,ii), var_mat(:,:,ii)] = correct(ukf,yk(:,ii),torque(ii),model_param);
    % Predict the states at next time step, k+1. This updates the State and
    % StateCovariance properties of the filter to contain x[k+1|k] and
    % P[k+1|k]. These will be utilized by the filter at the next time step.
    predict(ukf,torque(ii),model_param);
    
    
end

% Unpack States: MATLAB UKF
U_mat = states_mat(1,:);
w_mat = states_mat(2,:);
mu_mat = states_mat(3,:);

s_mat = model_param.r_e*w_mat./U_mat - 1;

%% Data Visualization
figure();subplot(4,1,1);
plot(t,mu_ukf,t,mu_mat,t,mu*ones(length(t),1));ylabel('mu');ylim([0,1]);
legend('UKF','MATLAB','Measurement');
subplot(4,1,2);
plot(t,U_ukf,t,U_mat,t,U);ylabel('U');
legend('UKF','MATLAB','Measurement');
subplot(4,1,3);
plot(t,w_ukf,t,w);ylabel('\omega');
subplot(4,1,4);
plot(t,s_ukf,t,s_mat,t,s);ylabel('slip');ylim([-2,0.2]);
legend('UKF','MATLAB','Measurement');

figure();
% plot(t,mu_ekf,t,mu_ukf,t,mu*ones(length(t),1)); 
plot(t,mu_mat,t,mu_ukf,t,mu*ones(length(t),1)); 
xlabel('Time [s]'); ylabel('\mu_{max}');
% legend('EKF','UKF','Measurement');
legend('MATLAB','UKF','Measurement');
grid on;

figure();
% plot(t,U_ukf,t,U_ekf,t,U);
plot(t,U_ukf,t,U_mat,t,U);
xlabel('Time [s]'); ylabel('U [m/s]');
% legend('UKF','EKF','Measurement');
legend('UKF','MATLAB','Measurement');

% Estimation Error
mu_err_ukf = mu - mu_ukf;
% mu_err_ekf = mu - mu_ekf;
mu_err_mat = mu - mu_mat;
figure
% plot(t,mu_err_ukf,t,mu_err_ekf);
plot(t,mu_err_ukf,t,mu_err_mat);
title('Estimation Error');
xlabel('Time');
ylabel('Error');
% legend('UKF','EKF');
legend('UKF','MATLAB');

% % Estimation Error
% mu_err = mu - mu_ukf; 
% figure
% plot(t,mu_err);
% title('Estimation Error');
% xlabel('Time');
% ylabel('Error');

% Plot distribution of errors
% PDF of Estimation Error
intv = 0.05;
xvals = -7:intv:7;
yvals_ukf = normpdf(xvals,0,sqrt(var_ukf(3,3,end)));
yvals_mat = normpdf(xvals,0,sqrt(var_mat(3,3,end)));
% yvals_ekf = normpdf(xvals,0,sqrt(var_ekf(3,3,end)));
% bins = 2 * xvals(end) / intv;
% newDat = histBins(mu_err, bins, xvals(end));
% newLen = length(xvals) - 1;

figure
% plot(xvals(1:newLen),newDat);
histogram(mu_err_ukf,'Normalization','pdf','DisplayStyle','stairs');
hold on
% histogram(mu_err_ekf,'Normalization','pdf','DisplayStyle','stairs');
histogram(mu_err_mat,'Normalization','pdf','DisplayStyle','stairs');
% plot(xvals(1:newLen),yvals(1:newLen));
plot(xvals(1:end-1),yvals_ukf(1:end-1));
% plot(xvals(1:end-1),yvals_ekf(1:end-1));
plot(xvals(1:end-1),yvals_mat(1:end-1));
% legend('UKF','EKF','Theoretical UKF PDF','Theoretical EKF PDF');
legend('UKF','MATLAB','Theoretical UKF PDF','Theoretical MATLAB PDF');
xlabel('Range');
ylabel('Frequency');
ylim([0 20]);
title('KF Estimation Error PDF');
grid on;

figure();
% plot(t,s_ukf,t,s_ekf,t,s);
plot(t,s_ukf,t,s_mat,t,s);
xlabel('Time [s]'); ylabel('\kappa');
ylim([-16,1.1]);
% legend('UKF','EKF','Measurement');
legend('UKF','MATLAB','Measurement');
grid on;
