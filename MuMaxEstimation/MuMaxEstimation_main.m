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

% Collect measurement data:
mu8 = matfile('mu0.80.mat');
t = mu8.t;  % time
U = mu8.U;  % longitudinal speed
s = mu8.s;  % long. tire slip
T = mu8.T;  % torque [INPUT]
w = mu8.w;  % wheel omega (ang. vel.)

% Determine sampling interval
ts = abs(t(2) - t(1)); 

% FOR EKF: the following are constant for all k 
C_pr = [1 0 0; 0 1 0];
E_pr = eye(3);
F_pr = eye(2);

% INITIAL values


%% Implement UKF, EKF

for k = 2:1:length(t)
    
    % Get measurement yk
    Vc(k) = exp(-del_t/tauc)*Vc(j) + Rc*(1-exp(-del_t/tauc))*I(j);
    yk = V(k) + Vc(k);
    
    % model prediction step, UKF
    j = k - 1;
    [xk_ukf,pk_ukf] = ukf_pred(model_param,states_ukf(j),var_ukf(j),I(j),@wheel_state_eqn);
    
    % model prediction step, EKF
    [xk_ekf,pk_ekf] = ekf_pred(model_param,states_ekf(j),var_ekf(j),I(j),E_pr,@wheel_state_eqn,ts);
    
    
    % measurement update step, UKF
    [states_ukf(k),var_ukf(k)] = ukf_upd(model_param, xk_ukf, pk_ukf, I(j), yk,...
                                @batt_output_eqn);
                            
    % measurement update step, EKF
    [states_ekf(k),var_ekf(k)] = ekf_upd(model_param, xk_ekf, pk_ekf, I(k), yk, C_pr, F_pr, @output_eqn);
    
end


%% Data Visualization 

