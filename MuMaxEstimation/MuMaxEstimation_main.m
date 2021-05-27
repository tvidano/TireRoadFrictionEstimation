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


%% Implement UKF

for k = 2:1:length(t)
    % model prediction step
    j = k - 1;
    [xk,pk] = ukf_pred(model_param,SOC(j),P_soc(j),I(j),@batt_state_eqn);
    
    % Get measurement yk
    Vc(k) = exp(-del_t/tauc)*Vc(j) + Rc*(1-exp(-del_t/tauc))*I(j);
    yk = V(k) + Vc(k);
    
    % measurement update step
    [SOC(k),P_soc(k)] = ukf_upd(model_param, xk, pk, I(j), yk,...
                                @batt_output_eqn);
    
end