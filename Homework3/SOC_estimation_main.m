%% Main MATLAB script for HW 3 in MAE 298: Estimation
% author: Trevor Vidano
% course: MAE 298: Estimation Spring Quarter 2021
% MATLAB implementation of kalman filter and extended kalman filter on
% battery equivalent circuit model. 

clearvars; close all;
%% Battery parameters
C_bat= 5*3600;  % A-s
R0 = 0.01; % ohm
Rc = 0.015; % ohm
Cc = 2400; % F
alpha = 0.65; % V
Vocv0 = 3.435; % V

Q = 2.5E-7; % Process noise covariance
R = 1.0E-4; % Sensor noise covariance
Ts = 0.1;  % sampling period

%%%%%%%%%%%%%%%%%%%%%%% State Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous time State Space Matrices
A_c = [-1/Cc/Rc,0;0,0]; 
B_c = [1/Cc;-1/C_bat]; 
C_c = [-1,alpha]; 
D_c = -R0;

% Discrete time State Space Matrics
F = expm(A_c*Ts); 
G = [Rc - Rc*exp(-Ts/Cc/Rc); -Ts/C_bat];
H = C_c; 
M = D_c;

% Helper functions:
%{
kalman_gain = @(P_priork,H_k,R_k) ...
    P_priork*H_k'*inv(H_k*P_priork*H_k' + R_k);
cov_prior = @(F_k0,P_postk0,Q_k0) ;
cov_post = @(H_k,K_k,P_priork,R_k,size) ...
    (eye(size) - K_k*H_k)*P_priork*(eye(size) - K-k*H_k)' + K_k*R_k*K_k';
prior_est = @(F_k0,x_postk0,G_k0,u_k0) F_k0*x_postk0 + G_k0*u_k0;
post_est = @(x_priork,K_k,y_k,H_k) x_priork + K_k*(y_k - H_k*x_priork);
%}
%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linearIV = matfile('IV_data_linear.mat');
u = linearIV.I';
t = linearIV.t';
SOC_act = linearIV.SOC_act';
V_measured = linearIV.V';

k_end = t(end)/Ts;

% scalar initilization:
% SOC_post = 1;
% SOC_open = SOC_post;
% Vc = 0;
% P_post = 0;
% matrix initialization:
x_post(:,1) = [0;1];
x_open(:,1) = x_post(:,1);
P_post{1} = eye(2)*0;
for k = 2:k_end
%     % Full system estimation:
%     % Model Predict:
%     P_prior{k} = F*P_post{k-1}*F' + [0,0;0,Q];
%     K(:,k) = P_prior{k}*H'*inv(H*P_prior{k}*H' + R);
%     x_prior(:,k) = F*x_post(:,k-1) + G*u(:,k-1);
%     % Measurement Update:
%     y(:,k) = V_measured(k) - Vocv0;
%     x_post(:,k) = x_prior(:,k) + K(:,k)*(y(:,k) - (H*x_prior(:,k) + M*u(:,k)));
% %     P_post{k} = (eye(2) - K(:,k)*H)*P_prior{k}*(eye(2) - K(:,k)*H)' + ...
% %                 K(:,k)*R*K(:,k)';
%     P_post{k} = P_prior{k} - P_prior{k}*H'*inv(H*P_prior{k}*H' + R)*H*P_prior{k};
%     x_open(:,k) = F*x_open(:,k-1) + G*u(k-1);
%     P(k) = P_post{k}(4);
    
    % Estimation of only SOC:
    % Model Predict:
    P_prior(k) = 1*P_post(k-1)*1' + Q;
    K(k) = P_prior(k)*alpha/(alpha*P_prior(k)*alpha + R);
    SOC_prior(k) = 1*SOC_post(k-1) - Ts/C_bat*u(k-1);
    Vc(k) = Vc(k-1)*exp(-Ts/Cc/Rc) + Rc*(1 - exp(-Ts/Cc/Rc))*u(k-1);
    % Measurement Update:
    y(k) = V_measured(k) - Vocv0 + Vc(k);
    SOC_post(k) = SOC_prior(k) + K(k)*(y(k) - (alpha*SOC_prior(k) - R0*u(k-1)));
    P_post(k) = (1 - K(k)*alpha)^2*P_prior(k) + K(k)^2*R;
    
    % Open loop estimate:
    SOC_open(k) = SOC_open(k-1) - Ts/C_bat*u(k-1);
end

%%%%%%%%%%%%%%%%%%%%%%% Generate Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
plot(t,SOC_post,t,SOC_act,t,SOC_open); 
legend('estimated','actual','open-loop');
ylim([0,1.0]);
figure();
plot(t,P_post);