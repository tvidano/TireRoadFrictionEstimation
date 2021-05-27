% Unscented KF Main --> using Nonlinear Battery Model

clear; close all; clc; 
%% Set Up

% Share copies of these constants amongst all functions
%global Q R R0 del_t Cbat

%global soc_intpts_OCV OCV_intpts     
%  NOTE: I don't like this, but don't want to 
% look up every time I call the function or pass in an entire array
% I can break the measurement update step up to make this easier to pass in
% as an input, but I was too lazy to do that tbh. 
% LMK if anyone has any proposals on how to deal with this. It's not a big
% deal. 
nonlinIV = matfile('IV_data_nonlinear.mat');
OCVTable = matfile('OCV_table.mat');

model_param.del_t = 0.1;    % sampling time, s
model_param.Rc = 0.015;     % Resistance, ohm
model_param.Cc = 2400;      % Capacitance, Farad
model_param.Cbat = 5 * 3600;% Battery Capacity, mAh?
model_param.R0 = 0.01;      % Battery Resistance, ohm
model_param.N = 1;          % Number of states to estimate
model_param.Q = 2.5E-7;     % process noise variance
model_param.R = 1E-4;       % sensor noice variance
model_param.t = nonlinIV.t; % time vector
model_param.I = nonlinIV.I;
model_param.V = nonlinIV.V;
model_param.soc_intpts_OCV = OCVTable.soc_intpts_OCV;
model_param.OCV_intpts = OCVTable.OCV_intpts;

% Unpack used model parameters:
tauc = model_param.Rc*model_param.Cc;
del_t = model_param.del_t;
Rc = model_param.Rc;
N = model_param.N;
t = model_param.t;
I = model_param.I;
V = model_param.V;
SOC_act = nonlinIV.SOC_act;

% Initial conditions / distributions
Vc = zeros(1,length(t));
SOC = Vc;
SOC(1) = 1;                 % initial value of SoC
OL(1) = SOC(1);             % open loop estimate, initial
% since initial SoC is a known value, it is its own mean and var = 0
P_soc = Vc;

%% Implement UKF


for k = 2:1:length(t)
    % model prediction step
    j = k - 1;
    [xk,pk] = ukf_pred(model_param,SOC(j),P_soc(j),I(j));
    
    % Get measurement yk
    Vc(k) = exp(-del_t/tauc)*Vc(j) + Rc*(1-exp(-del_t/tauc))*I(j);
    yk = V(k) + Vc(k);
    
    % measurement update step
    [SOC(k),P_soc(k)] = ukf_upd(model_param,xk,pk,I(j),yk);
    
end

%% Plots

% Compare Estimate vs Actual Data
figure 
plot(t,SOC_act);
hold on
plot(t, SOC);
xlabel('time');
ylabel('State of Charge');
legend('actual','UKF estimate');

fprintf('Final converged variance P[k|k] is %f\n', P_soc(end));

% estimation error
ek = SOC_act.' - SOC;
figure
plot(t,ek);
title('Estimation Error');
xlabel('Time');
ylabel('Error');

% PDF of Estimation Error
intv = 0.00005;
xvals = -0.02:intv:0.02;
yvals = normpdf(xvals,0,sqrt(P_soc(end)));
bins = 2 * xvals(end) / intv;
newDat = histBins(ek, bins, xvals(end));
newLen = length(xvals) - 1;

figure
plot(xvals(1:newLen),newDat);
hold on
plot(xvals(1:newLen),yvals(1:newLen));
legend('Estimation Error Data','Theoretical PDF');
xlabel('Range');
ylabel('Frequency');
title('Estimation Error PDF');