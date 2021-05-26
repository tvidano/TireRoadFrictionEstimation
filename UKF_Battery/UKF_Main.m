% Unscented KF Main --> using Nonlinear Battery Model

clear; close all; clc; 

load 'IV_data_nonlinear.mat'; 
load('OCV_table.mat');

%% Set Up

% Share copies of these constants amongst all functions
global Q R R0 del_t Cbat

global soc_intpts_OCV OCV_intpts     
%  NOTE: I don't like this, but don't want to 
% look up every time I call the function or pass in an entire array
% I can break the measurement update step up to make this easier to pass in
% as an input, but I was too lazy to do that tbh. 
% LMK if anyone has any proposals on how to deal with this. It's not a big
% deal. 

del_t = 0.1;   % sampling time, s
Rc = 0.015;
Cc = 2400;
Cbat = 5 * 3600;
R0 = 0.01;
tauc = Rc*Cc;

% Initial conditions / distributions
Q = 2.5E-7;  % process noise variance
R = 1E-4;    % sensor noice variance
SOC(1) = 1;  % initial value of SoC
OL(1) = SOC(1);  % open loop estimate, initial
Vc(1) = 0;

% since initial SoC is a known value, it is its own mean and var = 0
P_soc(1) = 0; 
N = 1;

%% Implement UKF

for k = 2:1:length(t)
    % model prediction step
    j = k - 1;
    [xk,pk] = ukf_pred(N,SOC(j),P_soc(j),I(j));
    
    % Get measurement yk
    Vc(k) = exp(-del_t/tauc)*Vc(j) + Rc*(1-exp(-del_t/tauc))*I(j);
    yk = V(k) + Vc(k);
    
    % measurement update step
    [SOC(k),P_soc(k)] = ukf_upd(N,xk,pk,I(j),yk);
    
end

%% Plots

% COmpare Estimate vs Actual Data
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