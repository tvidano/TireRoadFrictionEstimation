% Script to perform data processing for tire road friction estimation.
clear; clc; close all;

% Setup output directories:
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
relPath = fullfile('..','MuMaxEstimation','data');
dataPath = fullfile(pathstr,relPath);
if ~exist(dataPath, 'dir')
    mkdir(dataPath)
end

% Collect data from AEB with Warning mu=0.8,0.5,0.3:
mus = ["0.80","0.50","0.30"];
start_times = [15,15,20];
end_times = [20,22,27];
for i = 1:length(mus)
    muData = readtable("mu" + mus(i) + ".txt");
    t = muData{:,1};
    % Use only braking portions:
    iStart = find(t>start_times(i),1);
    iEnd = find(t>end_times(i),1);
    t = t(iStart:iEnd); % Start t at 0
    U = muData{iStart:iEnd,3};
    s = muData{iStart:iEnd,4};
    Fx = muData{iStart:iEnd,5};
    Tb = muData{iStart:iEnd,6};
    w = muData{iStart:iEnd,7};
    Tw = muData{iStart:iEnd,9};
    % Remove simulation artifacts:
    Tb(Tb > 1e4 | Tb < -1e4) = 0;
    Tw(Tw > 1e4 | Tw < -1e4) = 0;
    b = (1/30)*ones(1,30);
    Tw = filter(b,1,Tw);
    w(w < -40) = 0;
    % Plot signals:
    figure
    ax1 = subplot(4,1,1);
    plot(t,U); ylabel('Long. Velocity [m/s]');grid on;
    ax2 = subplot(4,1,2);
    plot(t,w); ylabel('\omega [rad/s]');grid on;
    ax3 = subplot(4,1,3);
    plot(t,Tb); ylabel('Torque [n-m]');grid on;
    ylim([-1000,5000]);
    ax4 = subplot(4,1,4);
    plot(t,s); ylabel('Long. Slip');grid on;
    ylim([-1.1,1.1]);
    linkaxes([ax1,ax2,ax3,ax4],'x');xlabel('t [s]'); 
    % Save signals without noise:
    save(fullfile(dataPath, "mu" + mus(i) + ".mat"),'t','U','s','Fx',...
    'Tb','w','Tw');
    % Add measurement noise:
    len = length(t);
    w_noise = normrnd(0,0.0834,[len,1]);
    T_noise = normrnd(0,1,[len,1]);
    U_noise = normrnd(0,0.2528,[len,1]);
    w = w + w_noise;
    Tb = Tb + T_noise;
    U = U + U_noise;
    % Save signals with noise:
    save(fullfile(dataPath, "noisy_mu" + mus(i) + ".mat"),'t','U','s','Fx',...
    'Tb','w','Tw');
end