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
start_times = [10,12.7,17.8];
for i = 1:length(mus)
    muData = readtable("mu" + mus(i) + ".txt");
    t = muData{:,1};
    iStart = find(t>start_times(i),1);
    t = t(iStart:end); % Start t at 0
    U = muData{iStart:end,3};
    s = muData{iStart:end,4};
    T = muData{iStart:end,6};
    save(fullfile(dataPath, "mu" + mus(i) + ".mat"),'t','U','s','T');
    figure();subplot(2,1,1);
    plot(t,s); ylabel('Long. Slip');
    ylim([-1.1,1.1]);
    subplot(2,1,2);
    plot(t,U); xlabel('t [s]'); ylabel('Long. Velocity [m/s]');
end