%% Clearing Variables
tic; clc; clear all; close all; more off; 
% profile off; profile viewer; 
% warning off;
format shortG; format compact % format longg %rat

%% Setting defaults for plots
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',3);
% set(0,'defaultAxesLineStyleOrder','remove')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesLineStyleOrder','-|:|-.|--');
% set(0,'DefaultLineLineWidth','remove')
set(0,'DefaultAxesXGrid','on');
set(0, 'DefaultAxesBox', 'on');
set(0,'DefaultAxesTickDir', 'in')

%% Color/Grayscale plot
% set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 0 0;0 1 0;0 1 0]);

% colormap(gray);
% set(0,'DefaultAxesColorOrder',1e-1*[1 1 1;2 2 2;2 2 2;3 3 3;4 4 4;5 5 5]);

% set(0,'defaultAxesColorOrder','remove')

%%
fprintf('Initialized (clearing data and closing figures) in %.2f seconds\n\n',toc)