%% Clearing Variables
clc;
clear;
% close all;
more off;
% profile off; profile viewer;
% warning off; 
format shortG; format compact % format longg %rat

homefolder = 'C:\Users\mkathana\Dropbox\Study\UWaterloo\Brian\MATLAB codes\Chemical-Master-Equation-for-Gillespies-SSA';
cd(homefolder); 
%% Step 1: Generate the biochemical model
% File which describe the Biochemical model

cd(homefolder);  cd('Models');  
% run('MM.m')
run('LinearModel_3Species.m')
% run('Gene_Toggle.m')
%% Step 2: Generate the CME

cd(homefolder);  cd('CME');
CME

%% Step 3: Generate sample paths via SSA

homefolder = 'C:\Users\mkathana\Dropbox\Study\UWaterloo\Brian\MATLAB codes\Chemical-Master-Equation-for-Gillespies-SSA';
cd(homefolder); cd('SSA')

SSA