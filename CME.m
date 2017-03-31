%% Set up biochemical system
clear variables; clc

Input_Data

%% Adding subdirectories to the search path
if ispc, b='\'; else, b='/'; end % defining forward/bckward slashes
addpath([pwd,b,'State_Space_Builder',b]);
addpath([pwd,b,'Dynamics_Builder',b]);

%% Identify State space of the CME
global number_species number_reactions Conservation
[number_species,number_reactions] = size(Stoichiometry);

[S,Conservation]=State_Builder(Stoichiometry);

global nstates
nstates = size(S,1)

%% Identify the transition matrix of the CME

[D]=Dynamics_Builder_new(Stoichiometry,S,Reactants_stoichiometry);
