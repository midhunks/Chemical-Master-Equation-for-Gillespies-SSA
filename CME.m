%% Set up biochemical system
Input_Data

%% Adding subdirectories to the search path
if ispc, b='\'; else, b='/'; end % defining forward/bckward slashes
addpath([pwd,b,'State_Space_Builder',b]);
addpath([pwd,b,'Dynamics_Builder',b]);

%% Identify State space of the CME

S = State_Builder(Stoichiometry);

%% Identify the transition matrix of the CME

[D, State_Transition_Index_Matrix,SSA_propensity_matrix] = Dynamics_Builder(S, Stoichiometry, Reactants_stoichiometry);
