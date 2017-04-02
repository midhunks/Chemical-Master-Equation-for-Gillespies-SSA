%% Identify State space of the CME 
S = State_Builder(Stoichiometry);

%% Identify the transition matrix of the CME
[D, State_Transition_Index_Matrix,SSA_propensity_matrix]...
        = Dynamics_Builder(S, Stoichiometry, Reactants_stoichiometry);