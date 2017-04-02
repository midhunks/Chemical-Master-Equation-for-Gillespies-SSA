function [D, State_Transition_Index_Matrix,SSA_propensity_matrix] = Dynamics_Builder(S, Stoichiometry, Reactants_stoichiometry)
tic; global nstates number_reactions Reaction_rates
%% Identifying the possible transitions from each state due to all reactions
State_Transition_Index_Matrix = State_Transition_Matrix_Finder(S,Stoichiometry);
 
%% Possible combination of interactions between species in all reactions
%(parallel pooling is implemented and can be used when number of reeaction is large)
[Combination, Combination_diag] = Combination_reaction_Finder(S,Reactants_stoichiometry);
 
%% Generating off-diagonal propensitites and index 
% Reaction rates array (if not defined in the Input_data)
Rate = Reaction_rate_Function(Reaction_rates);%,S);
% Preallocation
Reaction_propensity = zeros(nstates*(1+number_reactions),1);
Reaction_propensity_index = zeros(nstates*(1+number_reactions),2);
SSA_propensity_matrix = zeros(nstates,number_reactions);
j = 0;
for i = 1:number_reactions
    % Index of possible states from where ith reaction can occur
    State_index = find(State_Transition_Index_Matrix(:,i) ~= 0);
    % Index of states ith reaction occured
    Transition_state_index = State_Transition_Index_Matrix(State_index,i);
     
    l = length(Transition_state_index);    
    Reaction_propensity_index(j+1:j+l,:) = [Transition_state_index, State_index];    
    Reaction_propensity(j+1:j+l) = bsxfun(@times,Rate(Transition_state_index,i),Combination(State_index,i));    
         
    SSA_propensity_matrix(State_index,i) = Reaction_propensity(j+1:j+l);
     
    j = j + l;
end
 
%% Generating diagonal propensitites and index to off diagonal propensities
Reaction_propensity_index(j+1:j+nstates,:) = [1:nstates;1:nstates]';
Reaction_propensity(j+1:j+nstates) = -dot(Combination_diag',Rate');
j = j + nstates;
 
%% Removing extra preallocated area
Reaction_propensity_index(j+1:end,:) = [];
Reaction_propensity(j+1:end) = [];
 
%% Generating Dynamics matrix
D = sparse(Reaction_propensity_index(:,1),Reaction_propensity_index(:,2),...
           Reaction_propensity,nstates,nstates);
% nonzeros_in_D = length(Reaction_propensity);
fprintf('\nDynamics matrix generated in %.2f seconds\n',toc)
end