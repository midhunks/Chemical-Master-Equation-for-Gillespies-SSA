function [D] = Dynamics_Builder(S, Stoichiometry, Reactants_stoichiometry)
tic
global nstates Reaction_rates number_reactions

%% Identifying the possible transitions from each state due to all reactions
tic
State_Transition_Index_Matrix = State_Transition_Matrix_Finder(S,Stoichiometry);
toc
%% Possible combination of interactions between species in all reactions
%(parallel pooling is implemented and can be used when number of reeaction is large)
tic
[Combination, Combination_diag] = Combination_reaction_Finder(S,Reactants_stoichiometry);
toc
%% Transition Matrix Builder
tic
Reaction_propensity = [];
Reaction_propensity_index = [];
for i = 1:number_reactions
    % Index of possible states due to ith reaction 
    State_index = find(State_Transition_Index_Matrix(i,:) ~= 0);
    % Transitioned states after reaction i
    Transition_state_index = State_Transition_Index_Matrix(i,State_index);
    % Reaction rates
    Rate = Reaction_rate_Function(Reaction_rates,S(Transition_state_index,:));
    
    Reaction_propensity_index = cat(2,Reaction_propensity_index,...
                                      [Transition_state_index ; State_index]);
    Reaction_propensity_temp = bsxfun(@times,Rate(:,i),Combination(i,State_index)');
    Reaction_propensity = cat(1,Reaction_propensity,Reaction_propensity_temp);
end

Reaction_propensity_index = [Reaction_propensity_index, [1:nstates;1:nstates]];
Reaction_propensity_temp = -dot(Combination_diag',Reaction_rate_Function(Reaction_rates,S)');
Reaction_propensity = [Reaction_propensity;Reaction_propensity_temp'];

% nonzeros_in_D = length(Reaction_propensity);
D = sparse(Reaction_propensity_index(1,:),Reaction_propensity_index(2,:),...
    Reaction_propensity,nstates,nstates);

fprintf('\nTransition matrix generated in %.2f seconds\n',toc)
end