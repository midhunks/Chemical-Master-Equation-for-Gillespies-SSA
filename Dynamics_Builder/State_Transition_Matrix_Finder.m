function State_Transition_Index_Matrix = State_Transition_Matrix_Finder(S,Stoichiometry)
global nstates number_reactions
State_Transition_Index_Matrix = zeros(number_reactions,nstates);
for i = 1:number_reactions
    % Here the Logic is - What are the new states from the current state
    % due to different reactions. Another possible logic to code is -
    % From which all state the system can reach current state.
    B = bsxfun(@plus,S',Stoichiometry(:,i));    
    %% ismember function is (very) slow. which affects the speed of algorithm
    [~,State_Transition_Index_Matrix(i,:)] = ismember(B',S,'rows');
    
end
end
