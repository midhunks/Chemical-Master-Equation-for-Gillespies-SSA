function [Combination, Combination_diag] = Combination_reaction_Finder(S,Reactants_stoichiometry)
global nstates number_reactions
Combination = zeros(nstates,number_reactions);
Combination_diag = zeros(nstates,number_reactions);
for i = 1:number_reactions %(USE parfor for parallel pooling)
    %ith reaction stoichiometry
    Reaction_i = Reactants_stoichiometry(:,i);%to avoid overhead in parallel pool
    %positions of reactants
    Reactant_Position = find(Reaction_i ~= 0);
    
    % Dynamic matrix's diagonal element combinations
    Combination_diag(:,i) = prod(S(:,Reactant_Position),2);
    
    %To generate vectors for combination products
    Comb = zeros(nstates,length(Reactant_Position));
    %To avoid negative values in using nchoosek
    for j = 1:length(Reactant_Position)
        %Number of molecules (Reactants) needed for reaction i
        Reactant_change = -Reaction_i(Reactant_Position(j));
        %To avoids negative values in gammaln function
        S_non_negative_position = find(S(:,Reactant_Position(j)) - Reactant_change >= 0);
        %Possible states that the reaction i occurs
        S_temp = S(S_non_negative_position,Reactant_Position(j));
        %Finds nchoosek(S,R): number of ways reaction i can occur with S
        %molecule under R reactant consumption
        %gammaln function does the nchoosek function without overflow
        Comb(S_non_negative_position,j) = round(exp(gammaln(S_temp+1)...
            -gammaln(Reactant_change+1)-gammaln(S_temp-Reactant_change+1)));
    end
    %Possible Combinations of the occurance of reaction i from all states
    Combination(:,i) = prod(Comb,2);
end
end