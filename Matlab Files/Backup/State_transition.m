function S_transition = State_transition(S,Stoichiometry)
global number_reactions number_species nstates
S_transition = zeros(nstates,number_species,number_reactions);
for i=1:number_reactions
    for j= 1: number_species
        S_transition(:,j,i) = S(:,j) - Stoichiometry (j,i);
    end
end
end