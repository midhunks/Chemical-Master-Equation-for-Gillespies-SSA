%% Error check for Negative populations
[a,~] = find(S<0);
if ~isempty(a)
    error('found negative values in S and removed')
    % S(a,:) = [];
end

%% Error check for absence of species
%checking for zero column in S so that each species is present in atleast
%one state
Zero_species = find(sum(S ~= 0, 1) == 0);
if ~isempty(Zero_species)    
    fprintf('\n\nsNot all species states found. Check negative values in Conservation matrix\n');
    Conservation
    fprintf('\nOR check whether intial population has all necessary species\n');
    Initial_Molecular_Population
    error('\nSpecies %d are not present in state space\n\n',Zero_species)    
end
