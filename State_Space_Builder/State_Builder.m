function [S] = State_Builder(Stoichiometry)
%% State space Builder 
    % This function generates the state space for the chemical master
    % equation using stoichiometry of the system

tic; global number_species number_reactions Initial_Molecular_Population...
            Conservation_Sum dim_cons nstates minimal_CME

[number_species,number_reactions] = size(Stoichiometry);

%% Setting conservation Law. 
% For a closed system there are a number of conservartion laws exists. For
% open system, there may be some but not neccessary. We identify the
% possible conservation lwas in both systems using the following function

Conservation = Conservation_Law(Stoichiometry);

Conservation_Sum = Initial_Molecular_Population*Conservation;
dim_cons = min(size(Conservation)); % Need to confirm and optimize

%% Upper bound of molecular population of each species
% for a closed system, depending upon the conservation law and 
% initial molecular population, there is an upper bound on each species. 
% In open systems, we need to define the upper bound explicitely for
% generating a finite CME. For closed system algorithm  will automatically
% define the upper bopund of each species.

UB = Species_Upper_Bound(Conservation);

%% Re-indexing : 
% While using combvec, we have to keep the vectors generated using combvec
% less exploding. In every iteration, we will generate a set 
% of column vectors using combvec and then remove some rows which are 
% not satisfying the conservation laws. If the first vectors are large 
% enough, then by iteratively removing the rows using conservation, 
% the newer vectors formed with olde combvec output will be less expoloding. 
% The following re-indexing will help it.

Independent_Species_Index = Combvec_Optimization_index(Conservation,UB);

Conservation = Conservation(Independent_Species_Index,:);
UB = UB(Independent_Species_Index);
%% Initial molecular population for generating a minimal CME
% In a biochemical system depending upon the initial molecular population,
% dimension of the system varies. However, we can generate a mimal CME so
% that each species will take part in at least one of the reaction 
% using a smaller amount of initial molecular population. 
% The floowing lines does this job by idetifying the species which are 
% neccessary for this to happen.

if minimal_CME ~= 0
    Neccessary_Species = Independent_Species_Index(1:dim_cons);
    fprintf('To Produce a minimal CME use small values for species %d\n',Neccessary_Species);

    % Re-indexed Initial_Molecular_Population (in-built)
    Initial_Molecular_Population = zeros(1,number_species);
    Initial_Molecular_Population(1:dim_cons) = max(Conservation);
    
    UB = Species_Upper_Bound(Conservation);
end

%% All possible combinations of species
% The following function identifies al possible states the system can reach
% by means of conservation and upper bound.
S = Combination_Finder(Conservation,UB);
nstates = size(S,1);

%% Error check
Error_Check

%% Re-indexing back so that state space S is in the same order as input given.
[~,Reverse_Index] = sort(Independent_Species_Index,'ascend');
S = S(:,Reverse_Index);
if minimal_CME ~= 0
    Initial_Molecular_Population = Initial_Molecular_Population(Reverse_Index);
end

%% Re-ndex back only if needed
% Conservation = Conservation(Reverse_Index,:);
% UB = UB(Reverse_Index);

%%
fprintf('\n%.0f states generated in %.2f seconds\n',nstates,toc);
end


