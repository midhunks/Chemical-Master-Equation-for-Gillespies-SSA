function Minimal_CME(Independent_Species_Index)
%% Identify the species that are neccessary to be present initially 
  % in the system so that all reaction can occur
  
global minimal_set number_species Initial_Molecular_Population dim_cons
Neccessary_Species = 1:number_species;
Neccessary_Species = Neccessary_Species(Independent_Species_Index(1:dim_cons));
fprintf('To Produce a minimal CME use small values for species %d\n',Neccessary_Species)

if minimal_set ~= 0
  Initial_Molecular_Population = zeros(1,number_species);
  Initial_Molecular_Population(Neccessary_Species) = max(Conservation);
end

UB = Species_Upper_Bound(Conservation,Boundary_condition);
end