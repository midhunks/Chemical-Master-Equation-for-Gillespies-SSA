function UB = Species_Upper_Bound(Conservation)
%% Upper bound of molecular population of each species
    % for a closed system, depending upon the conservation law and 
    % initial molecular population, there is an upper bound on each species. 
    % In open systems, we need to define the upper bound explicitely for
    % generating a finite CME. For closed system algorithm  will automatically
    % define the upper bopund of each species.

global number_species Boundary_condition Conservation_Sum

% Default boundary condition
if isempty(Boundary_condition)
    Boundary_condition = 10*ones(1,size(Stoichiometry,1));
end

Species_Index = 1:number_species;
UB = zeros(1,number_species);
%checks whether there is any conservation law exists
if ~isempty(Conservation)
    % Open (Unconserved) Species' upper bound
    Open_Species_index = find(sum(Conservation,2) == 0);
    if ~isempty(Open_Species_index)
        %upper bound of ith species is set from the boundary condition
        UB(Open_Species_index) = Boundary_condition(Open_Species_index);
    end
    
    %Closed (Conserved) Species' upper bound
    Closed_Species_index = setdiff(Species_Index,Open_Species_index);
    % upper bound of ith species from the minimum of conservation sum
    % where ith species is present.
    if ~isempty(Closed_Species_index)
        Conservation_temp = Conservation(Closed_Species_index,:);
        Conservation_temp(Conservation==0)= nan;
        UB(Closed_Species_index)=min(floor(Conservation_Sum./Conservation_temp)');
    end
else
    %upper bound of all species is set from the boundary condition
    UB = Boundary_condition;
end
%%
end