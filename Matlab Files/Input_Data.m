%% Inpute data
%% Import from SBML file
% model = sbmlimport('filename.xml')
% Species=model.Species
% Stoichiometry=full(getstoichmatrix(model)); %confirm the reverse reactions
% %in stoichiometry. If not add it mannually by taking 
% %negative of the corresponding column
% Reactants_stoichiometry=Stoichiometry;
% Reactants_stoichiometry(find(Stoichiometry>0))=0;
% for i=1:size(model.Species,1)
% Initial_Molecular_Population(i)=model.Species(i).InitialAmount*model.Species(i).Parent.Capacity*6.623e23;
% end

%% Stoichiometry (Michalis-Menten mechansism)
Stoichiometry=[-1,  1,  0,  0;
               -1,  1,  1, -1;
                1, -1, -1,  1;
                0,  0,  1, -1];

%% Identify reactants: 
% In some reactions some species present as reactant and a product. 
% For e.g., A+B -> A+C. where A is a reactant and a product. 
% Net gain of A is zero molecules. This will remove A from original 
% Stoichiometry matrix and cannot understand the reactants properly 

Reactants_stoichiometry = [-1,  0,  0,  0;
                           -1,  0,  0, -1;
                            0, -1, -1,  0;
                            0,  0,  0, -1];
                      
%% Reaction rates (equations or constants)
% Usage: mention Species by S(:,i) for ith species.  
%        Mention parameters by Parm(i) for the ith parameter in Parm.
global Reaction_rate_Function Parm
Parm = [1e6,5e5,5e0,2e0]; %Parameters of reaction_rate
Reaction_rate_Function = @(Parm, S)...%S is the species popultation vector
                           [Parm(1)*ones(size(S,1),1) ... 
                            Parm(2)*ones(size(S,1),1) ...
                            Parm(3)*ones(size(S,1),1) ...
                            Parm(4)*ones(size(S,1),1)];

%% Initial molecular population
global Initial_Molecular_Population
Initial_Molecular_Population=[10 10 0 0];

%% Boundary condtion on molecular population of each species for open systems
%Usage: mention species values as a vector. 
%If the reaction system is closed, then keep the boundary condition as 

% Boundary_condition=inf*ones(1,size(Stoichiometry,1))

%if few species are conserved and others not, then keep conserved species 
%"inf" in the corresponding cell where xaxct upper bound will be identified
%in state space builder function. For unconserved species, mention the
%upper bound values in the corresponding cells

Boundary_condition=[inf inf inf inf]; %or Boundary_condition=inf*ones(1,size(Stoichiometry,1))

%% Conservation Laws
%keep this untouched as long as conservation laws identified from rref has
% integer values. Otherwise modify manually with integer conservation
% values
Conservation_temp=[];

