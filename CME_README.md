# Chemical Master equation


### Input the biochemical system to generate the CME
#### (Michaelis Menten Demo: Input_data.m)
=======================================================

##### To generate a Minimal CME (Choose 0 or 1)
	minimal_CME = 0;
    
##### Input Stoichiometry as per the system (demo: Michaelis Menten)

	Stoichiometry = [-1,  1,  0,  0;
                     -1,  1,  1, -1; 
                      1, -1, -1,  1;
                      0,  0,  1, -1];

##### Identify reactants: 
By setting stoichiometry in terms of reactants In some reactions some species present as reactant and product. For e.g., A + B -> A + C. where A is a reactant and a product. Net gain of A is zero molecules. This will remove A from original Stoichiometry matrix and cannot idnetify the reactants properly which is necassary for computing reaction propensities in CME

  	Reactants_stoichiometry = [-1,  0,  0,  0;
                               -1,  0,  0, -1;
                                0, -1, -1,  0;
                                0,  0,  0, -1];
                      
##### Input Reaction rates ( 4 constants for MM mechanism)
	Reaction_rates = [1e6,5e5,5e0,2e0];

The following function will generate an array of 'Reaction_rates'. However, when the 'Reaction_rates' are changing with respect to the current population of species, edit the following function accordingly. Leave it commented if the rates are constants as before.

	% global Rate
	% Rate = Reaction_rate_Function(Reaction_rates,S);
  	

##### Input Initial molecular population (4 species for MM mechanism)
  	Initial_Molecular_Population = [1, 2, 3, 4];

##### Boundary condtion on molecular population of each species for open systems.
In an open system, there will be species with unbounded population. If there exists a conservation of some species, they will be aleady bounded. For unconserved species, the upper bound values should be defined so that a finite CME can be generated. (Default) boundary condition is set to 10 molecules each for (unconserved) species. Conserved species will accept values irrespective of the input

  	Boundary_condition = [10 10 10 10];

##### Conservation Laws: 
Keep this untouched as long as the conservation law identified from rref has integer values. Otherwise, modify manually with integer conservation values. First identify conservation law by running (uncomment the following line)

  	Conservation = Conservation_Law (Stoichiometry)

Scale the result so that all columns have integer values.

  	Conservation_temp = [];
