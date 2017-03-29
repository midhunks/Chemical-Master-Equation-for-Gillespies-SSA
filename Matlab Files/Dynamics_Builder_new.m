function [D,State_Transition_Matrix,Reaction_propensity]=Dynamics_Builder_new(Stoichiometry,S,Reactants_stoichiometry)
global nstates Reaction_rate_Function Parm
% tic
number_reactions=size(Stoichiometry,2);
% Reaction_Matrix = zeros(number_reactions,nstates);
State_Transition_Matrix=zeros(number_reactions,nstates);
for i=1:number_reactions
    B=bsxfun(@plus,S',Stoichiometry(:,i));
    
    %[States_in_B_from_Reaction_i_occurs,~]=find(all(B'>=0,2));
    %Reaction_Matrix(i,States_in_B_from_Reaction_i_occurs)=i;
    
    [~,Common_States_exists_in_S_and_B]=ismember(B',S,'rows');
    State_Transition_Matrix(i,:)=Common_States_exists_in_S_and_B;
end
clear B Possible_reactions States_from_A
% fprintf('Reaction & States transition matrix generated in %s seconds\n',toc)

%% Another way (little slow, but with more information) which can be used
%  to sort states according to the reaction progress

% tic
% number_species=size(Stoichiometry,1);
% number_reactions=size(Stoichiometry,2);
% Reaction_State_Transition_Matrix=zeros(number_species,nstates,number_reactions);
% State_Transition_Matrix=zeros(number_reactions,nstates);
% for i=1:number_reactions
%     %transition of state space due to ith reaction (includes nonrealistic states)
%     Reaction_State_Transition_Matrix(:,:,i)=bsxfun(@plus,S',Stoichiometry(:,i));
%     %Position of realistic states in Reaction_State_Transition_Matrix due to ith reaction
%     Position_in_R= ismember(Reaction_State_Transition_Matrix(:,:,i)',S,'rows');
%     %Posiition of realistic reaction in S due to ith reaction
%     Position_in_S=find(ismember(S,Reaction_State_Transition_Matrix(:,:,i)','rows'));
%     %Every column represents Position of transitioned states from state
%     %with index of column
%     State_Transition_Matrix(i,Position_in_R)=Position_in_S;
%     %Every column represent the reaction number so that a realistic state
%     %transition occured
%     Reaction_Matrix(i,Position_in_R)=i;
% end
% fprintf('Reaction & States transition matrix generated in %s seconds\n',toc)


%% Combination of reactions
%(parallel pooling is implemented and can be used when number of reeaction
%is large)
% tic
Combination=zeros(number_reactions,nstates);
Combination_diag=zeros(nstates,number_reactions);
for i=1:number_reactions
    %ith reaction stoichiometry
    Reaction_i=Reactants_stoichiometry(:,i);%to avoid overhead in parallel pool
    %positions of reactants
    Reactant_Position=find(Reaction_i~=0);
    
    Comb=zeros(nstates,length(Reactant_Position));%To generate vectors for combination products
    Combination_diag(:,i)=prod(S(:,Reactant_Position),2);% diagonal combinations
    
    %To avoid negative values in using nchoosek
    for j=1:length(Reactant_Position)
        %Number of molecules (Reactants) needed for reaction i
        Reactant_change=-Reaction_i(Reactant_Position(j));
        %To avoids negative values in gammaln function
        S_non_zero_position=find(S(:,Reactant_Position(j))-Reactant_change+1>0);
        %Possible states that the reaction i occurs
        S_temp=S(S_non_zero_position,Reactant_Position(j));
        %Finds nchoosek(S,R): number of ways reaction i can occur with S
        %molecule under R reactant consumption
        %gammaln function does the nchoosek function without overflow
        Comb(S_non_zero_position,j)=round(exp(gammaln(S_temp+1)...
            -gammaln(Reactant_change+1)-gammaln(S_temp-Reactant_change+1)));
    end
    %Possible Combinations of the occurance of reaction i from all states
    Combination(i,:)=prod(Comb,2);
end
% fprintf('Combinations matrix generated in %s seconds\n',toc)
clear S_non_zero_position Comb

%% Transition Matrix Builder
Reaction_propensity=[];
Reaction_propensity_index=[];
for i=1:number_reactions
    %i
    %Index of states from where ith reaction can occur
    State_index=find(State_Transition_Matrix(i,:)~=0);
    %Transitioned states after reaction i
    Transition_state_index=State_Transition_Matrix(i,State_index);
    S_temp=S(Transition_state_index,:);
    Rate=Reaction_rate_Function(Parm,S_temp);
    
    Reaction_propensity_index=cat(2,Reaction_propensity_index,...
        [State_Transition_Matrix(i,State_index);State_index]);
    Reaction_propensity_temp=bsxfun(@times,Rate(:,i),Combination(i,State_index)');
    Reaction_propensity=cat(1,Reaction_propensity,Reaction_propensity_temp);
end

Reaction_propensity_index=[Reaction_propensity_index [1:nstates;1:nstates]];
Reaction_propensity_temp=-dot(Combination_diag',Reaction_rate_Function(Parm,S)');
Reaction_propensity=[Reaction_propensity;Reaction_propensity_temp'];

% nonzeros_in_D=length(Reaction_propensity);
D=sparse(Reaction_propensity_index(1,:),Reaction_propensity_index(2,:),...
    Reaction_propensity,nstates,nstates);

% clear Reaction_propensity_index...
%     Reaction_propensity...
%     Combination...
%     Combination_diag...
%     Reaction_Transition_Matrix...
%     State_Transition_Matrix...
%     Rate...
%     S_temp...
%     Reaction_propensity_temp...
%     State_index...
%     Transition_state_index...
%     i...
%     j
% 
% %%
% clear A B
% Index=zeros(1,nstates)
% global Initial_Molecular_Population
% A=ismember(Initial_Molecular_Population,S,'rows')
% length_B=1
% Index(1)=A
% while (length_B<nstates)
%     K=State_Transition_Matrix(:,A)
%     A=setdiff(K,Index)'
% %     A=A(A>0)
%     length_A=length(A);
%     Index(length_B+1:length_B+length_A)=A;
%     length_B=length_B+length_A;
%     pause
% end


end