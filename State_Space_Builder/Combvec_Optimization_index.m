function [Index] = Combvec_Optimization_index(Conservation,UB)
% Identify the index that combvec function generate a matrix is is minimal

global number_species dim_cons
Index = 1:number_species;
Logical_Conservation = (Conservation ~= 0); %Generate matrix with only ones and zeros

%% Re-ndexing as per the descending order of the upper bound
[~,m] = sort(UB,'descend');
Index = Index(m);
Logical_Conservation = Logical_Conservation(m,:);

%% Identify and remove the dependent sepcies
Dependent_Species_Index = find(sum(Logical_Conservation,2)~=1);
Logical_Conservation(Dependent_Species_Index,:)= [];
Index(Dependent_Species_Index)=[];

%% Identifying the independent rows with largest upperbound
for i = 1:dim_cons
    Temp = Logical_Conservation(i,:);    
    Temp_Repetition_index = ismember(Logical_Conservation, Temp,'rows');
    Temp_Repetition_index(i) = 0; % To avoid the removal of parent Temp

    Logical_Conservation(Temp_Repetition_index,:)= [];
    Index(Temp_Repetition_index) = [];
end

Index = [Index setdiff(1:number_species,Index)];
end