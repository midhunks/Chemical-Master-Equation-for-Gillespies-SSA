%% State space Builder
function [S,varargout]=State_Builder(Stoichiometry,Conservation_temp,...
    Initial_Molecular_Population,Boundary_condition)
number_species=size(Stoichiometry,1);
%% Setting conservation Law
Conservation=rref(null(Stoichiometry','r')')';

dim_Cons_Cols=size(Conservation,2);
if ~isempty(find(Conservation<0, 1)) && ~isempty(Conservation_temp)
    Conservation=Conservation_temp;
end

%While using combvec, we have to keep the first vectors in combvec as less
%exploding. This re-indexing will help it
Number_Species_conservation=sum(Conservation~=0,1);
[Number_Species_conservation,m]=sort(Number_Species_conservation,'descend');
Conservation=Conservation(:,m);

Conservation_Sum = Initial_Molecular_Population*Conservation;
dim_cons=min(size(Conservation));

%% Upper bound of molecular population of each species
UB=zeros(1,number_species);
%checks whether there is any conservation law exists
if ~isempty(Conservation)
    %Open (Unconserved) Species' upper bound
    Open_species_index=find(ismember(Conservation,...
        zeros(1,size(Conservation,2)),'rows'));
    if ~isempty(Open_species_index)
        %upper bound of ith species is set from the boundary condition
        UB(Open_species_index)=Boundary_condition(Open_species_index);
    end
    %Closed (Conserved) Species' upper bound
    for i=1:number_species
        %Index of ith species presence in conservation matrix
        Species_i_in_conservation = find(Conservation(i,:)~=0);
        %checks where ith species is present in conservation law
        if ~isempty(Species_i_in_conservation)
            if ismember(i,Open_species_index)
                error('open species is identified as conserved which is wrong')
            end
            % upper bound of ith species from minimum of conservation sum where ith
            % species is present.
            UB(i)=min(Conservation_Sum(Species_i_in_conservation));
        end
    end
else
    %upper bound of all species is set from the boundary condition
    UB=Boundary_condition;
end

%% First re-indexing according to upper bound (largest first) on species population
%While using combvec, we have to keep the vectors in combvec as small as
%possible. This reindexing will make sure last n-k species has lower upper bound

%finding index and re-indexing UB in descending order
[UB,Index_1]=(sort(UB,'descend'));

Conservation=Conservation(Index_1,:);%Re_index(Conservation,Index_1,'r');

%To identify Species that are neccessary to be present initially so that
%all reaction can occur
Neccessary_Species=[1:number_species];
Neccessary_Species=Neccessary_Species(Index_1);

%% Second Re-indexing Conservation matrices and Upper bound
%  This re-indexing will make sure, first square matrix formed with
%  'dim_cons' rows is an identity (or diagonal) matrix which will help to identify the
%  dependent species
K_eye=eye(dim_cons);
K_temp=Conservation(1:dim_cons,:);
if ~isequal(K_temp,K_eye);
    Idx=zeros(1,dim_cons);
    for i=1:dim_cons
        Idx(i)=find(ismember(Conservation,K_eye(i,:),'rows'),1);
%         Idx(i)=Index_Identity_Vec(1);
        K_temp(i,:)=Conservation(Idx(i),:);
    end
    %     [K_temp,Idx]=intersect(Conservation,K_eye,'rows','stable');
    if ~isequal(K_eye,K_temp)
        error('Problem with re-indexing of conservation law')
    end
else
    Idx=(1:dim_cons);
end
Index_2=1:number_species;
Index_2(ismember(Index_2,Idx))=[];
Index_2=[Idx Index_2];

Conservation=Conservation(Index_2,:);%Re_index(Conservation,Index_2,'r');

UB=UB(Index_2);%Re_index(UB,Index_2,'c');

%% (optional) Initial molecular population of the biochemical system so that
%  dimension of CME is minimal - Comment the code if not needed

% Neccessary_Species=Neccessary_Species(Index_2);
% Neccessary_Species=Neccessary_Species(1:dim_cons);
% fprintf('To reduce the dimension use small values for species %d\n',Neccessary_Species)
% 
% global minimal_set
% if minimal_set~=0
%     Initial_Molecular_Population=zeros(1,number_species);
%     Initial_Molecular_Population(Neccessary_Species)=max(Conservation);
% end

%% All possible combinations of species
%finding all combinations of linearly independent (n-k) species
A=0:UB(dim_Cons_Cols+1);
for i=dim_Cons_Cols+2:number_species
    A=combvec(A,0:UB(i));
    %fprintf('size of all combination is: %d x %d \n',size(A,1),size(A,2));

    %Removal of states that not satisfying conservation law
    K=bsxfun(@minus,Conservation_Sum,A'*Conservation(dim_Cons_Cols+1:i,:));
    A(:,any(K<0,2))=[];    
    %fprintf('size of all combination is: %d x %d \n',size(A,1),size(A,2));
end
%fprintf('size of all combination is: %d x %d \n',size(A,1),size(A,2));

%finding linearly dependent k species combination by solving linear problem
S=[linsolve(Conservation(1:dim_Cons_Cols,:)' ,bsxfun(@minus,...
    Conservation_Sum',Conservation(dim_Cons_Cols+1:end,:)'*A))' A'];

clear A K Conservation_Sum
%% Error check 1 (minute error)
[a,~]=find(S<0);
if ~isempty(a)
    warning('found negative values in S and removed')
    S(a,:)=[];
end
%% Error check 2 (Serious error)
%checking for zero column in s so that each species is present in atleast
%one state
Zero_species=find(sum(S~=0,1)==0);
if ~isempty(Zero_species)
    fprintf('Species which are not present in state space: %d \n',Zero_species)
    error('Not all species states found. Check negative values in Conservation matrix or check whether intial population has all necessary species')
end

%% Re-indexing back (2nd then 1st) so that state space S is in the same order as input given.
[~,Rev_Index_2]=sort(Index_2,'ascend');
[~,Rev_Index_1]=sort(Index_1,'ascend');

S=S(:,Rev_Index_2);%Re_index(S,Rev_Index_2,'c');
S=S(:,Rev_Index_1);%Re_index(S,Rev_Index_1,'c');

% %Re-ndex back only if needed
% Conservation=Conservation(Rev_Index_2,:);%Re_index(Conservation,Rev_Index_2,'r');
% Conservation=Conservation(Rev_Index_1,:);%Re_index(Conservation,Rev_Index_1,'r');
% UB=UB(Rev_Index_2);%Re_index(UB,Rev_Index_2,'c');
% UB=UB(Rev_Index_1);%Re_index(UB,Rev_Index_1,'c');

% clear Index_2 Index_1 Rev_Index_2 Rev_Index_1 UB
end


