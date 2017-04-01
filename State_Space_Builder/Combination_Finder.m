function S = Combination_Finder(Conservation,UB)
%% finding all combinations of linearly independent (n-k) species
% The function identifies al possible states the system can reach
% by means of conservation and upper bound.
global number_species Conservation_Sum
dim_Cons_Cols=size(Conservation,2);

%%
A = 0:UB(dim_Cons_Cols+1);
for i = dim_Cons_Cols+2:number_species
  A = combvec(A,0:UB(i));

  %Removal of states that not satisfying conservation law
  K = bsxfun(@minus,Conservation_Sum,A'*Conservation(dim_Cons_Cols+1:i,:));
  A(:,any(K<0,2)) = [];
end
% Finding linearly dependent k species combination by solving linear problem
S = [linsolve(Conservation(1:dim_Cons_Cols,:)',...
     bsxfun(@minus, Conservation_Sum',Conservation(dim_Cons_Cols+1:end,:)'*A))' A'];
end