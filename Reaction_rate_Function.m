function Rate = Reaction_rate_Function(rates)%, S)
% S will be useful if the rates are changing in terms of population
global number_reactions nstates
Rate = zeros(nstates,number_reactions);
for i=1:number_reactions
    % ones(n,1) can be changed to match with a variable rate if neccessary.
    % That is why we are inputing here instead of the length. But this
    % needs to be done by user depending upon the system. Enter the rate as
    % a column vector by replacing ones(n,1)
    Rate(:,i) =rates(i)*ones(nstates,1);
end
end