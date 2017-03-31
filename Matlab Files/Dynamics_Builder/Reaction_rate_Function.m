function Rate = Reaction_rate_Function(rates, S)
n = size(S,1);
m= length(rates);
Rate = zeros(n,m);
for i=1:m
    Rate(:,i) =rates(i)*ones(n,1);
end
end