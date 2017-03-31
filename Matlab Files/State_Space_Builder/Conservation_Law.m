%% Setting conservation Law
function Conservation = Conservation_Law (Stoichiometry)
global Conservation_temp
Conservation=rref(null(Stoichiometry','r')')';

if ~isempty(find(Conservation<0, 1)) && ~isempty(Conservation_temp)
    Conservation = Conservation_temp;
end
end