function Conservation = Conservation_Law (Stoichiometry)
%% Setting conservation Law
    % Conservation law of a biochemical system can be idetified by finding
    % the nul space of stoichiometry. So we tries to find an 
    % nonnegative-integer conservation laws using the following function. 
    
global Conservation_temp
Conservation = rref(null(Stoichiometry','r')')';

%% A negative or noninteger null space is not useful to generate the statespace. 
    % So if the above line of code fails, we need to enter conservation 
    % law manually which will be  easy by using the output of the previous 
    % line of codes. Save it to the variable 'Conservation_temp' in the
    % inpput_data.m file. Then run the algorithm to generate the CME
    
if ~isempty(find(Conservation < 0, 1)) && ~isempty(Conservation_temp)
    Conservation = Conservation_temp;
end
end