%% Download/Update CME codes
% If git is not configured, download it from the following repository
% https://github.com/midhunks/Chemical-Master-Equation

cd(SSA_Folder);
if exist('Chemical-Master-Equation','dir') == 7
    cd('Chemical-Master-Equation')
    % Update the repository if needed (Uncomment the following line)
    system('git pull https://github.com/midhunks/Chemical-Master-Equation.git')
else    
    system('git clone https://github.com/midhunks/Chemical-Master-Equation.git')    
end
