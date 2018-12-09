%% Set SSA and CME folder
SSA_Folder = uigetdir('','Choose SSA Folder');
CME_Folder = uigetdir('','Choose CME Folder');

% OR you can Change the path using the following line
% SSA_Folder = '\Chemical-Master-Equation-for-Gillespies-SSA';
% CME_Folder = '\Chemical-Master-Equation';

%% Generate CME:
% Model templates are available in 'Models' folder once the CME codes are downloded.
% Create or Change MATLAB files as needed. You can also modify 'CME' file
% in the folder to avoid asking model and default folders

cd(CME_Folder);
CME

%%
Propensity_CDF = cumsum(SSA_propensity_matrix,2);
Propensity_CDF_Normalized = bsxfun(@rdivide, Propensity_CDF,Propensity_CDF(:,end));
% Propensity_CDF_Reciprocal = 1./Propensity_Sum;

%% SSA inputs
Final_time = input('Enter Final time:');
Number_of_Sample_paths = input('Enter number of sample paths to be generated:');

%% SSA
n=1e4; TotalTime = 0;    Time_Full = [];     State_Full = [];
for k=1:Number_of_Sample_paths
    tic;
    State_Index = nan(n,1);    Time = nan(n,1);
    
    Current_state_index = 1;    n = 1;     t = 0;
    while ~isempty(Current_state_index) && t <= Final_time        
        %% Saving data for analysis
        State_Index(n) = Current_state_index;
        Time(n) = t;
        
        %% Updating the system
        n = n+1;
        % Next state's index
        Current_state_index = State_Transition_Index_Matrix(Current_state_index,...
            find(Propensity_CDF_Normalized(Current_state_index,:) >= rand,1));
        % Updating time
        t = t + 1/Propensity_CDF(Current_state_index,end)*log(1/rand);
        
    end
    % Remove the extra preallocated area
    Time(n:end)=[];
    State_Index(n:end)=[];

    % Identify the states from the index
    State = State_Space(State_Index,:);
    
    % Keep data for analysis (seperate different sample paths with NaN)
    Time_Full = cat(1,[Time_Full;nan],Time);
    State_Full = cat(1,[State_Full;nan*ones(1,size(State_Space,2))],State);

    TotalTime = TotalTime + toc;
    fprintf('A sample path is generated in %.2f seconds\n',toc)
end
fprintf('SSA generated %.0f sample paths in %.2f seconds\n',Number_of_Sample_paths,TotalTime)

clearvars Rand_numbers State_Index Time State

%% Plotting
figure(1)
plot(Time_Full, State_Full)
State_legendInfo = [];
for i = 1:size(State_Full,2)
    State_legendInfo =  cat(2, State_legendInfo, strcat('Species',{' '}, num2str(i)));
end
Legend = legend(State_legendInfo);
ax = gca; % current axes
% ax.Title.String = 'Dynamics';
ax.XLabel.String = 'Time';
ax.YLabel.String = 'Species population';
ax.XLim = [0, max(Time_Full)];
ax.YLim = [0, max(State_Full(:))];
cd(SSA_Folder);     cd('Generalized_Functions');
Figure_Setup % Common settings
