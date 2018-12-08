%% Set SSA and CME folder
SSA_Folder = uigetdir('','Choose SSA Folder');
CME_Folder = uigetdir('','Choose CME Folder');

% OR you can Change the path using the following line
% SSA_Folder = '\Chemical-Master-Equation-for-Gillespies-SSA';
% CME_Folder = '\Chemical-Master-Equation';

%% Generate CME:
% Model templates are available in folder once the CME codes are downloded.
% Create or Change MATLAB files as needed. You can also modify 'CME' file
% in the folder to avoid asking model and default folders

cd(CME_Folder);
CME

%%
% Propensity_Sum = sum(SSA_propensity_matrix,2);
% Propensity_Sum_Reciprocal = 1./Propensity_Sum;

Propensity_CDF = cumsum(SSA_propensity_matrix,2);
Propensity_CDF_Normalized = bsxfun(@rdivide, Propensity_CDF,Propensity_CDF(:,end));

%% Preallocation

State_Index = [State_Index; State_Index_temp];
Time = [Time; Time_temp];
Rand_numbers = [Rand_numbers; rand(n,2)];

global nstates
n=max(nstates,1e8); % minimum pre-allocation should be equal to number of states
State_Index_temp = zeros(n,1);      Time_temp = zeros(n,1);

%% SSA
TotalTime = 0;    Time_Full = [];     State_Full = [];
for k=1:Number_of_Sample_paths
    tic;
    State_Index = State_Index_temp;
    Time = Time_temp;
    Rand_numbers = rand(n,2);
    
    Current_state_index = 1;    i = 1;     t = 0;
    while ~isempty(Current_state_index) && t <= Final_time
        %         %% Preallocation updation
        %         % When 99% preallocation used, then add more area to the dynamics
        %         if rem(i,round(n*.99))==0
        %             PreAllocate_SSA
        %         end
        %% Saving data for analysis
        State_Index(i) = Current_state_index;
        Time(i) = t;
        %% Updating the system
        i = i+1;
        % Next state's index
        Current_state_index = State_Transition_Index_Matrix(Current_state_index,...
            find(Propensity_CDF_Normalized(Current_state_index,:) >= Rand_numbers(i,1),1));
        % Updating time
        t = t + 1/Propensity_CDF(Current_state_index,end)*log(1/Rand_numbers(i,2));
    end
    % Remove the extra preallocated area
    Time(i:end)=[];
    State_Index(i:end)=[];
    % Identify the states from the index
    State = S(State_Index,:);
    
    % Keep data for analysis (seperate different sample paths with NaN)
    Time_Full = cat(1,[Time_Full;nan],Time);
    State_Full = cat(1,[State_Full;nan*ones(1,size(S,2))],State);
    
    TotalTime = TotalTime + toc;
    fprintf('A sample path is generated in %.2f seconds\n',toc)
end
fprintf('SSA generated %.0f sample paths in %.2f seconds\n',Number_of_Sample_paths,TotalTime)

clear Rand_numbers State_Index Time State

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
