%% Boundary Consitions
Final_time = 1e2; % add inf if final time is not neccassary for the model
Number_of_Sample_paths = 2e1; % Number of sample paths to be generated

%% 
Propensity_Sum = sum(SSA_propensity_matrix,2);
Propensity_CDF = cumsum(SSA_propensity_matrix,2)./Propensity_Sum;
Propensity_Sum_Reciprocal = 1./Propensity_Sum;

%% Preallocation
global nstates number_species
n=max(nstates,1e8); % minimum pre-allocation should be equal to number of states
State_Index_temp = zeros(n,1);      Time_temp = zeros(n,1);

%% SSA
tic;    Time_Full = [];     State_Full = [];
for k=1:Number_of_Sample_paths
    State_Index = State_Index_temp;
    Time = Time_temp;
    Rand_numbers = rand(n,2);
    
    Current_state_index = 1;    i = 1;     t = 0;
    while ~isempty(Current_state_index) && t <= Final_time
        %% Preallocation updation
        % When 99% preallocation used, then add more area to the dynamics
        if rem(i,round(n*.99))==0
            PreAllocate_SSA
        end
        %% Saving data for analysis
        State_Index(i) = Current_state_index;
        Time(i) = t;
        %% Updating the system
        i = i+1;
        % Next state's index
        Current_state_index = State_Transition_Index_Matrix(Current_state_index,...
            find(Propensity_CDF(Current_state_index,:) >= Rand_numbers(i,1),1));
        % Updating time
        t = t + Propensity_Sum_Reciprocal(Current_state_index)*log(1/Rand_numbers(i,2));
    end    
    % Remove the extra preallocated area
    Time(i:end)=[];
    State_Index(i:end)=[];
    % Identify the states from the index
    State = S(State_Index,:);
    
    % Keep data for analysis (seperate different sample paths with NaN)
    Time_Full = cat(1,[Time_Full;nan],Time);
    State_Full = cat(1,[State_Full;nan*ones(1,number_species)],State);    
    
    fprintf('A sample path is generated in %.2f seconds\n',toc)
end
fprintf('SSA generated %.0f sample paths in %.2f seconds\n',Number_of_Sample_paths,toc)

clear Rand_numbers State_Index Time State

%% Plotting
figure(1)
plot(Time_Full, State_Full)
State_legendInfo =  {'Species 1','Species 2','Species 3','Species 4'};
Legend = legend(State_legendInfo);
ax = gca; % current axes
% ax.Title.String = 'Dynamics';
ax.XLabel.String = 'time';
ax.YLabel.String = 'molecular population';
ax.XLim = [0, max(Time_Full)];
ax.YLim = [0, max(State_Full(:))];
Figure_Setup % Common settings
