function [global_CE, regional_CE] = time_resolved_control_energy(Anorm, T, B, TS)
    %Anorm - normalized adjacency matrix
    %T - time-horizon
    %B - control strategy
    %TS - time-series overwhich to compute time-resolved control energy
    %TS should be size nparc x volumes
    % S. Parker Singleton, 2023
    
    transitions = size(TS,2)-1;
    nparc = size(TS,1);
    
    % Preallocate matrices for speed
    regional_CE = zeros(1,transitions,nparc); 
    global_CE = zeros(1,transitions);
    
    %final and initial states are adjacent volumes
    x0 = TS(:,1:size(TS,2)-1);
    xf = TS(:,2:size(TS,2));
    

    % Loop over each transition
    for transition = 1:transitions
        [~, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        regional_CE(1,transition,:) = sum(u.^2)*T/1001; % integrate over inputs for each region
        global_CE(1,transition) = sum(sum(u.^2))*T/1001; % integrate over regions
    end
end
