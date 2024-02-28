function [global_CE,regional_CE,global_CE_mean,regional_CE_mean] = get_weighted_control_energy_matrix(Anorm,T,states,B,normalizationType)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Anorm - normalized adjacency matrix
    % T - time-horizon
    % states - states to compute control energy matrix from
    % normalizationType - option to normalize each time-point. 'L2' is 
    % most common. (default = 'none')
    % S. Parker Singleton, 2023
    
    verbose=1; %set to 0 to suppress caution statments
    
    if nargin < 3
        normalizationType = 'none';
    end
    
    if iscell(states) % Check if states is a cell array
        
        [numRows, numCols] = size(states);
        global_CE = cell(numRows, numCols);
        regional_CE = cell(numRows, numCols);
        global_CE_mean = cell(numRows, numCols);
        regional_CE_mean = cell(numRows, numCols);
        
       
        
        if ~iscell(Anorm)
            if verbose
                disp('Caution: Multiple state sets provided but only one Anorm. Calculating using the same adjacency matrix for each set.')
            end
            for i = 1:numRows
                for j = 1:numCols
                    [global_CE{i, j}, regional_CE{i, j}, global_CE_mean{i,j}, regional_CE_mean{i,j}] = computeControlEnergy(Anorm, T, states{i,j}, B, normalizationType);
                end
            end
        else
            for i = 1:numRows
                for j = 1:numCols
                    [global_CE{i, j}, regional_CE{i, j}, global_CE_mean{i,j}, regional_CE_mean{i,j}] = computeControlEnergy(Anorm{i,j}, T, states{i,j}, B, normalizationType);
                end
            end
        end
        
        
        
    else
        [global_CE, regional_CE, global_CE_mean, regional_CE_mean] = computeControlEnergy(Anorm, T, states, B, normalizationType);
    end
end

function [g_CE, r_CE, g_CE_mean, r_CE_mean] = computeControlEnergy(Anorm, T, states_single, B, normalizationType)

    
    if size(states_single, 1) ~= size(Anorm, 1)
        states_single = states_single'; % states should be of size nparc x num_states
    end
    
    nparc = size(Anorm,1);
    num_states = size(states_single,2);

    Xf_ind = repmat(1:num_states,[1 num_states]); % final state order
    Xo_ind = repelem(1:num_states,num_states); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix

    x0 = states_single(:,Xo_ind);
    xf = states_single(:,Xf_ind); % now each column of x0 and xf represent state transitions

        % Apply normalization based on normalizationType
    switch normalizationType
        case 'L2'
            x0 = L2MAGNITUDENORM(x0);
            xf = L2MAGNITUDENORM(xf);
        case 'DISTANCE' 
            [x0,xf] = DISTANCENORM(x0,xf);
        case 'RADIAL'
            [x0,xf] = RADIALNORM(x0,xf);
        case 'DOUBLE'
            [x0,xf] = DOUBLENORM(x0,xf);
        otherwise
            % No normalization applied
    end

    r = NaN(1,num_states^2,nparc);
    g = NaN(1,num_states^2);
    % Loop over each transition
    for transition = 1:num_states^2
        [~, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        r(1,transition,:) = sum(u.^2)*T/1001; % integrate over inputs for each region
        g(1,transition) = sum(sum(u.^2))*T/1001; % integrate over regions
    end
    
    g_CE = reshape(g,[num_states num_states])';
    r_CE = reshape(r,[nparc num_states num_states]);
    r_CE = permute(r_CE,[1,3,2]);
    
    g_CE_mean = mean(g_CE,'all');
    r_CE_mean = mean(squeeze(mean(r_CE,2)),2);

    

end
