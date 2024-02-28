function [global_CE, regional_CE] = time_resolved_control_energy_fast(Anorm, T, TS,outliers,keeptimeseries,normalizationType)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this version uses the controlability gramian to quickly calculate
    % control energies and can only be used in scenarios where the identity
    % matrix is used as the control strategy. Often will be scaled from the
    % slower calculations, the amount of which depends on how the B matrix 
    % is constructed in those cases (with 1's or 2's along the diagonal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Anorm - normalized adjacency matrix
    % T - time-horizon
    % TS - time-series overwhich to compute time-resolved control energy
    % outliers - optional vector to remove specific columns (default=[])
    % keeptimeseries - option to keep the full time-series of control
    % energy rather than average over time. (default=false)
    % normalizationType - option to normalize each time-point. 'L2' is 
    % most common. (default = 'none')
    % S. Parker Singleton, 2023
    
    verbose=0; %set to 0 to suppress caution statments
    
    % Check for the existence of outliers and keeptimeseries
    if nargin < 4
        outliers = [];
    end
    if nargin < 5
        keeptimeseries = false;
    end
    if nargin < 6
        normalizationType = 'none';
    end
    
    if iscell(TS) % Check if TS is a cell array
        
        [numRows, numCols] = size(TS);
        global_CE = cell(numRows, numCols);
        regional_CE = cell(numRows, numCols);
        
        % Check for consistency between TS and outliers if outliers is given
        if ~isempty(outliers) 
            if ~iscell(outliers) || ~isequal(size(TS), size(outliers))
                error('Provide a cell array for outliers that matches dimensions of TS.');
            end
            
            if ~iscell(Anorm)
                if verbose
                    disp('Caution: Multiple TS provided but only one Anorm. Calculating using the same adjacency matrix for each scan.')
                end
                for i = 1:numRows
                    for j = 1:numCols
                        [global_CE{i, j}, regional_CE{i, j}] = computeControlEnergy(Anorm, T, TS{i,j}, outliers{i,j},keeptimeseries,normalizationType);
                    end
                end
            else
                for i = 1:numRows
                    for j = 1:numCols
                        [global_CE{i, j}, regional_CE{i, j}] = computeControlEnergy(Anorm{i,j}, T, TS{i,j}, outliers{i,j},keeptimeseries,normalizationType);
                    end
                end
            end
            
        else
        
            if ~iscell(Anorm)
                if verbose
                    disp('Caution: Multiple TS provided but only one Anorm. Calculating using the same adjacency matrix for each scan.')
                end
                for i = 1:numRows
                    for j = 1:numCols
                        [global_CE{i, j}, regional_CE{i, j}] = computeControlEnergy(Anorm, T, TS{i,j}, outliers,keeptimeseries,normalizationType);
                    end
                end
            else
                for i = 1:numRows
                    for j = 1:numCols
                        [global_CE{i, j}, regional_CE{i, j}] = computeControlEnergy(Anorm{i,j}, T, TS{i,j}, outliers,keeptimeseries,normalizationType);
                    end
                end
            end
        
        end
        
    else
        [global_CE, regional_CE] = computeControlEnergy(Anorm, T, TS, outliers,keeptimeseries,normalizationType);
    end
end

function [g_CE, r_CE] = computeControlEnergy(Anorm, T, TS_single, outliers,keeptimeseries,normalizationType)

    if size(TS_single, 1) ~= size(Anorm, 1)
        TS_single = TS_single'; % TS should be of size nparc x time
    end

    x0 = TS_single(:, 1:size(TS_single, 2)-1);
    xf = TS_single(:, 2:size(TS_single, 2));
    
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

    WcI = GRAMIAN_FAST(Anorm, T);

    [g_CE, r_CE] = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T, false);

    if ~isempty(outliers)
        if ~islogical(outliers) && ~all(ismember(outliers, [0, 1]))
            error('outliers vector must be logical or contain only 1s and 0s.');
        end
        if size(TS_single, 2) ~= length(outliers)
            error('Size of outliers vector must match the second dimension of TS.');
        end
        
        outliers_x0 = outliers(1:end-1); %outliers corresponding to x0 transitions
        outliers_xf = outliers(2:end); %outliers corresponding to xf transitions
        combined_outliers = outliers_x0 | outliers_xf;
        
        g_CE(:, combined_outliers) = [];
        r_CE(:, combined_outliers) = [];
    end
    
    if ~keeptimeseries
        g_CE = mean(g_CE, 2);
        r_CE = mean(r_CE, 2);
    end
    
    
    
end
