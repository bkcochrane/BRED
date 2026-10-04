function [T, N_J, propellant] = equilibriumSolve(MODE, TARGET, P, GasData, T, N_J, NEED_MU, IS_FROZEN)
    
    if nargin < 8; IS_FROZEN = false; end
    
    % If frozen, immediately calculate T and return (skip phase changes)
    if IS_FROZEN
        [T, N_J, propellant, ~] = findEquilibrium(MODE, TARGET, P, GasData, T, N_J, NEED_MU, true);
        return;
    end

    % Candidate condensed species
    if isfield(GasData,'CONDENSED_LIST'), CANDIDATES = GasData.CONDENSED_LIST(:).';
    else, CANDIDATES = [];
    end
 
    % Start with active set derived from the incoming guess
    active = CANDIDATES(N_J(CANDIDATES) > 1e-10);
    
    % Outer iteration loop
    for pass = 1:12
        GasData.ACTIVE_CONDENSED = active;
        [T, N_J, propellant, PI_I] = findEquilibrium(MODE, TARGET, P, GasData, T, N_J, NEED_MU);
 
        % Drop any active condensed species the solve pushed to zero
        drop = active(N_J(active) <= 1e-10);
        if ~isempty(drop)
            N_J(drop) = 0;
            active = setdiff(active, drop);
            continue; % Re-solve with the reduced set
        end
    
        % Test to see if species wants to condense
        [~, ~, ~, G_RT] = nasa9(GasData.A_NASA, GasData.B_NASA, GasData.T_HI, T);
        A_IJ = GasData.A_IJ;
        inactive = setdiff(CANDIDATES, active);
        add = [];

        for j = inactive
            % Eq 3.7, >0 means it wants to condense
            drive = A_IJ(:,j).' * PI_I - G_RT(j);
            if drive > 1e-8
                add(end+1) = j;
            end
        end

        if ~isempty(add)
            active = sort([active, add]);
            N_J(add) = 1e-6; % Prevents solver oscillation
            continue; % Re-solve with the enlarged set
        end
 
        break; % Stable equilibrium
        
    end
end