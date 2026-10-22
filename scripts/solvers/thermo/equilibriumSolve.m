function [temp, nj, Propellant] = equilibriumSolve(mode, target, p, PropData, temp, nj, needMu, isFrozen)

    if nargin < 8; isFrozen = false; end
    
    % If frozen, calculate T and return (skip phase changes)
    if isFrozen
        [temp, nj, Propellant, ~] = findEquilibrium(mode, target, p, PropData, temp, nj, needMu, true);
        return;
    end

    % Candidate condensed species
    if isfield(PropData, 'CONDENSED_LIST') && ~isempty(PropData.CONDENSED_LIST)
        candidates = PropData.CONDENSED_LIST(:).';
    else
        % No condensed species exist in database
        PropData.ACTIVE_CONDENSED = [];
        [temp, nj, Propellant, ~] = findEquilibrium(mode, target, p, PropData, temp, nj, needMu, false);
        return;
    end

    % Logical mask over candidates
    isActive = nj(candidates) > 1e-10;

    hasDropped = false(size(candidates));
    
    % Outer iteration loop
    for pass = 1:20
        active = candidates(isActive);
        PropData.ACTIVE_CONDENSED = active;
        [temp, nj, Propellant, piI] = findEquilibrium(mode, target, p, PropData, temp, nj, needMu, false);

        % Drop any active condensed species the solve pushed to zero
        if any(isActive)
            toDrop = isActive & (nj(candidates) <= 1e-10);
            if any(toDrop)
                nj(candidates(toDrop)) = 0;
                isActive(toDrop) = false;
                hasDropped(toDrop) = true;
                continue; % Re-solve with the reduced set
            end
        end
    
        % If all candidates are already active, none are inactive to test
        if all(isActive)
            break;
        end

        % Test inactive candidates to see if any species wants to condense
        gRT = Propellant.gRT;
        inactiveIdx = find(~isActive);
        inactSpecies = candidates(inactiveIdx);
        
        % Vectorized Eq 3.7 check across all inactive candidates at once
        drive = (piI.' * PropData.A_IJ(:, inactSpecies)) - gRT(inactSpecies).';
        toAdd = drive > 1e-8 & ~hasDropped(inactiveIdx).';
        
        if any(toAdd)
            addSpecies = inactSpecies(toAdd);
            isActive(inactiveIdx(toAdd)) = true;
            nj(addSpecies) = 1e-6; % Prevents solver oscillation
            continue; % Re-solve with the enlarged set
        end

        break; % Stable equilibrium
    end

end