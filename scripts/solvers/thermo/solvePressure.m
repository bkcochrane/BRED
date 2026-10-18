function p = solvePressure(s0R, h0, Engine, PropData, temp, nj, mdot, area, pRef, branch, isFrozen)

    massFluxTarget = mdot / area;
    fluxError = @(pVal) findMassFlux(pVal, s0R, h0, PropData, temp, nj, isFrozen) - massFluxTarget;
    
    if strcmp(branch, 'supersonic')
        pMin = 10;
        pMax = Engine.Throat.p;
        if pRef >= pMax
            p = pMax * 0.96;
        else
            p = min(max(pRef, pMin), pMax * 0.999);
        end
    else
        pMin = Engine.Throat.p;
        pMax = Engine.Chamber.P;
        if pRef <= pMin
            p = pMin * 1.04;
        else
            p = min(max(pRef, pMin * 1.001), pMax * 0.999);
        end
    end
    
    
    % Secant loop starting from pRef

    pPrev = p * 0.99; 
    fluxErrorPrev = fluxError(pPrev);
    
    for iter = 1:100
        fluxErrorCurr = fluxError(p);
        
        if abs(fluxErrorCurr) / massFluxTarget < 1e-6
            return;
        end

        if abs(fluxErrorCurr) > (1.20 * abs(fluxErrorPrev)) && iter > 2
            
            p = pPrev + 0.5 * (p - pPrev); % Cut step in half and try again
            
            if abs(p - pPrev) < 1e-4
                return;
            end

        end
        
        % Gradient calculation
        dpActual = p - pPrev;
        
        % Check for stagnant pressure
        if abs(dpActual) < 1e-6
            break; % Bail to fzero if pressure stops changing entirely
        end
        
        dFluxErrorDp = (fluxErrorCurr - fluxErrorPrev) / dpActual;
        
        % Check for flat gradient
        if abs(dFluxErrorDp) < 1e-16
            break; % Bail to fzero if mass flux becomes completely decoupled from pressure
        end
        
        % Calculate & clamp the step
        step = - (fluxErrorCurr / dFluxErrorDp);
        
        % Never let the solver change the pressure by more than 20% in one step
        maxStep = 0.20 * p; 
        clampedStep = max(-maxStep, min(maxStep, step));
        
        pNext = p + clampedStep;
        
        % Bounds clamping
        if pNext >= pMax
            pNext = p + 0.9 * (pMax - p);
            if (pMax - p) < 1e-3
                break; % Bail if crushed completely against the ceiling
            end
        elseif pNext <= pMin
            pNext = p + 0.9 * (pMin - p); % Push 90% of the way to the floor
            if (p - pMin) < 1e-3
                break;
            end
        end
        
        pPrev = p;
        fluxErrorPrev = fluxErrorCurr;
        p = pNext;
        
        if iter == 100
            return;
        end

    end
    
    % Fallback if Secant did not converge
    try
        p = fzero(fluxError, [pMin, pMax]);
        disp('Hello')
    catch
        % Only triggered at stations right next to throat where massFluxTarget >= G_max
        p = Engine.Throat.p;
    end

end