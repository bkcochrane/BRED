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
    pHigh = pMax;
    pLow = pMin;
    
    for iter = 1:200

        fluxErrorCurr = fluxError(p);
        
        if abs(fluxErrorCurr) / massFluxTarget < 1e-6 && abs(p - pPrev) / p < 1e-6
            return;
        end

        if strcmp(branch, 'supersonic')
            if fluxErrorCurr < 0
                pLow = p;
            else
                pHigh = p;
            end
        else
            if fluxErrorCurr < 0
                pHigh = p;
            else
                pLow = p;
            end
        end

        dP = p - pPrev;
        dMF = fluxErrorCurr - fluxErrorPrev;
        
        if abs(dMF) < 1e-12 || abs(dP) < 1e-12

            pNext = 0.5 * (pHigh + pLow);

        else

            pNext = p - fluxErrorCurr * (dP / dMF);

            if isnan(pNext) || pNext >= pHigh || pNext <= pLow
                pNext = 0.5 * (pHigh + pLow);
            end

        end

        if abs(pNext - p) < 1e-6 && (pNext >= pMax * 0.99 || pNext <= pMin * 1.01)
            p = pNext;
            return;
        end

        pPrev = p;
        p = pNext;
        fluxErrorPrev = fluxErrorCurr;

    end
    
    % Fallback if Secant did not converge
    try
        p = fzero(fluxError, [pMin, pMax]);
    catch
        p = Engine.Throat.p;
    end

end