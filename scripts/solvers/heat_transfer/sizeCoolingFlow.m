function Engine = sizeCoolingFlow(Engine, PropData, method)

    if nargin < 3
        method = 'Bartz';
    end

    TEMP_LIMIT_WALL = 1000;
    COOLANT_MARGIN = 40;
    TOLERANCE = 0.5;
    MAX_ITER = 20;

    % Saturation limit, evaluated at the jacket pressure
    if strcmp(Engine.Cooling.Coolant.NAME, 'WATER')
        pJacket = Engine.Cooling.Coolant.P_SUPPLY;
    else
        pJacket = Engine.Chamber.P * 1.2;
    end
    if ~isfield(Engine.Cooling.Coolant, 'tSat') || pJacket / 1e5 > max(Engine.Cooling.Coolant.P_REF)
        tempSat = Inf;
        tempLimitSat = Inf;
    else
        tempSat = Engine.Cooling.Coolant.tSat(pJacket);
        tempLimitSat = tempSat - COOLANT_MARGIN;
    end

    % Coking / thermal decomposition limit, only relevant for regen
    tempLimitDegradation = Inf;
    if Engine.Cooling.REGEN_COOLING && isfield(Engine.Cooling.Coolant, 'TEMP_DEGRADE')
        tempLimitDegradation = Engine.Cooling.Coolant.TEMP_DEGRADE - COOLANT_MARGIN;
    end

    tempLimitCoolant = min(tempLimitSat, tempLimitDegradation);


    tempIn = Engine.Cooling.Coolant.T_SUPPLY;

    filmFrac = Engine.Cooling.FILM_COOLING_FRACTION;
    mdotFuelAvailable = Engine.Fuel.mdot * (1 - filmFrac);

    if Engine.Cooling.REGEN_COOLING
        Engine.Cooling.Coolant.mdot = mdotFuelAvailable;
        Engine = calculateHeatFlux(Engine, PropData, method);
        Engine.Cooling.COOLANT_LIMITED = max(Engine.Wall.tempWG) > TEMP_LIMIT_WALL || max(Engine.Wall.tempCoolant) > tempLimitCoolant;
        Engine.Cooling.bypassFraction = 0;
        return;
    end

    converged = false;
    atCap = false;

    mdot = Engine.Cooling.Coolant.mdot;
    mdotCap = Inf;

    for iter = 1:MAX_ITER

        Engine.Cooling.Coolant.mdot = mdot;
        Engine = calculateHeatFlux(Engine, PropData, method);

        tempCoolantOut = max(Engine.Wall.tempCoolant);
        tempWallMax = max(Engine.Wall.tempWG);

        % Signed errors: positive means the limit is exceeded
        if isinf(tempLimitCoolant)
            errCoolant = -Inf;
            errCoolantNorm = -Inf;
        else
            errCoolant = tempCoolantOut - tempLimitCoolant;
            errCoolantNorm = errCoolant / max(tempLimitCoolant - tempIn, 10);
        end

        errWall = tempWallMax - TEMP_LIMIT_WALL;
        errWallNorm = errWall / 100;

        worstErr = max(errCoolant, errWall);
        worstErrNorm = max(errCoolantNorm, errWallNorm);

        if abs(worstErr) <= TOLERANCE
            converged = true;
            break;
        end

        % Already at the flow cap and still over a limit: nothing more to give
        if worstErr > 0 && mdot >= mdotCap * (1 - 1e-9)
            atCap = true;
            break;
        end

        % More flow if over a limit, less if comfortably under
        ALPHA = 0.5;
        correction = max(-0.3, min(0.5, ALPHA * worstErrNorm));
        mdot = min(mdot * (1 + correction), mdotCap);

    end

    if ~converged && ~atCap
        fprintf('[Warning] Coolant flow did not converge in %d iterations (worst error %.1f K).\n', ...
            MAX_ITER, worstErr);
    end

    Engine.Cooling.Coolant.mdot = mdot;
    Engine.Cooling.COOLANT_LIMITED = atCap;

    Engine.Cooling.Coolant.tempSat = tempSat;
    if isinf(tempSat)
        Engine.Cooling.Coolant.tempMargin = Inf;
    else
        Engine.Cooling.Coolant.tempMargin = tempSat - max(Engine.Wall.tempCoolant);
    end

    Engine.Cooling.bypassFraction = 0;

    if atCap
        fprintf(['[Warning] Regen cooling insufficient at full fuel flow (%.3f kg/s): ', ...
                 'wall %.0f K, coolant out %.0f K. Revisit channel geometry, ', ...
                 'wall material or cooling strategy.\n'], ...
                 mdot, max(Engine.Wall.tempWG), max(Engine.Wall.tempCoolant));
        Engine.Cooling.bypassFraction = 0;
    else
        Engine.Cooling.bypassFraction = 1 - (mdot / mdotFuelAvailable);
        fprintf('Cooling sized: %.1f%% of fuel bypasses the jacket.\n', ...
            Engine.Cooling.bypassFraction * 100);
    end

end