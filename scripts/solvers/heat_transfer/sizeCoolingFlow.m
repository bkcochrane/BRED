function Engine = sizeCoolingFlow(Engine, PropData)

    MARGIN = 40; % K below design margin

    % Test if boils
    P_MIN = Engine.Chamber.P * 1.2;
    if P_MIN > max(Engine.Cooling.Coolant.P_REF)
        % Coolant is supercritical
        tempLimitSat = Inf;
        tempSat = Inf;
    else
        tempSat = Engine.Cooling.Coolant.tSat(P_MIN);
        tempLimitSat = tempSat - MARGIN;
    end
    
    % Test if cokes
    tempLimitCoke = Inf;
    switch upper(PropData.FUEL)
        case 'RP-1'
            tempLimitCoke = 560 - MARGIN;
        case {'METHANE', 'CH4'}
            tempLimitCoke = 850 - MARGIN;
    end
    
    tempLimit = min(tempLimitSat, tempLimitCoke);
    tempIn = Engine.Cooling.Coolant.T_SUPPLY;
    TOLERANCE = 0.5;

    mdotNew = Engine.Cooling.Coolant.mdot;
    
    % Damped iterative loop to find required coolant mass flow
    for iter = 1:20

        mdot = mdotNew;
        Engine.Cooling.Coolant.mdot = mdot;
        Engine = calculateHeatFlux(Engine, PropData);
        tempOut = max(Engine.Wall.tempCoolant);

        tempError = tempOut - tempLimit;
        if abs(tempError) <= TOLERANCE
            break; 
        end

        % Damped proportional-integral correction step
        alpha = 0.3; % Relaxation factor
        correction = alpha * (tempError / (tempLimit - tempIn));
        mdotNew = mdot * (1 + correction);
        
        % Bound mdot within physical limits
        mdotNew = max(0.1 * mdot, min(3.0 * mdot, mdotNew));
    end

    Engine.Cooling.Coolant.tempSat = tempSat;
    Engine.Cooling.Coolant.tempMargin = tempSat - max(Engine.Wall.tempCoolant);

    if ~Engine.Cooling.REGEN_COOLING
        Engine.Cooling.Coolant.mdot = mdot;
        return;
    end

    % Guard against cooling flow exceeding total engine mass flow
    if mdot >= Engine.Specs.mdot
        warning('Required cooling mass flow (%.3f kg/s) exceeds total engine mdot. Clamping.\n', mdot);
        mdot = Engine.Specs.mdot * 0.95;
    end

    MIN_OF_RATIO = 0.8;

    if mdot <= Engine.Fuel.mdot
        Engine.Cooling.Coolant.mdot = Engine.Fuel.mdot;
        return;
    else
        requiredOxMdot = Engine.Specs.mdot - mdot;
        requiredOFRatio = requiredOxMdot / mdot;

        if requiredOFRatio < MIN_OF_RATIO
            warning('Required O/F ratio %.3f is below the minimum allowed limit: O/F ratio clamped to %.3f.\n', requiredOFRatio, MIN_OF_RATIO);
            Engine.Specs.ofRatio = MIN_OF_RATIO;
        else
            Engine.Specs.ofRatio = requiredOFRatio;
        end
    end

    Engine.Ox.mdot = mdot * Engine.Specs.ofRatio / (Engine.Specs.ofRatio + 1);
    Engine.Fuel.mdot = mdot / (Engine.Specs.ofRatio + 1);
    Engine.Cooling.Coolant.mdot = Engine.Fuel.mdot;
    
end