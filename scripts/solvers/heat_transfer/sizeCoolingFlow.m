function Engine = sizeCoolingFlow(Engine, PropData, method)
    if nargin < 3
        method = 'Bartz';
    end

    MARGIN = 40; % K below design margin

   % Saturation / boiling limit
    P_MIN = Engine.Chamber.P * 1.2;
    if P_MIN > max(Engine.Cooling.Coolant.P_REF)
        tempLimitSat = Inf;
        tempSat = Inf;
    else
        tempSat = Engine.Cooling.Coolant.tSat(P_MIN);
        tempLimitSat = tempSat - MARGIN;
    end
    
    % Coking / decomposition limit
    tempLimitCoke = Inf;
    if contains(upper(PropData.FUEL), 'RP') || contains(upper(PropData.FUEL), 'KEROSENE')
        tempLimitCoke = 560 - MARGIN;
    elseif contains(upper(PropData.FUEL), 'CH4') || contains(upper(PropData.FUEL), 'METHANE')
        tempLimitCoke = 850 - MARGIN;
    end
    
    tempLimit = min(tempLimitSat, tempLimitCoke);
    tempIn = Engine.Cooling.Coolant.T_SUPPLY;
    TOLERANCE = 0.5;

    % Start by testing 100% of the nominal fuel flow
    if Engine.Cooling.REGEN_COOLING
        mdotNominal = Engine.Fuel.mdot;
    else
        mdotNominal = Engine.Cooling.Coolant.mdot;
    end

    Engine.Cooling.Coolant.mdot = mdotNominal;
    Engine = calculateHeatFlux(Engine, PropData, method);
    tempOut = max(Engine.Wall.tempCoolant);

    % Store saturation margin at nominal flow
    Engine.Cooling.Coolant.tempSat = tempSat;
    if isinf(tempSat)
        Engine.Cooling.Coolant.tempMargin = Inf;
    else
        Engine.Cooling.Coolant.tempMargin = tempSat - tempOut;
    end

    if Engine.Cooling.REGEN_COOLING && (tempOut <= tempLimit)
        return;
    end

    % Iterate to find the required coolant mass flow if current flow not sufficienet
    mdotNew = mdotNominal;
    for iter = 1:20
        mdot = mdotNew;
        Engine.Cooling.Coolant.mdot = mdot;
        Engine = calculateHeatFlux(Engine, PropData, method);
        tempOut = max(Engine.Wall.tempCoolant);

        tempError = tempOut - tempLimit;
        if abs(tempError) <= TOLERANCE
            break; 
        end

        alpha = 0.3;
        correction = alpha * (tempError / max(tempLimit - tempIn, 10));
        mdotNew = mdot * (1 + correction);
        
        % Never drop below nominal fuel flow if regen cooling is active
        if Engine.Cooling.REGEN_COOLING
            mdotNew = max(Engine.Fuel.mdot, min(3.0 * mdot, mdotNew));
        else
            mdotNew = max(0.1 * mdot, min(3.0 * mdot, mdotNew));
        end
    end

    % Update margin after loop
    if ~isinf(tempSat)
        Engine.Cooling.Coolant.tempMargin = tempSat - max(Engine.Wall.tempCoolant);
    end

    if ~Engine.Cooling.REGEN_COOLING
        Engine.Cooling.Coolant.mdot = mdot;
        return;
    end

    % Guard against cooling flow exceeding total engine mass flow
    if mdot >= Engine.Specs.mdot
        warning('Required cooling mass flow (%.3f kg/s) exceeds total engine mdot. Clamping.\n', mdot);
        mdot = Engine.Specs.mdot * 0.95;
    end

    MIN_OF_RATIO = 0.5;
    requiredOxMdot = Engine.Specs.mdot - mdot;
    requiredOFRatio = requiredOxMdot / mdot;

    if requiredOFRatio < MIN_OF_RATIO
        warning('Required O/F ratio %.3f is below the minimum allowed limit: O/F ratio clamped to %.3f.\n', requiredOFRatio, MIN_OF_RATIO);
        Engine.Specs.ofRatio = MIN_OF_RATIO;
    else
        Engine.Specs.ofRatio = requiredOFRatio;
    end

    Engine.Ox.mdot = Engine.Specs.mdot * Engine.Specs.ofRatio / (Engine.Specs.ofRatio + 1);
    Engine.Fuel.mdot = Engine.Specs.mdot / (Engine.Specs.ofRatio + 1);
    Engine.Cooling.Coolant.mdot = Engine.Fuel.mdot;

    % Re-run calculateHeatFlux if mdot was clamped
    Engine = calculateHeatFlux(Engine, PropData, method);
end