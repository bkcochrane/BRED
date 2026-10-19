function Engine = calculateProfile(Engine)

    volumeChamber = Engine.Geometry.LStar * Engine.Throat.area; % m3
    areaChamber = Engine.Throat.area * Engine.Geometry.CONTRACTION_RATIO_TARGET; % m2
    diameterChamber = 2 * sqrt(areaChamber/pi); % m
    diameterThroat = 2 * sqrt(Engine.Throat.area/pi); % m
    diameterExit = 2 * sqrt(Engine.Exit.area/pi); % m

    % ------------------------------------------------------------------------
    % Calculate converging section profile 
    % using two blends and a converging angle beta
    % ------------------------------------------------------------------------
    
    betaAngle = Engine.Geometry.BETA_ANGLE;
    radiusThroat = diameterThroat / 2;
    radiusChamber = diameterChamber / 2;
    
    throatFactor = Engine.Geometry.BLEND_RADIUS_THROAT_FACTOR;
    chamberFactor = Engine.Geometry.BLEND_RADIUS_CHAMBER_FACTOR;

    radialSpan = radiusChamber - radiusThroat;
    blendRadiusThroat = throatFactor * radiusThroat;
    if blendRadiusThroat * (1-cos(betaAngle)) > 0.45 * radialSpan
        blendRadiusThroat = 0.45 * radialSpan / (1-cos(betaAngle));
    end
    
    % Blend to throat
    theta = linspace(0, betaAngle, 20);
    xThroatArc = -blendRadiusThroat * sin(theta);
    rThroatArc = radiusThroat + blendRadiusThroat * (1 - cos(theta));
    x1 = xThroatArc(end);
    r1 = rThroatArc(end);

    contractionRatio = Engine.Geometry.CONTRACTION_RATIO_TARGET;
    contractionRatioFloor = 2;
    
    for outer = 1:20
        
        areaChamber = Engine.Throat.area * contractionRatio;
        diameterChamber = 2*sqrt(areaChamber/pi);
        radiusChamber = diameterChamber / 2;

        maxChamberDrop = 0.5 * radialSpan;
        blendRadiusChamber = chamberFactor * radiusChamber;
        blendDrop = blendRadiusChamber * (1 - cos(betaAngle));
        
        if blendDrop > maxChamberDrop
            blendRadiusChamber = maxChamberDrop / (1 - cos(betaAngle));
        end
        
        r3 = radiusChamber - blendRadiusChamber * (1 - cos(betaAngle));
        if r3 <= r1
            blendRadiusChamber = 0.8 * (radiusChamber - r1) / (1 - cos(betaAngle));
            r3 = radiusChamber - blendRadiusChamber * (1 - cos(betaAngle));
        end
        x2 = x1 - (r3 - r1) / tan(betaAngle);
        
        % Blend to chamber
        phi = linspace(betaAngle, 0, 10);
        xChamberCentre = x2 - blendRadiusChamber * sin(betaAngle);
        xChamberArc = xChamberCentre + blendRadiusChamber * sin(phi);
        rChamberArc = radiusChamber - blendRadiusChamber * (1 - cos(phi));

        numConePoints = 10;
        xCone = linspace(x1, x2, numConePoints);
        rCone = linspace(r1, r3, numConePoints);
        
        % Combine profile
        xConverging = [xThroatArc, xCone, xChamberArc];
        rConverging = [rThroatArc, rCone, rChamberArc];
        
        % Calculate volumes
        [xs, idx] = sort(xConverging);
        volumeConverging = pi * trapz(xs, rConverging(idx).^2);
        volumeChamberCylindrical = volumeChamber - volumeConverging;
        
        % Adjust parameters if invalid geometry
        if volumeChamberCylindrical > 0
            break;
        end

        if contractionRatio > contractionRatioFloor
            contractionRatio = max(contractionRatioFloor, contractionRatio - 0.5);
        else
            chamberFactor = chamberFactor * 0.8;
        end
        
    end
    
    if volumeChamberCylindrical <= 0

        warning('Cannot close chamber volume even at CR=%.1f, min blend (L*=%.2f) — raise L*\n', contractionRatio, Engine.Geometry.LStar);
        lengthChamberCylindrical = 0;

    else
        
        lengthChamberCylindrical = volumeChamberCylindrical / areaChamber;

    end
    
    Engine.Geometry.contractionRatio = contractionRatio;
    lengthChamber = min(xConverging) - lengthChamberCylindrical;
       
    % ------------------------------------------------------------------------
    
    % Diverging profile calculator - using Rao's TOP bell design
    
    % ------------------------------------------------------------------------

    % ------------------------------------------------------------------------
    % Rao data from https://github.com/ravi4ram/Bell-Nozzle/blob/master/bell_nozzle.py
    % ------------------------------------------------------------------------
    EXP_RATIO_CASES = [4, 5, 10, 20, 30, 40, 50, 100];
    THETA_N_60 = [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 38.0, 40.0];
    THETA_N_80 = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5];
    THETA_N_90 = [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0];
    THETA_E_60 = [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2];
    THETA_E_80 = [14.0, 13.0, 11.0, 9.0, 8.5, 8.0, 7.5, 7.0];
    THETA_E_90 = [11.5, 10.5, 8.0, 7.0, 6.5, 6.0, 6.0, 6.0];
    % ------------------------------------------------------------------------
    
    switch Engine.Geometry.CONE_PERCENT
        case 60, thetaN = THETA_N_60; thetaE = THETA_E_60;
        case 80, thetaN = THETA_N_80; thetaE = THETA_E_80;
        case 90, thetaN = THETA_N_90; thetaE = THETA_E_90;
        otherwise, error('CONE_PERCENT must be 60, 80, or 90');
    end
    
    thetaN = deg2rad(interp1(EXP_RATIO_CASES, thetaN, Engine.Specs.areaRatio, 'pchip'));
    thetaE = deg2rad(interp1(EXP_RATIO_CASES, thetaE, Engine.Specs.areaRatio, 'pchip'));
    Engine.Geometry.thetaE = thetaE;
    
    % Throat arc downstream
    rArcD = 0.382 * radiusThroat;
    angleDiverging = linspace(-pi/2, thetaN - pi/2, 20);
    xDivergingArc = rArcD * cos(angleDiverging);
    rDivergingArc = radiusThroat + rArcD + rArcD * sin(angleDiverging);
    
    % Inflection point between blend and parabola
    xInflection = xDivergingArc(end);
    rInflection = rDivergingArc(end);
    
    % Nozzle length
    rExit = sqrt(Engine.Specs.areaRatio) * radiusThroat;
    lengthCone15 = (rExit - radiusThroat) / tan(deg2rad(15));
    lengthNozzle = Engine.Geometry.CONE_PERCENT * lengthCone15 / 100;
    
    % ------------------------------------------------------------------------
    % Parabola
    % ------------------------------------------------------------------------
    
    xExit = lengthNozzle;
    mN = tan(thetaN);
    mE = tan(thetaE);
    
    % Control points
    xControlPoints = (rExit - rInflection + mN * xInflection - mE * xExit) / (mN - mE);
    rControlPoints = rInflection + mN*(xControlPoints - xInflection);
    
    % Create bezier curve
    t = linspace(0, 1, 80).';
    xParabolic = (1-t).^2 * xInflection + 2*(1-t).*t * xControlPoints + t.^2 * xExit;
    rParabolic = (1-t).^2 * rInflection + 2*(1-t).*t * rControlPoints + t.^2 * rExit;
    
    % Combine profile
    Engine.Geometry.xDiverging = [xDivergingArc, xParabolic.'];
    Engine.Geometry.rDiverging = [rDivergingArc, rParabolic.'];

    % ------------------------------------------------------------------------
    % Full
    % ------------------------------------------------------------------------
    
    numCylidricalPoints = 2;
    xCylindrical = linspace(lengthChamber, min(xConverging), numCylidricalPoints);
    rCylindrical = radiusChamber * ones(1, numCylidricalPoints);

    Engine.Geometry.xCoordinates = [xCylindrical, flip(xConverging), Engine.Geometry.xDiverging(:,2:end)];
    Engine.Geometry.rCoordinates = [rCylindrical, flip(rConverging), Engine.Geometry.rDiverging(:,2:end)];

    Engine.Chamber.diameter = diameterChamber;
    Engine.Throat.diameter = diameterThroat;
    Engine.Exit.diameter = diameterExit;
    Engine.Chamber.volume = volumeChamber;
    Engine.Chamber.area = areaChamber;
    Engine.Geometry.volumeConverging = volumeConverging;
    Engine.Geometry.lengthChamber = -lengthChamber;
    Engine.Geometry.lengthNozzle = max(Engine.Geometry.xDiverging);
    Engine.Geometry.lengthEngine = Engine.Geometry.lengthNozzle + Engine.Geometry.lengthChamber;
    
end