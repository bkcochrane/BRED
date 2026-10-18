function Engine = calculateHeatFlux(Engine, PropData)

    % --------------------------------------------------------------------

    % Using the Bartz model at chamber conditions

    % --------------------------------------------------------------------

    muChamber = Engine.Stream.mu(1);
    cpChamber = Engine.Chamber.cp;
    kChamber = Engine.Stream.k(1);
    prandtlChamber = muChamber * cpChamber / kChamber;
    r = PropData.R;
    gammaChamber = cpChamber / (cpChamber - r/Engine.Chamber.mbar);

    cStar = Engine.Specs.cStar;
    P_CHAMBER = Engine.Chamber.P;
    diameterThroat = Engine.Throat.diameter;
    radiusCurvature = 0.382 * diameterThroat/2;
    areaThroat = Engine.Throat.area;
    tempChamber = Engine.Chamber.temp;

    x = Engine.Stream.x;
    r = Engine.Stream.r;
    mach = Engine.Stream.mach;
    numX = numel(x);

    ds = hypot(diff(x), diff(r));
    rMean = 0.5 * (r(1:end-1) + r(2:end));
    areaSegment = 2*pi .* rMean .* ds;

    filmPercent = Engine.Cooling.FILM_COOLING_FRACTION;
    fuelMdot = Engine.Fuel.mdot;

    if filmPercent > 0

        Engine.Cooling.Film.mdot = fuelMdot * filmPercent;

    end
    
    if Engine.Cooling.REGEN_COOLING == true

        Engine.Cooling.Coolant = PropData.Fuel;
        Engine.Cooling.Coolant.mdot = fuelMdot * (1 - filmPercent);
        tempCoolant = PropData.Fuel.T_SUPPLY * ones(numX,1);

    else

        Engine.Cooling.Coolant.P_REF = [1, 2, 5, 10, 20, 50, 100, 150, 220.6]; 
        Engine.Cooling.Coolant.T_REF = [373.1, 393.4, 425.0, 453.2, 485.6, 537.1, 584.3, 615.4, 647.1]; 
        
        Engine.Cooling.Coolant.tSat = @(P) interp1(P_REF, T_REF, P/1e5, 'pchip', 'extrap');
        
        tempSupply = 298; % K
        Engine.Cooling.Coolant.T_SUPPLY = tempSupply;
        tempCoolant = tempSupply * ones(numX,1);
        
        Engine.Cooling.Coolant.RHO = 1000.0; % kg/m3
        Engine.Cooling.Coolant.MU = 8.9e-4; % Pa*s
        Engine.Cooling.Coolant.CP = 4184.0; % J/(kg*K)
        Engine.Cooling.Coolant.K = 0.606; % W/(m*K)

    end

    % Size cooling channels based on geometry
    Engine = sizeCooling(Engine);
    widthStation = Engine.Cooling.widthChannelStation;
    heightStation = Engine.Cooling.heightChannelStation;

    hG = zeros(numX,1);
    tempAW = zeros(numX,1);
    q = zeros(numX,1);
    tempWG = 800 * ones(numX,1);
    tempWI = 800 * ones(numX,1);

    for k = 1:20

        tempCoolantOld = tempCoolant;
    
        for i = 1:numX
    
            areaLocal = pi * r(i)^2;
            areaRatioLocal = areaThroat / areaLocal;
            mach2 = (gammaChamber-1)/2 * mach(i)^2;
            tempAW(i) = tempChamber * (1 + prandtlChamber^0.33*mach2) / (1 + mach2);
            
            hCoolant = coolantH(Engine, widthStation(i), heightStation(i), tempCoolant(i));

            RELAX = 0.5;
    
            for iter = 1:50
    
                sigma = 1 / ( (0.5 * (max(tempWG(i), 1) / tempChamber) * (1 + mach2) + 0.5)^0.68 * (1 + mach2)^0.12);
                
                hG(i) = sigma * (0.026 / diameterThroat^0.2) * (muChamber^0.2 * cpChamber / prandtlChamber^0.6) ...
                             * (P_CHAMBER / cStar)^0.8 * (diameterThroat / radiusCurvature)^0.1 * areaRatioLocal^0.9;
                
                q(i) = hG(i) * (tempAW(i) - tempWG(i));
                
                tempWGNew = tempCoolant(i) + q(i) * (1 / hCoolant + Engine.Wall.THICKNESS / Engine.Wall.K);
                
                if abs(tempWGNew - tempWG(i)) < 1e-3, tempWG(i) = tempWGNew; break; end
                
                tempWGNew = max(tempWGNew, Engine.Cooling.Coolant.T_SUPPLY);
    
                tempWG(i) = (1-RELAX) * tempWG(i) + RELAX * tempWGNew;
            
            end

            tempWI(i) = tempCoolant(i) + q(i) * (1 / hCoolant);
    
        end
    
        for i = numX-1:-1:1
            
            qSegment = 0.5*(q(i) + q(i+1)) * areaSegment(i);
            tempCoolant(i) = tempCoolant(i+1) + qSegment / (Engine.Cooling.Coolant.mdot * Engine.Cooling.Coolant.CP);
        
        end

        if max(abs(tempCoolant - tempCoolantOld)) < 0.1
            break;
        end

    end
        
    for i = 1:numX
        hIdx = coolantH(Engine, widthStation(i), heightStation(i), tempCoolant(i));
        tempWG(i) = tempCoolant(i) + q(i)*(1/hIdx + Engine.Wall.THICKNESS/Engine.Wall.K);
        tempWI(i) = tempCoolant(i) + q(i)*(1/hIdx);
    end

    Engine.Wall.q = q;
    Engine.Wall.hG = hG;
    Engine.Wall.tempAW = tempAW;
    Engine.Wall.tempWG = tempWG;
    Engine.Wall.tempWI = tempWI;
    Engine.Wall.tempCoolant = tempCoolant;

end