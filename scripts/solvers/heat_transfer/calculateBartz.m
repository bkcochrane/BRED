function Bartz = calculateBartz(Engine)

    muChamber = Engine.Stream.mu(1);
    cpChamber = Engine.Chamber.cp;
    kChamber = Engine.Stream.k(1);
    prandtlChamber = muChamber * cpChamber / kChamber;
    cStar = Engine.Specs.cStar;
    P_CHAMBER = Engine.Chamber.P;
    diameterThroat = Engine.Throat.diameter;
    radiusCurvature = 0.382 * diameterThroat / 2;
    areaThroat = Engine.Throat.area;
    tempChamber = Engine.Chamber.temp;
    
    x = Engine.Stream.x;
    r = Engine.Stream.r;
    tempCore = Engine.Stream.temp;
    mach = Engine.Stream.mach;
    gammaGas = Engine.Stream.gamma;
    cpGas = Engine.Stream.cp;
    muGas = Engine.Stream.mu;
    numX = numel(x);
    
    % tempCoolant = Engine.Cooling.Coolant.T_SUPPLY * ones(numX, 1);
    % prandtlGas = Engine.Stream.prandtl;
    % rhoGas = Engine.Stream.rho;
    % v = Engine.Stream.v;
    % 
    % filmCoolingFraction = Engine.Cooling.FILM_COOLING_FRACTION;
    % 
    % ds = hypot(diff(x), diff(r));
    % rMean = 0.5 * (r(1:end-1) + r(2:end));
    % areaSegment = 2 * pi .* rMean .* ds;
    % 
    % % Size cooling channels based on geometry
    % Engine = sizeCooling(Engine);
    % widthStation = Engine.Cooling.widthChannelStation;
    % heightStation = Engine.Cooling.heightChannelStation;
    % 
    % % Pre-calculate adiabatic wall temperature (tempAW)
    % mach2Vec = (gammaChamber - 1) / 2 .* mach.^2;
    % tempAW = tempChamber .* (1 + prandtlChamber^0.33 .* mach2Vec) ./ (1 + mach2Vec);
    % Engine.Wall.tempAW = tempAW;
    
    hG = zeros(numX, 1);
    tempWG = Engine.Wall.tempWG;

    for i = 1:numX

        mach2 = (gammaGas(i) - 1) / 2 * mach(i)^2;
        areaRatioLocal = areaThroat / (pi * r(i)^2);
    
        sigma = 1 / ((0.5 * (max(tempWG(i), 1) / tempChamber) * (1 + mach2) + 0.5)^0.68 * (1 + mach2)^0.12);
                        
        hG(i) = sigma * (0.026 / diameterThroat^0.2) * (muChamber^0.2 * cpChamber / prandtlChamber^0.6) ...
              * (P_CHAMBER / cStar)^0.8 * (diameterThroat / radiusCurvature)^0.1 * areaRatioLocal^0.9;

    end

    Bartz.hG = hG;

end