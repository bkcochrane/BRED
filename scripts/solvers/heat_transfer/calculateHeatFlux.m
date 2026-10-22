function Engine = calculateHeatFlux(Engine, PropData, method)

    muChamber = Engine.Stream.mu(1);
    cpChamber = Engine.Chamber.cp;
    kChamber = Engine.Stream.k(1);
    prandtlChamber = muChamber * cpChamber / kChamber;
    rGas = PropData.R;
    gammaChamber = cpChamber / (cpChamber - rGas / Engine.Chamber.mbar);
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
    
    tempCoolant = Engine.Cooling.Coolant.T_SUPPLY * ones(numX, 1);
    prandtlGas = Engine.Stream.prandtl;
    rhoGas = Engine.Stream.rho;
    v = Engine.Stream.v;

    filmCoolingFraction = Engine.Cooling.FILM_COOLING_FRACTION;

    if filmCoolingFraction > 0

        EngineNoFilm = Engine;
        EngineNoFilm.Cooling.FILM_COOLING_FRACTION = 0;
        ReferenceEngine = calculateHeatFlux(EngineNoFilm, PropData, method);

        Film = calculateFilmCooling(ReferenceEngine, PropData, filmCoolingFraction);
        qRef = ReferenceEngine.Wall.q;
        SRatio = Film.SRatio;
        Engine.Cooling.Film.STable = Film.STable;

    end
        
    ds = hypot(diff(x), diff(r));
    rMean = 0.5 * (r(1:end-1) + r(2:end));
    areaSegment = 2 * pi .* rMean .* ds;
    
    % Size cooling channels based on geometry
    Engine = sizeCooling(Engine);
    widthStation = Engine.Cooling.widthChannelStation;
    heightStation = Engine.Cooling.heightChannelStation;
    
    % Pre-calculate adiabatic wall temperature (tempAW)
    mach2Vec = (gammaChamber - 1) / 2 .* mach.^2;
    tempAW = tempChamber .* (1 + prandtlChamber^0.33 .* mach2Vec) ./ (1 + mach2Vec);
    Engine.Wall.tempAW = tempAW;
    tempAWClean = tempAW;
    
    hL = zeros(numX, 1);
    hG = zeros(numX, 1);
    q = zeros(numX, 1);
    if ~isfield(Engine.Wall, 'tempWG') || numel(Engine.Wall.tempWG) ~= numX
        Engine.Wall.tempWG = 800 * ones(numX, 1);
    end
    if ~isfield(Engine.Wall, 'tempWI') || numel(Engine.Wall.tempWI) ~= numX
        Engine.Wall.tempWI = 800 * ones(numX, 1);
    end

    [~, PropData.bCirc] = chemistrySolve(Engine.Specs.ofRatio, PropData);
    [~, nj1500, Propellant1500K] = equilibriumSolve('TP', 1500, Engine.Chamber.P, PropData, 1500, Engine.Chamber.nj, 0);
    Engine.Chamber.rSpecific1500K = Propellant1500K.rSpecific;
    Engine.Chamber.nj1500 = nj1500;
    
    for k = 1:20
        tempCoolantOld = tempCoolant;

        if strcmp(method, 'Bartz')
            
            TempWall = calculateBartz(Engine);

        end
        
        if strcmp(method, 'Ievlev')

            TempWall = calculateIevlev(Engine, PropData);
            if k == 1
                qIevlev = TempWall.q;
            else
                qIevlev = 0.5 * qIevlev + 0.5 * TempWall.q;
            end

        end

        hG = TempWall.hG;
        tempWG = Engine.Wall.tempWG;
        tempWI = Engine.Wall.tempWI;
    
        for i = 1:numX

            hCoolant = coolantH(Engine, widthStation(i), heightStation(i), tempCoolant(i));
            R_wall = Engine.Wall.THICKNESS / Engine.Wall.K;
            R_total = (1 / hCoolant) + R_wall;
            RELAX = 0.5;

            for iter = 1:50

                if filmCoolingFraction > 0
                    tempWGNew = tempCoolant(i) + qRef(i) * SRatio(i) * R_total;
                elseif strcmp(method, 'Ievlev')
                    tempWGNew = tempCoolant(i) + qIevlev(i) * R_total;
                else
                    tempWGNew = (tempCoolant(i) + hG(i) * R_total * tempAW(i)) / (1 + hG(i) * R_total);
                end

                if abs(tempWGNew - tempWG(i)) < 1e-3
                    tempWG(i) = tempWGNew;
                    break;
                end

                tempWG(i) = (1 - RELAX) * tempWG(i) + RELAX * tempWGNew;
            end

            if filmCoolingFraction > 0
                q(i) = qRef(i) * SRatio(i);
            elseif strcmp(method, 'Ievlev')
                q(i) = qIevlev(i);
            else
                q(i) = hG(i) * (tempAW(i) - tempWG(i));
            end

            tempWI(i) = tempCoolant(i) + q(i) * (1 / hCoolant);

        end
    
        % March coolant temperature from nozzle exit (numX) to injector (1)
        for i = numX-1:-1:1
            qSegment = 0.5 * (q(i) + q(i+1)) * areaSegment(i);
            tempCoolant(i) = max(tempCoolant(i+1) + qSegment / (Engine.Cooling.Coolant.mdot * Engine.Cooling.Coolant.CP), Engine.Cooling.Coolant.T_SUPPLY);
        end

        Engine.Wall.hG = hG;
        Engine.Wall.tempAW = tempAW;
        Engine.Wall.tempWG = tempWG;
        Engine.Wall.tempWI = tempWI;
        
        if max(abs(tempCoolant - tempCoolantOld)) < 0.1
            break;
        end
        
    end

    hCoolant = zeros(numX, 1);
    % Final consistent pass with converged coolant temperatures
    for i = 1:numX
        hIdx = coolantH(Engine, widthStation(i), heightStation(i), tempCoolant(i));
        R_total = (1 / hIdx) + (Engine.Wall.THICKNESS / Engine.Wall.K);

        if filmCoolingFraction > 0
            q(i) = qRef(i) * SRatio(i);
            tempWG(i) = tempCoolant(i) + q(i) * R_total;
        elseif strcmp(method, 'Ievlev')
            tempWG(i) = tempCoolant(i) + qIevlev(i) * R_total;
            q(i) = qIevlev(i);
        else
            tempWG(i) = (tempCoolant(i) + hG(i) * R_total * tempAW(i)) / (1 + hG(i) * R_total);
            q(i) = hG(i) * (tempAW(i) - tempWG(i));
        end

        tempWI(i) = tempCoolant(i) + q(i) * (1 / hIdx);
        hCoolant(i) = hIdx;
    end

    Engine.Wall.q = q;
    Engine.Wall.hG = hG;
    Engine.Wall.tempAW = tempAW;
    Engine.Wall.tempWG = tempWG;
    Engine.Wall.tempWI = tempWI;
    Engine.Wall.tempCoolant = tempCoolant;
    Engine.Wall.hCoolant = hCoolant;
        
end