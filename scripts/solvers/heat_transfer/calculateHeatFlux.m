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
    mach = Engine.Stream.mach;
    gammaGas = Engine.Stream.gamma;
    cpGas = Engine.Stream.cp;
    numX = numel(x);
    
    tempCoolant = Engine.Cooling.Coolant.T_SUPPLY * ones(numX, 1);
    prandtlGas = Engine.Stream.prandtl;
    rhoGas = Engine.Stream.rho;
    v = Engine.Stream.v;
    
    ds = hypot(diff(x), diff(r));
    rMean = 0.5 * (r(1:end-1) + r(2:end));
    areaSegment = 2 * pi .* rMean .* ds;
    
    % Size cooling channels based on geometry
    Engine = sizeCooling(Engine);
    widthStation = Engine.Cooling.widthChannelStation;
    heightStation = Engine.Cooling.heightChannelStation;
    
    % Pre-calculate adiabatic wall temperature (tempAW)
    mach2_vec = (gammaChamber - 1) / 2 .* mach.^2;
    tempAW = tempChamber .* (1 + prandtlChamber^0.33 .* mach2_vec) ./ (1 + mach2_vec);
    Engine.Wall.tempAW = tempAW;
    
    hG = zeros(numX, 1);
    hL = zeros(numX, 1);
    hEff = zeros(numX, 1);
    q = zeros(numX, 1);
    tempWG = 800 * ones(numX, 1);
    tempWI = 800 * ones(numX, 1);
    
    for k = 1:20
        tempCoolantOld = tempCoolant;
        
        if strcmp(method, 'Ievlev')
            Engine.Wall.tempWG = tempWG;
            Ievlev = calculateIevlev(Engine, PropData);
            zT = Ievlev.zT;
            beta2 = Ievlev.beta2;
            rSpecific = Ievlev.rSpecific;
            rSpecific1500K = Ievlev.rSpecific1500K;
        end

        if Engine.Cooling.FILM_COOLING_FRACTION > 0

            filmCoolingFraction = Engine.Cooling.FILM_COOLING_FRACTION;
            mdotFuel = Engine.Fuel.mdot;
            mdotFilmCooling = mdotFuel * filmCoolingFraction;
            cpFilm = PropData.Fuel.CP;
            muFilm = PropData.Fuel.MU;
            kFilm = PropData.Fuel.K;
            prandtlFilm = muFilm * cpFilm / kFilm;

        end
    
        for i = 1:numX
            areaLocal = pi * r(i)^2;
            areaRatioLocal = areaThroat / areaLocal;
            mach2 = mach2_vec(i);
            
            hCoolant = coolantH(Engine, widthStation(i), heightStation(i), tempCoolant(i));
            R_wall = Engine.Wall.THICKNESS / Engine.Wall.K;
            R_total = (1 / hCoolant) + R_wall;
            RELAX = 0.5;
    
            for iter = 1:50
                if strcmp(method, 'Bartz')
                    sigma = 1 / ((0.5 * (max(tempWG(i), 1) / tempChamber) * (1 + mach2) + 0.5)^0.68 * (1 + mach2)^0.12);
                    
                    hG(i) = sigma * (0.026 / diameterThroat^0.2) * (muChamber^0.2 * cpChamber / prandtlChamber^0.6) ...
                          * (P_CHAMBER / cStar)^0.8 * (diameterThroat / radiusCurvature)^0.1 * areaRatioLocal^0.9;
                    
                elseif strcmp(method, 'Ievlev')
                    tempE = rSpecific(i) * tempChamber / rSpecific1500K;
                    tempWGBar = min(max(tempWG(i) / tempE, 0.01), 0.95); 
                    
                    denomzTz = max(1 - tempWGBar + 0.1 * beta2(i), 1e-3);
                    numerzTz = 1 - beta2(i) + beta2(i) * (1 - 0.08696 * (1 - beta2(i)) / denomzTz);
                    zTz = (1.769 * max(numerzTz / denomzTz, 1e-3))^0.54;
                    
                    z = max(abs(zT(i) / zTz), 10);
                    
                    baseNum = max(1 - 0.21 * ((1 - prandtlGas(i)) / prandtlGas(i)^(4/3)) * beta2(i) / max(1 - tempWGBar, 1e-3), 1e-3);
                    numerHG = (1 / zTz)^(0.089 * prandtlGas(i)^-0.56) * baseNum^0.9225;
                    
                    denomHG = (307.8 + 54.8 * log10(prandtlGas(i) / 19.5)^2) * prandtlGas(i)^0.45 * z^0.08 - 650;
                    denomHGclamped = max(abs(denomHG), 25); 
                    
                    alphaT = numerHG / denomHGclamped;
                    hG(i) = alphaT * rhoGas(i) * v(i) * cpChamber;
                end

                if exist('mdotFilmCooling', 'var') && mdotFilmCooling > 0

                    % From Stechman et al. (1969) doi.org/10.2514/3.29545
                    
                    dx = max(x(i) - x(1), 1e-4) + 0.0015;

                    turbulentHeatTransfer = 0.0288 * cpFilm / (prandtlFilm^(1/3) * muFilm^0.2 * dx^0.2);

                    massFluxFilm = mdotFilmCooling / (2 * pi * r(i));
                    reynoldsFilm = massFluxFilm * dx / muFilm;
                    eta = 0.5 + 0.5 * exp(- (reynoldsFilm / 1800).^1.3);
                    
                    turbulentVelocityNumer = eta * mdotFilmCooling * v(i) * hG(i) * prandtlGas(i)^(2/3) * rhoGas(i);
                    turbulentVelocityDenom = pi * r(i) * cpGas(i);
                    
                    turbulentVelocity = (turbulentVelocityNumer / turbulentVelocityDenom)^0.4;

                    hL(i) = turbulentHeatTransfer * turbulentVelocity;

                    transitionMdot = 0.01 * mdotFuel;
                    alpha = 1 / (1 + exp(-50 * (mdotFilmCooling - transitionMdot)));
                    
                    hEff(i) = (1 - alpha) * hG(i) + alpha * hL(i);

                else

                    hL(i) = 0;

                    hEff(i) = hG(i);

                end

                tempWGNew = (tempCoolant(i) + hEff(i) * R_total * tempAW(i)) / (1 + hEff(i) * R_total);
                
                if abs(tempWGNew - tempWG(i)) < 1e-3
                    tempWG(i) = tempWGNew; 
                    break; 
                end
                
                tempWG(i) = (1 - RELAX) * tempWG(i) + RELAX * tempWGNew;
            end
            
            % Synchronize final station heat flux and coolant-side wall temperature
            q(i) = hEff(i) * (tempAW(i) - tempWG(i));
            tempWI(i) = tempCoolant(i) + q(i) * (1 / hCoolant);
        end
    
        % March coolant temperature from nozzle exit (numX) to injector (1)
        for i = numX-1:-1:1
            qSegment = 0.5 * (q(i) + q(i+1)) * areaSegment(i);
            tempCoolant(i) = tempCoolant(i+1) + qSegment / (Engine.Cooling.Coolant.mdot * Engine.Cooling.Coolant.CP);
        end
        
        if max(abs(tempCoolant - tempCoolantOld)) < 0.1
            break;
        end
    end
        
    % Final consistent pass with converged coolant temperatures
    for i = 1:numX
        hIdx = coolantH(Engine, widthStation(i), heightStation(i), tempCoolant(i));
        R_total = (1 / hIdx) + (Engine.Wall.THICKNESS / Engine.Wall.K);
        tempWG(i) = (tempCoolant(i) + hEff(i) * R_total * tempAW(i)) / (1 + hEff(i) * R_total);
        q(i) = hEff(i) * (tempAW(i) - tempWG(i));
        tempWI(i) = tempCoolant(i) + q(i) * (1 / hIdx);
    end
    
    Engine.Wall.q = q;
    Engine.Wall.hEff = hEff;
    Engine.Wall.tempAW = tempAW;
    Engine.Wall.tempWG = tempWG;
    Engine.Wall.tempWI = tempWI;
    Engine.Wall.tempCoolant = tempCoolant;
end