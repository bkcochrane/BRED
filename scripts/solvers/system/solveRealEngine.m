function Engine = solveRealEngine(Ideal, PropData)

    % ------------------------------------------------------------------------
    % Combustion efficiency
    % ------------------------------------------------------------------------

    % Estimate time to vapourise the propellant - Spalding Ranz-Marshall
    
    vGas = Ideal.Specs.mdot / (Ideal.Chamber.rho * Ideal.Chamber.area);
    vJet = Ideal.Fuel.v;
    vDroplet = 0.5 * (vGas + vJet);
    vRel = abs(vGas - vDroplet);
    
    Droplets = {PropData.Fuel, PropData.Ox};
    d0s = [Ideal.Fuel.d0, Ideal.Ox.d0];
    timesVapourise = [0, 0];

    prandtlNumber = (Ideal.Chamber.cp * Ideal.Chamber.mu) / Ideal.Chamber.k; % Prandtl number

    for i = 1:2
        Prop = Droplets{i}; 
        d0 = d0s(i);
        
        cpFilm = Ideal.Chamber.cp;
        tempFilm = 0.5 * (Ideal.Chamber.temp + Prop.TEMP_BOIL); % K
        kFilm = Prop.calcK(tempFilm);
        
        reynoldsDroplet = (Ideal.Chamber.rho * vRel * d0) / Ideal.Chamber.mu;
        
        % Ranz-Marshall Nusselt Number
        nuStar = 2 + 0.552 * sqrt(reynoldsDroplet) * prandtlNumber^(1/3);
        
        B = cpFilm * (tempFilm - Prop.TEMP_BOIL) / Prop.QS; % Spalding transfer number
        
        K = (4 * nuStar) * kFilm / (Prop.RHO * cpFilm) * log(1 + B);
        
        timesVapourise(i) = d0^2 / K; % s

    end
    
    timeVapourise = max(timesVapourise);

    areaThroat = Ideal.Throat.area;
    volumeConverging = Ideal.Geometry.volumeConverging;
    areaChamber = Ideal.Chamber.area;
    
    LStarRequired = (timeVapourise * vDroplet * areaChamber + volumeConverging) / areaThroat;
    LSTAR_TARGET = Ideal.Geometry.LSTAR_TARGET;
    maxLStar = 1.35 * LSTAR_TARGET;
    
    LStar = max(LSTAR_TARGET, min(LStarRequired, maxLStar));
    
    volumeCylinder = LStar * areaThroat - volumeConverging;
    timeResidence = max(volumeCylinder, 0) / areaChamber / vDroplet; % converging section ignored, conservative
    fracVapour = 1 - (1 - min(timeResidence / timeVapourise, 1))^1.5;

    Ideal.Geometry.LStar = LStar;
    
    Engine = Ideal;

    radiusThroat = Engine.Throat.diameter / 2;
    radiusExit = Engine.Exit.diameter / 2;
    
    % ------------------------------------------------------------------------
    
    for nozzle_it = 1:3
    
        % ------------------------------------------------------------------------
        % Nozzle efficiency
        % ------------------------------------------------------------------------
        
        if ~isfield(Engine.Geometry, 'lock') || Engine.Geometry.lock == false
            Engine = calculateProfile(Engine);
        end
        
        % ------------------------------------------------------------------------
        % Divergence losses
        
        radiusBarExit = radiusExit / radiusThroat;
        LBar = Engine.Geometry.lengthNozzle / radiusThroat;
        radiusBarChamber = 1 + (radiusBarExit - 1) / LBar;
        nDiv = 1.45 * (radiusBarChamber^0.25) - 0.005 * radiusBarChamber;
        
        [~, ~, Propellant] = equilibriumSolve('SP', Engine.Chamber.s0R, Engine.Exit.p, PropData, Engine.Throat.temp, Engine.Throat.nj, 0);
        rSpecificExit = PropData.R / Propellant.mbar;
        rSpecificChamber = PropData.R / Engine.Chamber.mbar;
        kAverage = log(Engine.Chamber.P / Engine.Exit.p) / log((Engine.Chamber.P / Engine.Exit.p) * (rSpecificExit * Propellant.temp) / (rSpecificChamber * Engine.Chamber.temp));
        
        B = 1.52 * (exp(-30 * (kAverage - 1)) + 0.1);
        radiusRatio = max(radiusBarChamber, radiusBarExit) / min(radiusBarExit, radiusBarChamber);
        zetaD = 1 - B * (exp(nDiv * (1 - radiusRatio)) - 1) / (exp(nDiv) - 1);
        if zetaD <= 0 || zetaD > 1.0
            zetaD = 0.985;
        end
        
        % ------------------------------------------------------------------------
        % Frictional losses
        [~, PropData.bCirc] = chemistrySolve(Engine.Specs.ofRatio, PropData);
        
        Engine = calculateStreamwise(Engine, PropData, 'diverging');
        
        reynolds = Engine.Stream.rho .* Engine.Stream.v .* max(Engine.Stream.x,1e-6) ./ Engine.Stream.mu; % Reynolds number
        cfSkin = 0.026 ./ reynolds.^(1/7); % Turbulent skin friction (1/7th Power Law)
        tauW  = 0.5 * cfSkin .* Engine.Stream.rho .* Engine.Stream.v.^2; % Pa - Wall shear stress
        
        dXAxial = diff(Engine.Stream.x);
        dR = diff(Engine.Stream.r);
        dsSurf = hypot(dXAxial, dR);
        cosTheta = dXAxial ./ dsSurf;
        
        rAverage = 0.5 * (Engine.Stream.r(1:end-1) + Engine.Stream.r(2:end));
        tauAverage = 0.5 * (tauW(1:end-1) + tauW(2:end));
        
        integralVal = sum(tauAverage .* rAverage .* cosTheta .* dXAxial);
        cfIdeal = Engine.Specs.thrust / (Engine.Chamber.P * Engine.Throat.area); 
        
        zetaF = 1 - (2 * pi * integralVal) / (cfIdeal * Engine.Throat.area * Engine.Chamber.P);
        
        % ------------------------------------------------------------------------
        % Multi-phase flow

        zCondensed = Propellant.numCondensed * 12.011;
        
        if zCondensed > 1e-5
            cS = 2000; % Approx specific heat of solid soot J/(kg K)
            wE = Engine.Stream.v(end);
            
            zetaZW = 1 - (zCondensed * cS / wE^2) * (Engine.Chamber.temp - Propellant.temp * (1 + log(Engine.Chamber.temp / Propellant.temp)));
            zetaZT = 1 - (zCondensed / 2);
            zetaZ = 0.2 * zetaZW + 0.8 * zetaZT;
        else
            zetaZ = 1.0;
        end
        
        Engine.Efficiency.nozzle = zetaD * zetaF * zetaZ;

        % ------------------------------------------------------------------------
        % Discharge coefficient
        Engine.Efficiency.cdThroat = 0.97;
        
        % ------------------------------------------------------------------------
        % Combustion efficiency
        % ------------------------------------------------------------------------
        % Finite rate kinematics loss
        h0Kin = 6.7e-4; 
        aKin = 0.35; 
        bKin = 0.8;
        
        xiRPrime = (h0Kin / radiusThroat)^aKin * (Engine.Specs.P_AMBIENT / Engine.Chamber.P)^bKin * log10(radiusExit / radiusThroat);
        xiRPrimePrime = max(0, (2.1 - log(Engine.Chamber.P / 2e6)) / 100);
        
        zetaR = (1 - xiRPrime) * (1 - xiRPrimePrime);

        Engine.Efficiency.combustion = fracVapour * Engine.Efficiency.MIXING * zetaR;
      
        % ------------------------------------------------------------------------
        
        Engine = solveEngine(Engine, PropData);
    end

    if LStarRequired > maxLStar
        fprintf('[Warning] Required L* = %.2f m is greater than the maximum L* %.1f m. Only %.1f%% of propellant is vapourised: combustion efficiency has dropped to %.3f.\n', ...
            LStarRequired, maxLStar, fracVapour * 100, Engine.Efficiency.combustion);
    end  

end