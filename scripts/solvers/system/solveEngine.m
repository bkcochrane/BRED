function Engine = solveEngine(Engine, PropData)

    THRUST_TARGET = Engine.Specs.TARGET_THRUST; 
    NUM_SPECIES = PropData.NUM_SPECIES;
    
    UNIT_MDOT = 1.0; % Unit mass flow
    [hReactants, bCirc] = chemistrySolve(Engine.Specs.ofRatio, PropData);
    
    PropData.bCirc = bCirc;

    P_CHAMBER = Engine.Chamber.P;
    
    [tempChamber, njChamber, PropellantChamber] = solveChamber(P_CHAMBER, NUM_SPECIES, hReactants, PropData);
    h0 = PropellantChamber.h;
    s0R = PropellantChamber.sR;
    
    % Solve per unit mass flow (ideal physics)
    [pThroat, tempThroat, njThroat, vThroat, areaThroatPerMdot, cStarIdeal, machThroat] = solveThroat(P_CHAMBER, h0, s0R, tempChamber, njChamber, PropData, UNIT_MDOT);
    Engine.Throat.p = pThroat;

    if strcmp(Engine.Specs.EXIT_MODE, 'pressure')

        [pExit, tempExit, njExit, vExitIdeal, areaExitPerMdot] = solveExit(h0, s0R, Engine.Exit.p, PropData, tempThroat, njThroat, UNIT_MDOT);

    else

        areaExitPerMdot = Engine.Specs.AREA_RATIO_TARGET * areaThroatPerMdot;
        pExit = solvePressure(s0R, h0, Engine, PropData, tempThroat, njThroat, UNIT_MDOT, areaExitPerMdot, pThroat, 'supersonic', 1);
        [tempExit, njExit, propellant_exit] = equilibriumSolve('SP', s0R, pExit, PropData, tempThroat, njThroat, 0);
        vExitIdeal = sqrt(max(2*(h0 - propellant_exit.h), 0));

    end
    
    % Apply efficiencies
    cStar = cStarIdeal * Engine.Efficiency.combustion;
    vExit = vExitIdeal * Engine.Efficiency.nozzle;
    
    % Size to target thrust
    thrustPerMdot = vExit + (pExit - Engine.Specs.P_AMBIENT) * areaExitPerMdot;

    if thrustPerMdot <= 0
        error(['Ideal thrust is negative (F/mdot = %.0f N per kg/s): nozzle is ', ...
               'over-expanded at eps=%.0f, Pc=%.1f bar, Pa=%.0f Pa. ', ...
               'Raise Pc, or lower the area ratio, Pa or exit pressure ratio.'], ...
               thrustPerMdot, Engine.Specs.areaRatio, Engine.Chamber.P/1e5, Engine.Specs.P_AMBIENT);
    end
    
    if ~isfield(Engine.Geometry, 'LOCK') || Engine.Geometry.LOCK == false

        mdot = THRUST_TARGET / thrustPerMdot; % kg/s
        areaThroat = areaThroatPerMdot * mdot / Engine.Efficiency.combustion / Engine.Efficiency.cdThroat; % m 2
        areaExit = areaExitPerMdot * mdot / Engine.Efficiency.combustion / Engine.Efficiency.cdThroat; % m 2
        Engine.Throat.area = areaThroat;
        Engine.Exit.area = areaExit;
        Engine.Throat.diameter = 2*sqrt(areaThroat/pi);
        Engine.Exit.diameter = 2*sqrt(areaExit/pi);
        Engine.Specs.areaRatio = areaExit / areaThroat;

    else

        mdot = (Engine.Throat.area * Engine.Efficiency.combustion * Engine.Efficiency.cdThroat) / areaThroatPerMdot;

    end

    thrust = mdot * thrustPerMdot; % N
    Isp = thrust / (mdot * 9.80665); % s
    
    Engine.Specs.mdot = mdot;
    Engine.Ox.mdot = mdot * Engine.Specs.ofRatio / (Engine.Specs.ofRatio + 1);
    Engine.Fuel.mdot = mdot / (Engine.Specs.ofRatio + 1);
    Engine.Specs.thrust = thrust;
    Engine.Specs.Isp = Isp;
    Engine.Specs.cStar = cStar;
    Engine.Specs.cStarIdeal = cStarIdeal;
    
    Engine.Chamber.temp = tempChamber;
    Engine.Chamber.cp = PropellantChamber.cp;
    Engine.Throat.temp = tempThroat;
    Engine.Exit.temp = tempExit;
    Engine.Throat.p = pThroat;
    Engine.Exit.p = pExit;
    Engine.Throat.v = vThroat;
    Engine.Exit.v = vExit;
    Engine.Throat.mach = machThroat;
    
    Engine.Chamber.mbar = PropellantChamber.mbar;
    Engine.Chamber.rho = Engine.Chamber.P / ((PropData.R / PropellantChamber.mbar) * tempChamber);
    Engine.Chamber.h0 = h0;
    Engine.Chamber.s0R = s0R;
    Engine.Chamber.mu = PropellantChamber.mu;
    Engine.Chamber.k = PropellantChamber.k;
    Engine.Chamber.nj = njChamber;
    Engine.Throat.nj = njThroat;
    Engine.Exit.nj = njExit;

end
