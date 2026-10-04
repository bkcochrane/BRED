function Engine = solveEngine(Engine, GasData)

    F_TARGET = Engine.TARGET_THRUST; 
    NUM_SPECIES = GasData.NUM_SPECIES;
    
    MDOT = 1.0; % Unit mass flow
    [H_REACTANTS, GasData.B_I_CIRC] = chemistrySolve(Engine.OF_RATIO, GasData);

    P_CHAMBER = Engine.P_CHAMBER;
    
    [T_CHAMBER, N_J_CHAMBER, propellant_chamber] = solveChamber(P_CHAMBER, NUM_SPECIES, H_REACTANTS, GasData);
    H0 = propellant_chamber.H;
    S0_R = propellant_chamber.S_R;
    
    % Solve per unit mass flow (ideal physics)
    [P_THROAT, T_THROAT, N_J_THROAT, V_THROAT, AREA_THROAT_PER_MDOT, C_STAR_IDEAL, MACH_T] = solveThroat(P_CHAMBER, H0, S0_R, T_CHAMBER, N_J_CHAMBER, GasData, MDOT);
    Engine.P_THROAT = P_THROAT;

    if strcmp(Engine.EXIT_MODE, 'pressure')

        [P_EXIT, T_EXIT, N_J_EXIT, V_EXIT_IDEAL, AREA_EXIT_PER_MDOT] = solveExit(H0, S0_R, Engine.P_EXIT, GasData, T_THROAT, N_J_THROAT, MDOT);

    else

        AREA_EXIT_PER_MDOT = Engine.AREA_RATIO_TARGET * AREA_THROAT_PER_MDOT;
        P_EXIT = solvePressure(S0_R, H0, Engine, GasData, T_THROAT, N_J_THROAT, MDOT, AREA_EXIT_PER_MDOT, P_THROAT, 'supersonic', 1);
        [T_EXIT, N_J_EXIT, propellant_exit] = equilibriumSolve('SP', S0_R, P_EXIT, GasData, T_THROAT, N_J_THROAT, 0);
        V_EXIT_IDEAL = sqrt(max(2*(H0 - propellant_exit.H), 0));

    end
    
    % Apply efficiencies
    C_STAR = C_STAR_IDEAL * Engine.ETA_COMBUSTION;
    V_EXIT = V_EXIT_IDEAL * Engine.ETA_NOZZLE;
    
    % Size to target thrust
    F_PER_MDOT = V_EXIT + (P_EXIT - Engine.P_AMBIENT) * AREA_EXIT_PER_MDOT;

    if F_PER_MDOT <= 0
        error(['Ideal thrust is negative (F/mdot = %.0f N per kg/s): nozzle is ', ...
               'over-expanded at eps=%.0f, Pc=%.1f bar, Pa=%.0f Pa. ', ...
               'Raise Pc, or lower the area ratio, Pa or exit pressure ratio.'], ...
               F_PER_MDOT, Engine.AREA_RATIO, Engine.P_CHAMBER/1e5, Engine.P_AMBIENT);
    end

    MDOT = F_TARGET / F_PER_MDOT; % kg/s
    AREA_THROAT = AREA_THROAT_PER_MDOT * MDOT / Engine.ETA_COMBUSTION / Engine.DISCHARGE_COEFF; %m2
    AREA_EXIT = AREA_EXIT_PER_MDOT * MDOT / Engine.ETA_COMBUSTION / Engine.DISCHARGE_COEFF; %m2
    F = MDOT * F_PER_MDOT; % N
    ISP = F / (MDOT * 9.80665); % s
    
    Engine.MDOT = MDOT;
    Engine.MDOT_OX = MDOT * Engine.OF_RATIO / (Engine.OF_RATIO + 1);
    Engine.MDOT_FUEL = MDOT / (Engine.OF_RATIO + 1);
    Engine.THRUST = F;
    Engine.ISP = ISP;
    Engine.C_STAR = C_STAR;
    Engine.C_STAR_IDEAL = C_STAR_IDEAL;
    Engine.AREA_RATIO = AREA_EXIT / AREA_THROAT;
    
    Engine.T_CHAMBER = T_CHAMBER;
    Engine.CP_CHAMBER = propellant_chamber.CP;
    Engine.T_THROAT = T_THROAT;
    Engine.T_EXIT = T_EXIT;
    Engine.P_THROAT = P_THROAT;
    Engine.P_EXIT = P_EXIT;
    Engine.V_THROAT = V_THROAT;
    Engine.V_EXIT = V_EXIT;
    Engine.MACH_THROAT = MACH_T;
    
    Engine.AREA_THROAT = AREA_THROAT;
    Engine.AREA_EXIT = AREA_EXIT;
    Engine.DIAMETER_THROAT = 2*sqrt(AREA_THROAT/pi);
    Engine.DIAMETER_EXIT = 2*sqrt(AREA_EXIT/pi);
    
    Engine.MBAR = propellant_chamber.MBAR;
    Engine.RHO_CHAMBER = Engine.P_CHAMBER / ((GasData.R / propellant_chamber.MBAR) * T_CHAMBER);
    Engine.H0 = H0;
    Engine.S0_R = S0_R;
    Engine.N_J_CHAMBER = N_J_CHAMBER;
    Engine.N_J_THROAT = N_J_THROAT;
    Engine.N_J_EXIT = N_J_EXIT;

end
