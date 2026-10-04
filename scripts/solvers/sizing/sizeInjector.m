function Engine = sizeInjector(Engine)

    RHO_FUEL = Engine.RHO_FUEL; % kg/m3
    RHO_OX = Engine.RHO_OX; % kg/m3
    Cd = Engine.INJECTOR_DISCHARGE_COEFFICIENT;
    t = Engine.INJECTOR_PLATE_THICKNESS; % m
    MDOT_FUEL = Engine.MDOT_FUEL; % kg/s
    MDOT_OX = Engine.MDOT_OX; % kg/s
    DIAMETER_MIN = Engine.INJECTOR_PORT_MIN_DIAMETER; % m
    GAMMA = Engine.INJECTOR_IMPINGEMENT_ANGLE; % radians

    P_CHAMBER = Engine.P_CHAMBER; % Pa

    dP_FUEL = Engine.P_FUEL - P_CHAMBER; % Pa
    dP_OX = Engine.P_OX - P_CHAMBER; % Pa

    AREA_FUEL = MDOT_FUEL / (Cd * sqrt(2 * RHO_FUEL * dP_FUEL)); % m2
    AREA_OX = MDOT_OX / (Cd * sqrt(2 * RHO_OX * dP_OX)); % m2

    AREA_RATIO = 0.5 * AREA_OX / AREA_FUEL;
    DIAMETER_RATIO = sqrt(AREA_RATIO);

    if DIAMETER_RATIO > 1
        DIAMETER_FUEL = DIAMETER_MIN;
        DIAMETER_OX = DIAMETER_MIN * DIAMETER_RATIO;
    else
        DIAMETER_FUEL = DIAMETER_MIN / DIAMETER_RATIO;
        DIAMETER_OX = DIAMETER_MIN;
    end

    N_PORTS_FUEL = ceil(AREA_FUEL / (pi / 4 * DIAMETER_FUEL^2));
    N_PORTS_OX = 2 * N_PORTS_FUEL;

    JET_DIAMETERS = 5; % diameters
    L_JET = JET_DIAMETERS * DIAMETER_OX; % m - along the ox jet
    L_IMPINGE = L_JET * cos(GAMMA); % m - axial from the face
    SPACING = L_IMPINGE * tan(GAMMA);
    
    WEB = SPACING - 0.5*(DIAMETER_FUEL + DIAMETER_OX);
    if WEB < 0.001
        warning('Land between fuel and ox ports is %.2f mm - too thin to drill', WEB*1000);
    end

    V_FUEL = MDOT_FUEL / (RHO_FUEL * AREA_FUEL); % m/s
    V_OX = MDOT_OX / (RHO_OX * AREA_OX); % m/s

    MOMENTUM_FUEL = MDOT_FUEL * V_FUEL;
    MOMENTUM_OX = MDOT_OX * V_OX * cos(GAMMA);
    RATIO = MOMENTUM_FUEL / MOMENTUM_OX;

    STIFFNESS_FUEL = dP_FUEL / P_CHAMBER;
    STIFFNESS_OX = dP_OX / P_CHAMBER;

    LD_FUEL = t / DIAMETER_FUEL;
    LD_OX = t / DIAMETER_OX / cos(GAMMA);

    if RATIO > 1
        warning('Fuel momentum greater than oxidiser momentum - mixing efficiency will drop.')
    end
    if RATIO < 0.67
        warning('Oxidiser momentum far greater than fuel momentum - mixing efficiency will drop.')
    end
    if STIFFNESS_FUEL < 0.15 || STIFFNESS_OX < 0.15
        warning('Injector stiffness below 0.15 - combustion instability risk.');
    end
    if LD_FUEL < 3 || LD_OX < 3
        warning('Orifice L/D below 3 - discharge coefficient will be unstable.');
    end
    
    Engine.INJECTOR_AREA_FUEL = AREA_FUEL;
    Engine.INJECTOR_AREA_OX = AREA_OX;
    Engine.INJECTOR_DIAMETER_FUEL = DIAMETER_FUEL;
    Engine.INJECTOR_DIAMETER_OX = DIAMETER_OX;
    Engine.INJECTOR_N_PORTS_FUEL = N_PORTS_FUEL;
    Engine.INJECTOR_N_PORTS_OX = N_PORTS_OX;
    Engine.INJECTOR_N_PORT_SETS = N_PORTS_FUEL;
    Engine.INJECTOR_V_FUEL = V_FUEL;
    Engine.INJECTOR_V_OX = V_OX;
    Engine.INJECTOR_MOMENTUM_RATIO = RATIO;
    Engine.INJECTOR_STIFFNESS_FUEL = STIFFNESS_FUEL;
    Engine.INJECTOR_STIFFNESS_OX = STIFFNESS_OX;
    
end
