function Engine = engine_RAPTOR

    Engine.Fuel.NAME = 'CH4';
    Engine.Ox.NAME = 'LOX';
    
    Engine.Specs.TARGET_THRUST = 2500e3; % N
    
    Engine.Specs.P_AMBIENT = 101325; % Pa - Design ambient pressure
    
    Engine.Specs.EXIT_MODE = 'area'; % Set to 'area' or 'pressure'
    if strcmpi(Engine.Specs.EXIT_MODE, 'pressure')
        Engine.Specs.EXIT_PR = 0.5;
        Engine.Exit.p = Engine.Specs.EXIT_PR * Engine.Specs.P_AMBIENT; % Pa
    else
        Engine.Specs.AREA_RATIO_TARGET = 34.34;
    end
    
    Engine.Specs.OF_RATIO_TARGET = 3.6;
    
    Engine.Chamber.P = 330e5; % Pa
    
    Engine.Geometry.LSTAR_TARGET = 5; % m
    Engine.Geometry.CONTRACTION_RATIO_TARGET = 6; % AREA_CHAMBER / AREA_THROAT
    
    Engine.Geometry.CONE_PERCENT = 90;
    Engine.Geometry.BETA_ANGLE = deg2rad(30); % radians - Converging angle
    Engine.Geometry.BLEND_RADIUS_THROAT_FACTOR = 1.5;
    Engine.Geometry.BLEND_RADIUS_CHAMBER_FACTOR = 2;
    
    Engine.Efficiency.MIXING = 0.94;
    Engine.Efficiency.KINETIC = 0.96;
    
    Engine.Cooling.REGEN_COOLING = true;
    Engine.Cooling.FILM_COOLING_FRACTION = 0.15;
    
    Engine.Cooling.CHANNEL_WIDTH_MIN = 0.00150641; % m
    Engine.Cooling.CHANNEL_HEIGHT_MIN = 0.02; % m
    Engine.Cooling.CHANNEL_HEIGHT_MAX = 0.0398964; % m
    Engine.Cooling.CHANNEL_ROUGHNESS = 5e-6; % m
    Engine.Cooling.LAND_WIDTH = 0.001; % m
    Engine.Cooling.N_STARTS = 0;
    Engine.Cooling.N_TURNS = 1;
    
    Engine.Wall.THICKNESS = 0.001; % m
    Engine.Wall.K = 230; % W/(m K)
    
    Engine.Fuel.P = 400e5; % Pa
    Engine.Ox.P = 400e5; % Pa
    Engine.Injector.Cd = 0.8;
    Engine.Injector.PLATE_THICKNESS = 0.010; % m
    Engine.Injector.PORT_MIN_DIAMETER = 0.0015; % m

end