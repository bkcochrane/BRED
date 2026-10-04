function Engine = calculateHeatFlux(Engine, GAS_DATA)

    % --------------------------------------------------------------------

    % Using the Bartz model at chamber conditions

    % --------------------------------------------------------------------

    MU_CHAMBER = Engine.STREAM_MU(1);
    CP_CHAMBER = Engine.CP_CHAMBER;
    k_CHAMEBER = Engine.STREAM_k(1);
    PR_C = MU_CHAMBER * CP_CHAMBER / k_CHAMEBER;
    R = GAS_DATA.R;
    GAMMA_CHAMBER = CP_CHAMBER / (CP_CHAMBER - R/Engine.MBAR);

    C_STAR = Engine.C_STAR;
    PRESSURE_CHAMBER = Engine.P_CHAMBER;
    DIAMETER_THROAT = Engine.DIAMETER_THROAT;
    RAD_CURV = 0.382 * DIAMETER_THROAT/2;
    AREA_THROAT = Engine.AREA_THROAT;
    T_CHAMBER = Engine.T_CHAMBER;

    X = Engine.STREAM_X;
    R = Engine.STREAM_R;
    M = Engine.STREAM_M;
    N = numel(X);

    ds = hypot(diff(X), diff(R));
    R_MEAN = 0.5 * (R(1:end-1) + R(2:end));
    AREA_SEG = 2*pi .* R_MEAN .* ds;

    H_G = zeros(N,1);
    T_AW = zeros(N,1);
    Q = zeros(N,1);
    T_WG = 800*ones(N,1);
    T_WI = 800*ones(N,1);

    % Size cooling channels based on geometry
    Engine = sizeCooling(Engine);
    WIDTH_STATION = Engine.CHANNEL_WIDTH_STATION;
    HEIGHT_STATION = Engine.CHANNEL_HEIGHT_STATION;
    
    for i = 1:N

        AREA_LOCAL = pi * R(i)^2;
        AREA_AREA = AREA_THROAT / AREA_LOCAL;
        MACH2 = (GAMMA_CHAMBER-1)/2 * M(i)^2;
        T_AW(i) = T_CHAMBER * (1 + PR_C^0.33*MACH2) / (1 + MACH2);
        
        H_COOLANT = coolantH(Engine, WIDTH_STATION(i), HEIGHT_STATION(i), Engine.T_COOLANT_IN);

        T_WG(i) = 800;

        for iter = 1:50

            SIGMA = 1 / ( (0.5 * (max(T_WG(i), 1) / T_CHAMBER) * (1 + MACH2) + 0.5)^0.68 * (1 + MACH2)^0.12);
            
            H_G(i) = SIGMA * (0.026 / DIAMETER_THROAT^0.2) * (MU_CHAMBER^0.2 * CP_CHAMBER / PR_C^0.6) ...
                         * (PRESSURE_CHAMBER / C_STAR)^0.8 * (DIAMETER_THROAT / RAD_CURV)^0.1 * AREA_AREA^0.9;
            
            Q(i) = H_G(i) * (T_AW(i) - T_WG(i));
            
            T_WG_NEW = Engine.T_COOLANT_IN + Q(i)*(1 / H_COOLANT ...
                        + Engine.WALL_THICKNESS / Engine.WALL_K);
            
            if abs(T_WG_NEW - T_WG(i)) < 1e-3, T_WG(i) = T_WG_NEW; break; end
            
            T_WG_NEW = max(T_WG_NEW, Engine.T_COOLANT_IN);

            RELAX = 0.5;
            T_WG(i) = (1-RELAX) * T_WG(i) + RELAX * T_WG_NEW;
        
        end

    end

    % Coolant enters at exit, flows to nozzle
    T_COOL = zeros(N,1);
    T_COOL(N) = Engine.T_COOLANT_IN;

    for i = N-1:-1:1
        
        Q_SEG = 0.5*(Q(i) + Q(i+1)) * AREA_SEG(i);
        T_COOL(i) = T_COOL(i+1) + Q_SEG / (Engine.MDOT_COOLANT * Engine.CP_COOLANT);
    
    end
        
    for i = 1:N
        h_i = coolantH(Engine, WIDTH_STATION(i), HEIGHT_STATION(i), T_COOL(i));
        T_WG(i) = T_COOL(i) + Q(i)*(1/h_i + Engine.WALL_THICKNESS/Engine.WALL_K);
        T_WI(i) = T_COOL(i) + Q(i)*(1/h_i);
    end

    Engine.WALL_Q = Q;
    Engine.WALL_HG = H_G;
    Engine.WALL_TAW = T_AW;
    Engine.WALL_TWG = T_WG;
    Engine.WALL_TWI = T_WI;
    Engine.WALL_TCOOLANT = T_COOL;

end