function Engine = sizeCoolingFlow(Engine, GAS_DATA)

    T_SAT = waterTSat(Engine.P_COOLANT);
    MARGIN = 40; % K below saturation, design margin
    T_LIMIT = T_SAT - MARGIN;

    mdot = Engine.MDOT_COOLANT; % kg/s - Starting guess

    for it = 1:30
        Engine.MDOT_COOLANT = mdot;
        Engine = calculateHeatFlux(Engine, GAS_DATA);
        T_OUT = max(Engine.WALL_TCOOLANT);

        if T_OUT <= T_LIMIT, break; end

        % Delta T roughly 1/mDot
        rise_now = T_OUT - Engine.T_COOLANT_IN;
        rise_max = T_LIMIT - Engine.T_COOLANT_IN;
        mdot = mdot * (rise_now / rise_max) * 1.05;
    end

    Engine.COOLANT_TSAT = T_SAT;
    Engine.COOLANT_MARGIN = T_SAT - max(Engine.WALL_TCOOLANT);
    
end