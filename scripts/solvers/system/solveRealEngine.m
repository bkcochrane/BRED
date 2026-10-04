function Engine = solveRealEngine(Ideal, GasData)

    % ------------------------------------------------------------------------
    % Combustion efficiency
    % ------------------------------------------------------------------------
    
    % Calculate how long the propellant remains in the chamber 
    t_RES = Ideal.RHO_CHAMBER * Ideal.VOLUME_CHAMBER / Ideal.MDOT; % s
    
    % ------------------------------------------------------------------------
    % Estimate time to vapourise the propellant
    % From: https://www.researchgate.net/publication/336456135_Vaporization-controlled_simplified_model_for_liquid_propellant_rocket_engine_combustion_chamber_design
    
    CP_FILM = Ideal.CP_CHAMBER;
    RHO_RP1 = 807; % kg/m3
    QS = 246e3; % J/kg Heat of vapourisation
    T_BOIL_RP1 = 490; % K - RP-1 value >> LOX of 153K
    T_FILM = 0.5 * (Ideal.T_CHAMBER + T_BOIL_RP1); % K
    k_FILM = 0.002 + 0.0001 * T_FILM; % W/(m K)
    
    B = CP_FILM * (T_FILM - T_BOIL_RP1) / QS; % Spalding number
    K = 8 * k_FILM / (RHO_RP1 * CP_FILM) * log(1 + B);
    t_VAP = Ideal.D0^2 / K; % s 
    % ------------------------------------------------------------------------
    
    if t_RES >= t_VAP
        FRAC_VAPOUR = 1.0;
    else
        FRAC_VAPOUR = 1 - (1 - t_RES / t_VAP)^1.5;
    end
    
    Engine = Ideal;

    RADIUS_THROAT = Engine.DIAMETER_THROAT / 2;
    RADIUS_EXIT = Engine.DIAMETER_EXIT / 2;
    
    % ------------------------------------------------------------------------
    
    for nozzle_it = 1:3
    
        % ------------------------------------------------------------------------
        % Nozzle efficiency
        % ------------------------------------------------------------------------
        
        Engine = calculateProfile(Engine);
        
        % ------------------------------------------------------------------------
        % Divergence losses
        
        RADIUS_BAR_EXIT = RADIUS_EXIT / RADIUS_THROAT;
        L_BAR = Engine.NOZZLE_LENGTH / RADIUS_THROAT;
        RADIUS_BAR_CHAMBER = 1 + (RADIUS_BAR_EXIT - 1) / L_BAR;
        n_div = 1.45 * (RADIUS_BAR_CHAMBER^0.25) - 0.005 * RADIUS_BAR_CHAMBER;

        [~, ~, propellant] = equilibriumSolve('SP', Engine.S0_R, Engine.P_EXIT, GasData, Engine.T_THROAT, Engine.N_J_THROAT, 0);
        R_SPECIFIC_EXIT = GasData.R / propellant.MBAR;
        R_SPECIFIC_CHAMBER = GasData.R / Engine.MBAR;
        k_AVE = log(Engine.P_CHAMBER / Engine.P_EXIT) / log((Engine.P_CHAMBER / Engine.P_EXIT) * (R_SPECIFIC_EXIT * propellant.T) / (R_SPECIFIC_CHAMBER * Engine.T_CHAMBER));
        
        B = 1.52 * (exp(-30 * (k_AVE - 1)) + 0.1);
        ZETA_D = 1 - B * (exp(n_div * (1 - RADIUS_BAR_EXIT / RADIUS_BAR_CHAMBER)) - 1) / (exp(n_div) - 1);
    
        % ------------------------------------------------------------------------
        % Frictional losses
        [~, GasData.B_I_CIRC] = chemistrySolve(Engine.OF_RATIO, GasData);
        
        Engine = calculateStreamwise(Engine, GasData, 'diverging');
        
        Re = Engine.STREAM_RHO .* Engine.STREAM_V .* max(Engine.STREAM_X,1e-6) ./ Engine.STREAM_MU; % Reynolds number
        Cf_SKIN = 0.026 ./ Re.^(1/7); % Turbulent skin friction (1/7th Power Law)
        TAU_W  = 0.5 * Cf_SKIN .* Engine.STREAM_RHO .* Engine.STREAM_V.^2; % Pa - Wall shear stress
        
        dX_AXIAL = diff(Engine.STREAM_X);
        dR = diff(Engine.STREAM_R);
        ds_surf = hypot(dX_AXIAL, dR);
        cos_THETA = dX_AXIAL ./ ds_surf;
        
        R_AVE   = 0.5 * (Engine.STREAM_R(1:end-1) + Engine.STREAM_R(2:end));
        TAU_AVE = 0.5 * (TAU_W(1:end-1) + TAU_W(2:end));
        
        integral_val = sum(TAU_AVE .* R_AVE .* cos_THETA .* dX_AXIAL);
        Cf_IDEAL = Engine.THRUST / (Engine.P_CHAMBER * Engine.AREA_THROAT); 
        
        ZETA_F = 1 - (2 * pi * integral_val) / (Cf_IDEAL * Engine.AREA_THROAT * Engine.P_CHAMBER);
        
        % ------------------------------------------------------------------------
        % Multi-phase flow

        Z_CONDENSED = propellant.N_CONDENSED * 12.011;
        
        if Z_CONDENSED > 1e-5
            C_s = 2000; % Approx specific heat of solid soot J/(kg K)
            w_e = Engine.STREAM_V(end);
            
            ZETA_ZW = 1 - (Z_CONDENSED * C_s / w_e^2) * (Engine.T_CHAMBER - propellant.T * (1 + log(Engine.T_CHAMBER / prop_e.T)));
            ZETA_zT = 1 - (Z_CONDENSED / 2);
            ZETA_Z = 0.2 * ZETA_ZW + 0.8 * ZETA_zT;
        else
            ZETA_Z = 1.0;
        end
        
        Engine.ETA_NOZZLE = ZETA_D * ZETA_F * ZETA_Z;

        % ------------------------------------------------------------------------
        % Discharge coefficient
        Engine.DISCHARGE_COEFF = 0.97;
        
        % ------------------------------------------------------------------------
        % Combustion efficiency
        % ------------------------------------------------------------------------
        % Finite rate kinematics loss
        h0_kin = 6.7e-4; 
        a_kin = 0.35; 
        b_kin = 0.8;
        
        XI_r_prime = (h0_kin / RADIUS_THROAT)^a_kin * (Engine.P_AMBIENT / Engine.P_CHAMBER)^b_kin * log10(RADIUS_EXIT / RADIUS_THROAT);
        XI_r_double_prime = max(0, (2.1 - log(Engine.P_CHAMBER / 2e6)) / 100);
        
        ZETA_R = (1 - XI_r_prime) * (1 - XI_r_double_prime);

        Engine.ETA_COMBUSTION = FRAC_VAPOUR * Engine.ETA_MIXING * ZETA_R;
        % ------------------------------------------------------------------------
        
        Engine = solveEngine(Engine, GasData);
    end

end