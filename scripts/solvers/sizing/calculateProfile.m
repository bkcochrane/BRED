function Engine = calculateProfile(Engine)

    VOLUME_CHAMBER = Engine.L_STAR * Engine.AREA_THROAT; % m3
    AREA_CHAMBER = Engine.AREA_THROAT * Engine.CONTRACTION_RATIO; % m2
    DIAMETER_CHAMBER = 2*sqrt(AREA_CHAMBER/pi); % m
    DIAMETER_THROAT = 2*sqrt(Engine.AREA_THROAT/pi); % m
    DIAMETER_EXIT = 2*sqrt(Engine.AREA_EXIT/pi); % m

    % ------------------------------------------------------------------------
    % Calculate converging section profile 
    % using two blends and a converging angle beta
    % ------------------------------------------------------------------------
    
    BETA = Engine.BETA;
    RADIUS_THROAT = DIAMETER_THROAT / 2;
    RADIUS_CHAMBER = DIAMETER_CHAMBER / 2;
    
    THROAT_FACTOR = Engine.BLEND_RADIUS_THROAT_FACTOR;
    CHAMBER_FACTOR = Engine.BLEND_RADIUS_CHAMBER_FACTOR;

    RADIAL_SPAN = RADIUS_CHAMBER - RADIUS_THROAT;
    BLEND_RADIUS_THROAT = THROAT_FACTOR * RADIUS_THROAT;
    if BLEND_RADIUS_THROAT * (1-cos(BETA)) > 0.45 * RADIAL_SPAN
        BLEND_RADIUS_THROAT = 0.45 * RADIAL_SPAN / (1-cos(BETA));
    end
    
    % Blend to throat
    THETA = linspace(0, BETA, 50);
    X_T_ARC = -BLEND_RADIUS_THROAT * sin(THETA);
    R_T_ARC = RADIUS_THROAT + BLEND_RADIUS_THROAT * (1 - cos(THETA));
    X_1 = X_T_ARC(end);
    R_1 = R_T_ARC(end);

    CR = Engine.CONTRACTION_RATIO;
    CR_FLOOR = 2;
    
    for outer = 1:20
        AREA_CHAMBER = Engine.AREA_THROAT * CR;
        DIAMETER_CHAMBER = 2*sqrt(AREA_CHAMBER/pi);
        RADIUS_CHAMBER = DIAMETER_CHAMBER / 2;

        MAX_CHAMBER_DROP = 0.5 * RADIAL_SPAN;
        BLEND_RADIUS_CHAMBER = CHAMBER_FACTOR * RADIUS_CHAMBER;
        BLEND_DROP = BLEND_RADIUS_CHAMBER * (1 - cos(BETA));
        
        if BLEND_DROP > MAX_CHAMBER_DROP
            BLEND_RADIUS_CHAMBER = MAX_CHAMBER_DROP / (1 - cos(BETA));
        end

        R_3 = RADIUS_CHAMBER - BLEND_RADIUS_CHAMBER * (1 - cos(BETA));
        if R_3 <= R_1
            BLEND_RADIUS_CHAMBER = 0.8 * (RADIUS_CHAMBER - R_1) / (1 - cos(BETA));
            R_3 = RADIUS_CHAMBER - BLEND_RADIUS_CHAMBER * (1 - cos(BETA));
        end
        X_2 = X_1 - (R_3 - R_1) / tan(BETA);
        
        % Blend to chamber
        PHI = linspace(BETA, 0, 50);
        X_C_CENTRE = X_2 - BLEND_RADIUS_CHAMBER * sin(BETA);
        X_C_ARC = X_C_CENTRE + BLEND_RADIUS_CHAMBER * sin(PHI);
        R_C_ARC = RADIUS_CHAMBER - BLEND_RADIUS_CHAMBER * (1 - cos(PHI));

        N_CONE = 30;
        X_CONE = linspace(X_1, X_2, N_CONE);
        R_CONE = linspace(R_1, R_3, N_CONE);
        
        % Combine profile
        X_CONVERGING = [X_T_ARC, X_CONE, X_C_ARC];
        R_CONVERGING = [R_T_ARC, R_CONE, R_C_ARC];
        
        % Calculate volumes
        [Xs, idx] = sort(X_CONVERGING);
        VOLUME_CONVERGING = pi * trapz(Xs, R_CONVERGING(idx).^2);
        VOLUME_CHAMBER_CYLINDRICAL = VOLUME_CHAMBER - VOLUME_CONVERGING;
        
        % Adjust parameters if invalid geometry
        if VOLUME_CHAMBER_CYLINDRICAL > 0
            break;
        end

        if CR > CR_FLOOR
            CR = max(CR_FLOOR, CR - 0.5);
        else
            CHAMBER_FACTOR = CHAMBER_FACTOR * 0.8;
        end
    end
    
    if VOLUME_CHAMBER_CYLINDRICAL <= 0

        warning('Cannot close chamber volume even at CR=%.1f, min blend (L*=%.2f) — raise L*', CR, Engine.L_STAR);
        LENGTH_CHAMBER_CYLINDRICAL = 0;

    else
        
        LENGTH_CHAMBER_CYLINDRICAL = VOLUME_CHAMBER_CYLINDRICAL / AREA_CHAMBER;

    end
    
    Engine.CONTRACTION_RATIO_USED = CR;
    LENGTH_CHAMBER = min(X_CONVERGING) - LENGTH_CHAMBER_CYLINDRICAL;
        
    % ------------------------------------------------------------------------
    
    % Diverging profile calculator - using Rao's TOP bell design
    
    % ------------------------------------------------------------------------

    % ------------------------------------------------------------------------
    % Rao data from https://github.com/ravi4ram/Bell-Nozzle/blob/master/bell_nozzle.py
    % ------------------------------------------------------------------------
    EXP_RATIO_CASES = [4, 5, 10, 20, 30, 40, 50, 100];
    THETA_N_60 = [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 38.0, 40.0];
    THETA_N_80 = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5];
    THETA_N_90 = [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0];
    THETA_E_60 = [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2];
    THETA_E_80 = [14.0, 13.0, 11.0, 9.0, 8.5, 8.0, 7.5, 7.0];
    THETA_E_90 = [11.5, 10.5, 8.0, 7.0, 6.5, 6.0, 6.0, 6.0];
    % ------------------------------------------------------------------------
    
    switch Engine.CONE_PERCENT
        case 60, THETA_N_SWITCH = THETA_N_60; THETA_E_SWITCH = THETA_E_60;
        case 80, THETA_N_SWITCH = THETA_N_80; THETA_E_SWITCH = THETA_E_80;
        case 90, THETA_N_SWITCH = THETA_N_90; THETA_E_SWITCH = THETA_E_90;
        otherwise, error('CONE_PERCENT must be 60, 80, or 90');
    end
    
    THETA_N = deg2rad(interp1(EXP_RATIO_CASES, THETA_N_SWITCH, Engine.AREA_RATIO, 'pchip'));
    THETA_E = deg2rad(interp1(EXP_RATIO_CASES, THETA_E_SWITCH, Engine.AREA_RATIO, 'pchip'));
    Engine.THETA_E = THETA_E;
    
    % Throat arc downstream
    R_ARC_D = 0.382 * RADIUS_THROAT;
    ANGLE = linspace(-pi/2, THETA_N - pi/2, 50);
    X_DIV_ARC = R_ARC_D * cos(ANGLE);
    R_DIV_ARC = RADIUS_THROAT + R_ARC_D + R_ARC_D * sin(ANGLE);
    
    % Inflection point between blend and parabola
    X_N = X_DIV_ARC(end);
    R_N = R_DIV_ARC(end);
    
    % Nozzle length
    R_EXIT = sqrt(Engine.AREA_RATIO) * RADIUS_THROAT;
    L_CONE15 = (R_EXIT - RADIUS_THROAT) / tan(deg2rad(15));
    L_NOZZLE = Engine.CONE_PERCENT * L_CONE15 / 100;
    
    % ------------------------------------------------------------------------
    % Parabola
    % ------------------------------------------------------------------------
    
    X_EXIT = L_NOZZLE;
    m_N = tan(THETA_N);
    m_E = tan(THETA_E);
    
    % Control points
    X_CP = (R_EXIT - R_N + m_N * X_N - m_E * X_EXIT) / (m_N - m_E);
    R_CP = R_N + m_N*(X_CP - X_N);
    
    % Create bezier curve
    t = linspace(0, 1, 80).';
    X_PARA = (1-t).^2 * X_N + 2*(1-t).*t * X_CP + t.^2 * X_EXIT;
    R_PARA = (1-t).^2 * R_N + 2*(1-t).*t * R_CP + t.^2 * R_EXIT;
    
    % Combine profile
    Engine.X_DIVERGING = [X_DIV_ARC, X_PARA.'];
    Engine.R_DIVERGING = [R_DIV_ARC, R_PARA.'];

    % ------------------------------------------------------------------------
    % Full
    % ------------------------------------------------------------------------
    
    N_CYL = 20;
    X_CYL = linspace(LENGTH_CHAMBER, min(X_CONVERGING), N_CYL);
    R_CYL = RADIUS_CHAMBER * ones(1, N_CYL);

    Engine.X_COORDINATES = [X_CYL, flip(X_CONVERGING), Engine.X_DIVERGING(:,2:end)];
    Engine.R_COORDINATES = [R_CYL, flip(R_CONVERGING), Engine.R_DIVERGING(:,2:end)];

    Engine.DIAMETER_CHAMBER = DIAMETER_CHAMBER;
    Engine.DIAMETER_THROAT = DIAMETER_THROAT;
    Engine.DIAMETER_EXIT = DIAMETER_EXIT;
    Engine.VOLUME_CHAMBER = VOLUME_CHAMBER;
    Engine.LENGTH_CHAMBER = -LENGTH_CHAMBER;
    Engine.NOZZLE_LENGTH = max(Engine.X_DIVERGING);

end