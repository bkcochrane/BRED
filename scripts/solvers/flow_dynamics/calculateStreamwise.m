function Engine = calculateStreamwise(Engine, GasData, region)

    switch region

        case 'diverging', X = Engine.X_DIVERGING(:); R = Engine.R_DIVERGING(:);
        case 'full', X = Engine.X_COORDINATES(:); R = Engine.R_COORDINATES(:);
        otherwise, error('region must be ''diverging'' or ''full''');

    end

    [~, GasData.B_I_CIRC] = chemistrySolve(Engine.OF_RATIO, GasData);
    H0 = Engine.H0;
    S0_R = Engine.S0_R;

    % Sort ascending and drop the duplicate throat point
    [X, iu] = unique(X);
    R = R(iu);
    N = numel(X);

    RHOs = zeros(N,1); % Densities
    Vs = RHOs; % Velocities
    MUs = RHOs; % Dynamic viscosities
    Ps = RHOs; % Pressures
    Ts = RHOs; % Temperatures
    Ms = RHOs; % Mach numbers
    ks = RHOs; % Thermal conductivities

    % Geometric areas are inflated by 1/(eta*Cd) — the isentropic march needs
    % the effective flow area, or the throat never reaches M = 1
    AREA_SCALE = Engine.ETA_COMBUSTION * Engine.DISCHARGE_COEFF;

    % Throat - reference point for both marches
    [~, iThroat] = min(abs(X));
    [Tt, N_Jt] = store(iThroat, Engine.P_THROAT, Engine.T_THROAT, Engine.N_J_THROAT, 0);

    % March from throat
    % Suddenly frozen flow principle at A/At = 1.3
    march(iThroat+1 : 1 : N, 'supersonic', Tt, N_Jt); % To exit
    march(iThroat-1 : -1 : 1, 'subsonic', Tt, N_Jt); % To chamber

    Engine.STREAM_X = X;
    Engine.STREAM_R = R;
    Engine.STREAM_RHO = RHOs;
    Engine.STREAM_V = Vs;
    Engine.STREAM_MU = MUs;
    Engine.STREAM_P = Ps;
    Engine.STREAM_T = Ts;
    Engine.STREAM_M = Ms;
    Engine.STREAM_k = ks;

    Engine = applySeparation(Engine);

    % ------------------------------------------------------------------------
    % March function
    % ------------------------------------------------------------------------
    function march(idx, branch, T, N_J)

        P_PREV = Engine.P_THROAT;
        N_J_FROZEN = N_J;
        
        for i = idx
            AREA_LOCAL = pi * R(i)^2 * AREA_SCALE;
            
            % Freeze chemistry when supersonic and A/At > 1.3
            IS_FROZEN = strcmp(branch, 'supersonic') && ((AREA_LOCAL / Engine.AREA_THROAT) > 1.3);
            
            if IS_FROZEN
                P = solvePressure(S0_R, H0, Engine, GasData, T, N_J_FROZEN, Engine.MDOT, AREA_LOCAL, P_PREV, branch, true);
                [T, N_J_FROZEN] = store(i, P, T, N_J_FROZEN, true);
            else
                P = solvePressure(S0_R, H0, Engine, GasData, T, N_J, Engine.MDOT, AREA_LOCAL, P_PREV, branch, false);
                [T, N_J] = store(i, P, T, N_J, false);
                N_J_FROZEN = N_J; % Keep syncing until the freeze point
            end
           
            P_PREV = P;
        end
    end


    
    % ------------------------------------------------------------------------
    % Solve and update, then seed value to next point
    % ------------------------------------------------------------------------
    function [T, N_J] = store(i, P, T_guess, N_J_guess, IS_FROZEN)

        [T, N_J, prop] = equilibriumSolve('SP', S0_R, P, GasData, T_guess, N_J_guess, 1, IS_FROZEN);

        RHOs(i) = prop.RHO;
        Vs(i) = sqrt(max(2*(H0 - prop.H), 0));
        MUs(i) = prop.MU;
        Ps(i) = P;
        Ts(i) = prop.T;
        [Ms(i), ~] = findMach(P, S0_R, H0, GasData, T, N_J);
        ks(i) = prop.k;

    end

end