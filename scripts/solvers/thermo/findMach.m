function [MACH, a] = findMach(P, S0_R, H0, GAS_DATA, T0, N_J0)

    [T, N_J, Propellant] = equilibriumSolve('SP', S0_R, P, GAS_DATA, T0, N_J0, 0);

    V  = sqrt(max(2*(H0 - Propellant.H), 0));
    
    A_IJ = GAS_DATA.A_IJ;
    NUM_ELEM = size(A_IJ, 1);
    
    % Separate gaseous components
    if isfield(GAS_DATA, 'CONDENSED_LIST')
        CANDIDATES = GAS_DATA.CONDENSED_LIST(:).';
    else
        CANDIDATES = [];
    end
    C = CANDIDATES(N_J(CANDIDATES) > 1e-10);
    G = setdiff(1:GAS_DATA.NUM_SPECIES, CANDIDATES);
    NC = numel(C);
    
    N_g = sum(N_J(G));
    
    % Build Jacobian matrix
    iN = NUM_ELEM + NC + 1;
    M = zeros(iN, iN);
    M(1:NUM_ELEM, 1:NUM_ELEM) = (A_IJ(:, G) .* N_J(G).') * A_IJ(:, G).';
    M(1:NUM_ELEM, 1:NUM_ELEM) = M(1:NUM_ELEM, 1:NUM_ELEM) + eye(NUM_ELEM) * 1e-15;

    for q = 1:NC
        M(1:NUM_ELEM, NUM_ELEM + q) = A_IJ(:, C(q));
        M(NUM_ELEM + q, 1:NUM_ELEM) = A_IJ(:, C(q)).';
    end

    M(1:NUM_ELEM, iN) = A_IJ(:, G) * N_J(G);
    M(iN, 1:NUM_ELEM) = M(1:NUM_ELEM, iN).';
    M(iN, iN) = 0;
    
    % --------------------------------------------------------------------
    % Derivative with respect to ln(P)
    B_P = zeros(iN, 1);
    B_P(1:NUM_ELEM) = A_IJ(:, G) * N_J(G);
    B_P(iN) = N_g;

    X_P = M \ B_P;
    dlnNg_dlnP = X_P(iN);
    dlnV_dlnP_T = dlnNg_dlnP - 1;
    % --------------------------------------------------------------------
    % Derivative with respect to ln(T) at constant P
    [CP_R, H_RT, ~, ~] = nasa9(GAS_DATA.A_NASA, GAS_DATA.B_NASA, GAS_DATA.T_HI, T);
    
    B_T = zeros(iN, 1);
    B_T(1:NUM_ELEM) = -A_IJ(:, G) * (N_J(G) .* H_RT(G));
    for q = 1:NC
        B_T(NUM_ELEM + q) = -H_RT(C(q));
    end
    B_T(iN) = -sum(N_J(G) .* H_RT(G));

    X_T = M \ B_T;
    pi_T = X_T(1:NUM_ELEM);
    dnC_dlnT_P = X_T(NUM_ELEM + 1 : NUM_ELEM + NC); % Condensed derivatives
    dlnNg_dlnT = X_T(iN);
    dlnV_dlnT_P = dlnNg_dlnT + 1;
    % --------------------------------------------------------------------

    dlnNj_dlnT = dlnNg_dlnT + A_IJ(:, G).' * pi_T + H_RT(G);

    R = GAS_DATA.R;

    % Equilibrium Cp
    Cp_FROZEN = sum(N_J .* CP_R) * R;
    Cp_G = sum(N_J(G) .* H_RT(G) .* dlnNj_dlnT) * R; 
    Cp_C = sum(H_RT(C) .* dnC_dlnT_P) * R;
    Cp_EQUILIBRIUM = Cp_FROZEN + Cp_G + Cp_C;
    
    % Equilibrium Cv
    Cv_EQUILIBRIUM = Cp_EQUILIBRIUM + (N_g * R) * (dlnV_dlnT_P^2) / dlnV_dlnP_T;
    
    % Isentropic Exponent (gamma)
    GAMMA_EQUILIBRIUM = Cp_EQUILIBRIUM / Cv_EQUILIBRIUM;
    
    % Analytical Speed of Sound
    % a^2 = - gamma / dlnV_dlnP_T * (P/rho) --> where (P/rho) = N_g * R * T
    a_SQUARED = -(GAMMA_EQUILIBRIUM / dlnV_dlnP_T) * N_g * R * T;
    a = sqrt(a_SQUARED);
    
    MACH = V / a;

end