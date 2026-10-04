function [T, N_J, Propellant, PI_I] = findEquilibrium(MODE, TARGET, P, GasData, T, N_J, NEED_MU, IS_FROZEN)

    if nargin < 8; IS_FROZEN = false; end

    % HP mode solves for fixed enthalpy (chamber pressure)
    % HP targets H_REACTANTS
    %
    % SP mode solves for insentropic expansion
    % SP targets S0_R

    A_IJ = GasData.A_IJ;
    B_I_CIRC = GasData.B_I_CIRC;
    A_NASA = GasData.A_NASA;
    B_NASA = GasData.B_NASA;
    T_HI = GasData.T_HI;
    MOLAR_MASSES = GasData.M;
    R = GasData.R;
    V_COEFFS = GasData.V_COEFFS;
    V_RANGES_MAX = GasData.V_RANGES_MAX;
    C_COEFFS = GasData.C_COEFFS;
    C_RANGES_MAX = GasData.C_RANGES_MAX;
    LOG_P = log(P / 1e5);
 
    NUM_SPECIES = numel(N_J);
    NUM_ELEM = size(A_IJ,1);
 
    % Partition species into gas (G) and active condensed (C)
    if isfield(GasData,'CONDENSED_LIST')
        COND_ALL = GasData.CONDENSED_LIST(:).';
    else
        COND_ALL = [];
    end

    if isfield(GasData,'ACTIVE_CONDENSED')
        C = GasData.ACTIVE_CONDENSED(:).';
    else
        C = [];
    end

    G  = setdiff(1:NUM_SPECIES, COND_ALL);
    NC = numel(C);
    N_J(G) = max(N_J(G), 1e-250);
 
    MAX_ITERATIONS = 100;
 
    for iter = 1:MAX_ITERATIONS
        
        % Skip matrix if equilibrium is frozen
        if IS_FROZEN && strcmp(MODE, 'SP')

            [CP_R, ~, S_R, ~] = nasa9(A_NASA, B_NASA, T_HI, T);
            N_TOTAL = sum(N_J(G));
            
            S_MU_R = S_R;
            S_MU_R(G) = S_R(G) - log(max(N_J(G), 1e-250) / N_TOTAL) - LOG_P;
            S_MIX_R = sum(N_J .* S_MU_R);
            
            % 1D Newton step for Temperature
            delta_ln_T = (TARGET - S_MIX_R) / sum(N_J .* CP_R);
            T = max(T * exp(delta_ln_T), 150);
            
            if abs(delta_ln_T) < 1e-5
                PI_I = zeros(NUM_ELEM, 1);
                break;
            end

            continue;
        end
 
        N_TOTAL = sum(N_J(G)); % Gas - moles only
        [CP_R, H_RT, ~, G_RT] = nasa9(A_NASA, B_NASA, T_HI, T);
 
        % Chemical potentials: gases carry the ln terms, condensed do not
        MU_RT = G_RT; % Condensed: mu = G_RT
        MU_RT(G) = G_RT(G) + log(N_J(G) / N_TOTAL) + LOG_P; % Gases add ln(x)+ln P
 
        NH  = N_J .* H_RT;
        NMU = N_J .* MU_RT;
 
        sz = NUM_ELEM + NC + 2; % Elements + condensed + dlnN + dlnT
        A  = zeros(sz, sz);
        B  = zeros(sz, 1);
 
        iN = NUM_ELEM + NC + 1; % row/col index of dln N
        iT = NUM_ELEM + NC + 2; % row/col index of dln T
 
        % Element rows (gas contribution + condensed columns)
        A(1:NUM_ELEM, 1:NUM_ELEM) = (A_IJ(:,G) .* N_J(G).') * A_IJ(:,G).';
        A(1:NUM_ELEM, 1:NUM_ELEM) = A(1:NUM_ELEM, 1:NUM_ELEM) + eye(NUM_ELEM) * 1e-15;
        for q = 1:NC
            A(1:NUM_ELEM, NUM_ELEM+q) = A_IJ(:, C(q)); % Condensed column = atom vector
        end
        A(1:NUM_ELEM, iN) = A_IJ(:,G) * N_J(G);
        A(1:NUM_ELEM, iT) = A_IJ(:,G) * NH(G);
        B(1:NUM_ELEM) = B_I_CIRC - A_IJ*N_J + A_IJ(:,G) * NMU(G);
 
        % Condensed rows: phase equilibrium sum a_ij pi_i + (H_c/RT)dlnT = mu_c/RT
        for q = 1:NC
            j = C(q);
            A(NUM_ELEM+q, 1:NUM_ELEM) = A_IJ(:,j).';
            A(NUM_ELEM+q, iT) = H_RT(j);
            B(NUM_ELEM+q) = MU_RT(j);
        end
 
        % Total-moles row (gases only)
        A(iN, 1:NUM_ELEM) = (A_IJ(:,G) * N_J(G)).';
        A(iN, iN) = sum(N_J(G)) - N_TOTAL;
        A(iN, iT) = sum(NH(G));
        B(iN) = N_TOTAL - sum(N_J(G)) + sum(NMU(G));
 
        % Energy / entropy row 
        if strcmp(MODE,'HP')
            H_MIX_RT = sum(NH); 
            A(iT, 1:NUM_ELEM) = A(1:NUM_ELEM, iT).';
            for q = 1:NC, A(iT, NUM_ELEM+q) = H_RT(C(q)); end
            A(iT, iN) = sum(NH(G));
            
            % FIXED: Substitution artifacts must only sum over Gas (G)
            A(iT, iT) = sum(N_J .* CP_R) + sum(NH(G) .* H_RT(G)); 
            B(iT) = TARGET/(R*T) - H_MIX_RT + sum(NMU(G) .* H_RT(G)); 
        else
            [~, ~, S_R_full, ~] = nasa9(A_NASA, B_NASA, T_HI, T);
            S_MU_R = S_R_full - MU_RT + G_RT; 
            S_MIX_R = sum(N_J .* S_MU_R);
            A(iT, 1:NUM_ELEM) = (A_IJ(:,G) * (N_J(G) .* S_MU_R(G))).';
            for q = 1:NC, A(iT, NUM_ELEM+q) = S_MU_R(C(q)); end
            
            % FIXED: Substitution artifacts must only sum over Gas (G)
            A(iT, iN) = sum(N_J(G) .* S_MU_R(G)); 
            A(iT, iT) = sum(N_J .* CP_R) + sum(N_J(G) .* H_RT(G) .* S_MU_R(G)); 
            B(iT) = TARGET - S_MIX_R + (N_TOTAL - sum(N_J(G))) + sum(N_J(G) .* S_MU_R(G) .* MU_RT(G)); 
        end
 
        % Solve
        if rcond(A) < 1e-15
            % Use Moore-Penrose pseudo-inverse for singular/badly scaled matrices
            X = pinv(full(A)) * B; 
        else
            X = A \ B;
        end
        PI_I = X(1:NUM_ELEM);
        DELTA_N_C = X(NUM_ELEM + (1:NC)); % Linear update for condensed
        DELTA_LN_N = X(iN);
        DELTA_LN_T = X(iT);
 
        % Gas log-updates derived from the multipliers
        DELTA_LN_NJ = zeros(NUM_SPECIES,1);
        DELTA_LN_NJ(G) = A_IJ(:,G).'*PI_I + DELTA_LN_N + H_RT(G)*DELTA_LN_T - MU_RT(G);
 
        % Convergence
        conv_gas = max(abs(N_J(G).*DELTA_LN_NJ(G)))/N_TOTAL;
        conv_con = 0; if NC>0, conv_con = max(abs(DELTA_N_C))/max(N_TOTAL,eps); end
        if conv_gas < 0.5e-5 && conv_con < 0.5e-5 && ...
           abs(DELTA_LN_N) < 0.5e-5 && abs(DELTA_LN_T) < 1e-4
            break;
        end
 
        % Step limiting
        big = (N_J(G)/N_TOTAL) > 1e-8;
        Gbig = G(big);
        l1 = 2 / max([5*abs(DELTA_LN_T), 5*abs(DELTA_LN_N), ...
                      max(abs(DELTA_LN_NJ(Gbig)))]);
        lambda = min([1, l1]);
 
        % Apply
        delta_T_safe = max(-10, min(10, lambda*DELTA_LN_T));
        gas_update_safe = max(-50, min(50, lambda*DELTA_LN_NJ(G)));
        
        T_new = T * exp(delta_T_safe);
        T = max(T_new, 150); 
        
        N_J(G) = max(N_J(G) .* exp(gas_update_safe), 1e-250);
        
        if NC > 0
            N_J(C) = max(N_J(C) + lambda*DELTA_N_C, 0); % Additive, can't be negative
        end
    end
 
    % Derived properties
    [~, H_RT, S_R, ~] = nasa9(A_NASA, B_NASA, T_HI, T);
    N_TOTAL = sum(N_J(G));
    S_MU_R = S_R;
    S_MU_R(G) = S_R(G) - log(N_J(G)/N_TOTAL) - LOG_P;
 
    MBAR = sum(N_J(G).*MOLAR_MASSES(G)) / N_TOTAL;
    Propellant.H = R * T * sum(N_J .* H_RT);
    Propellant.S_R = sum(N_J .* S_MU_R);
    Propellant.MBAR = MBAR;
    Propellant.RHO = P * MBAR / (R * T);
    Propellant.T = T;
    Propellant.CP = sum(N_J .* CP_R) * R;
    Propellant.N_CONDENSED = sum(N_J(C));
 
    % Viscosity / conductivity over the gaseous species only
    if NEED_MU == 1
        ng = numel(G);
        MUs = zeros(NUM_SPECIES, 1);
        ks = zeros(NUM_SPECIES, 1);
        
        lnT = log(T);
        invT = 1 / T;
        invT2 = 1 / (T^2);
        
        % 1. Evaluate NASA polynomials for Viscosity and Conductivity
        for idx = 1:ng
            i = G(idx);
            
            % Find the correct Viscosity interval for temperature T
            v_int = find(T <= V_RANGES_MAX(i, :), 1, 'first');
            if isempty(v_int), v_int = size(V_RANGES_MAX, 2); end % Fallback to highest
            
            rowV = (i - 1) * size(V_RANGES_MAX, 2) + v_int;
            
            % Calculate viscosity
            A_V = V_COEFFS(rowV, 1);
            B_V = V_COEFFS(rowV, 2);
            C_V = V_COEFFS(rowV, 3);
            D_V = V_COEFFS(rowV, 4);
            
            MUs(i) = exp(A_V * lnT + B_V * invT + C_V * invT2 + D_V) * 1e-7; % Convert microPoise to Pa*s
            
            % Find the correct Conductivity interval for temperature T
            c_int = find(T <= C_RANGES_MAX(i, :), 1, 'first');
            if isempty(c_int), c_int = size(C_RANGES_MAX, 2); end
            
            row_C = (i - 1) * size(C_RANGES_MAX, 2) + c_int;
            
            % Calculate conductivity
            A_C = C_COEFFS(row_C, 1);
            B_C = C_COEFFS(row_C, 2);
            C_C = C_COEFFS(row_C, 3);
            D_C = C_COEFFS(row_C, 4);
            
            ks(i) = exp(A_C * lnT + B_C * invT + C_C * invT2 + D_C) * 1e-4; % W / (m K)
        end
        
        % Wilke interaction matrix (Phi) for mixture calculations
        PHI = zeros(ng);
        for a = 1:ng
            for b = 1:ng
                i = G(a); j = G(b);
                NUMER = (1 + sqrt(MUs(i)/MUs(j))*(MOLAR_MASSES(j)/MOLAR_MASSES(i))^0.25)^2;
                DENOM = sqrt(8)*sqrt(1 + MOLAR_MASSES(i)/MOLAR_MASSES(j));
                PHI(a,b) = NUMER/DENOM;
            end
        end
        
        % Calculate final mixture properties
        Xg = N_J(G) ./ N_TOTAL; % Mole fractions
        Propellant.MU = sum((Xg .* MUs(G)) ./ (PHI * Xg));
        Propellant.k  = sum((Xg .* ks(G))  ./ (PHI * Xg));

    end
    
end