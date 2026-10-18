function [temp, nj, Propellant, piI] = findEquilibrium(mode, target, p, PropData, temp, nj, needMu, isFrozen)

    if nargin < 8; isFrozen = false; end

    % HP mode solves for fixed enthalpy (chamber pressure)
    % HP targets H_REACTANTS
    %
    % SP mode solves for insentropic expansion
    % SP targets S0_R

    A_IJ = PropData.A_IJ;
    bCirc = PropData.bCirc;
    A_NASA = PropData.A_NASA;
    B_NASA = PropData.B_NASA;
    T_HI = PropData.T_HI;
    MOLAR_MASSES = PropData.M;
    R = PropData.R;
    V_COEFFS = PropData.V_COEFFS;
    V_RANGES_MAX = PropData.V_RANGES_MAX;
    C_COEFFS = PropData.C_COEFFS;
    C_RANGES_MAX = PropData.C_RANGES_MAX;
    logP = log(p / 1e5);
 
    NUM_SPECIES = numel(nj);
    NUM_ELEM = size(A_IJ, 1);
 
    % Partition species into gas (G) and active condensed (C)
    if isfield(PropData, 'CONDENSED_LIST') && any(PropData.CONDENSED_LIST)
        condensedAll = PropData.CONDENSED_LIST(:).';
    else
        condensedAll = [];
    end

    if isfield(PropData, 'ACTIVE_CONDENSED') && any(PropData.ACTIVE_CONDENSED)
        condensedActive = PropData.ACTIVE_CONDENSED(:).';
    else
        condensedActive = [];
    end

    gasIdx = 1:NUM_SPECIES;
    gasIdx(condensedAll) = [];
    numC = numel(condensedActive);
    nj(gasIdx) = max(nj(gasIdx), 1e-30);
 
    MAX_ITERATIONS = 100;
 
    for iter = 1:MAX_ITERATIONS
        
        % Skip matrix if equilibrium is frozen
        if isFrozen && strcmp(mode, 'SP')

            [cpR, ~, sR, ~] = nasa9(A_NASA, B_NASA, T_HI, temp);
            nTotal = max(sum(nj(gasIdx)), 1e-30);
            
            sMuR = sR;
            sMuR(gasIdx) = sR(gasIdx) - log(max(nj(gasIdx), 1e-250) / nTotal) - logP;
            sMixR = sum(nj .* sMuR);
            
            % 1D Newton step for Temperature
            deltaLnT = (target - sMixR) / sum(nj .* cpR);
            temp = max(temp * exp(deltaLnT), 150);
            
            if abs(deltaLnT) < 1e-5
                piI = zeros(NUM_ELEM, 1);
                break;
            end

            continue;
        end
 
        nTotal = max(sum(nj(gasIdx)), 1e-30); % Gas - moles only
        [cpR, hRT, ~, gRT] = nasa9(A_NASA, B_NASA, T_HI, temp);
 
        % Chemical potentials: gases carry the ln terms, condensed do not
        muRT = gRT; % Condensed: mu = G_RT
        muRT(gasIdx) = gRT(gasIdx) + log(nj(gasIdx) / nTotal) + logP; % Gases add ln(x)+ln P
 
        nH  = nj .* hRT;
        nMu = nj .* muRT;
        nMu(isnan(nMu) | ~isfinite(nMu)) = 0.0;
 
        sz = NUM_ELEM + numC + 2; % Elements + condensed + dlnN + dlnT
        A  = zeros(sz, sz);
        B  = zeros(sz, 1);
 
        iN = NUM_ELEM + numC + 1; % row/col index of dln N
        iTemp = NUM_ELEM + numC + 2; % row/col index of dln T
 
        % Element rows (gas contribution + condensed columns)
        A(1:NUM_ELEM, 1:NUM_ELEM) = (A_IJ(:,gasIdx) .* nj(gasIdx).') * A_IJ(:,gasIdx).';
        A(1:NUM_ELEM, 1:NUM_ELEM) = A(1:NUM_ELEM, 1:NUM_ELEM) + eye(NUM_ELEM) * 1e-15;
        for q = 1:numC
            A(1:NUM_ELEM, NUM_ELEM+q) = A_IJ(:, condensedActive(q)); % Condensed column = atom vector
        end
        A(1:NUM_ELEM, iN) = A_IJ(:,gasIdx) * nj(gasIdx);
        A(1:NUM_ELEM, iTemp) = A_IJ(:,gasIdx) * nH(gasIdx);
        B(1:NUM_ELEM) = bCirc - A_IJ*nj + A_IJ(:,gasIdx) * nMu(gasIdx);
 
        % Condensed rows: phase equilibrium sum a_ij pi_i + (H_c/RT)dlnT = mu_c/RT
        for q = 1:numC
            j = condensedActive(q);
            A(NUM_ELEM+q, 1:NUM_ELEM) = A_IJ(:,j).';
            A(NUM_ELEM+q, iTemp) = hRT(j);
            B(NUM_ELEM+q) = muRT(j);
        end
 
        % Total-moles row (gases only)
        A(iN, 1:NUM_ELEM) = (A_IJ(:,gasIdx) * nj(gasIdx)).';
        A(iN, iN) = sum(nj(gasIdx)) - nTotal;
        A(iN, iTemp) = sum(nH(gasIdx));
        B(iN) = nTotal - sum(nj(gasIdx)) + sum(nMu(gasIdx));
 
        % Energy / entropy row 
        if strcmp(mode,'HP')
            hMixRT = sum(nH); 
            A(iTemp, 1:NUM_ELEM) = A(1:NUM_ELEM, iTemp).';
            for q = 1:numC, A(iTemp, NUM_ELEM+q) = hRT(condensedActive(q)); end
            A(iTemp, iN) = sum(nH(gasIdx));
            
            % Substitution artifacts must only sum over Gas (G)
            A(iTemp, iTemp) = sum(nj .* cpR) + sum(nH(gasIdx) .* hRT(gasIdx)); 
            B(iTemp) = target/(R*temp) - hMixRT + sum(nMu(gasIdx) .* hRT(gasIdx)); 
        else
            [~, ~, S_R_full, ~] = nasa9(A_NASA, B_NASA, T_HI, temp);
            sMuR = S_R_full - muRT + gRT; 
            sMixR = sum(nj .* sMuR);
            A(iTemp, 1:NUM_ELEM) = (A_IJ(:,gasIdx) * (nj(gasIdx) .* sMuR(gasIdx))).';
            for q = 1:numC, A(iTemp, NUM_ELEM+q) = sMuR(condensedActive(q)); end
            
            % Substitution artifacts must only sum over Gas (G)
            A(iTemp, iN) = sum(nj(gasIdx) .* sMuR(gasIdx)); 
            A(iTemp, iTemp) = sum(nj .* cpR) + sum(nj(gasIdx) .* hRT(gasIdx) .* sMuR(gasIdx)); 
            B(iTemp) = target - sMixR + (nTotal - sum(nj(gasIdx))) + sum(nj(gasIdx) .* sMuR(gasIdx) .* muRT(gasIdx)); 
        end
        
        % Solve
        if rcond(A) < 1e-15
            % Use Moore-Penrose pseudo-inverse for singular/badly scaled matrices
            x = pinv(full(A)) * B; 
        else
            x = A \ B;
        end

        piI = x(1:NUM_ELEM);
        deltaNumC = x(NUM_ELEM + (1:numC)); % Linear update for condensed
        deltaLnN = x(iN);
        deltaLnTemp = x(iTemp);
 
        % Gas log-updates derived from the multipliers
        deltaLnNJ = zeros(NUM_SPECIES,1);
        deltaLnNJ(gasIdx) = A_IJ(:,gasIdx).'*piI + deltaLnN + hRT(gasIdx)*deltaLnTemp - muRT(gasIdx);
 
        % Convergence
        convGas = max(abs(nj(gasIdx).*deltaLnNJ(gasIdx)))/nTotal;
        convCondensed = 0; if numC>0, convCondensed = max(abs(deltaNumC))/max(nTotal,eps); end
        if convGas < 0.5e-5 && convCondensed < 0.5e-5 && ...
            abs(deltaLnN) < 0.5e-5 && abs(deltaLnTemp) < 1e-4
            break;
        end
 
        % Step limiting
        big = (nj(gasIdx)/nTotal) > 1e-8;
        gBig = gasIdx(big);
        lambda1 = 2 / max([5*abs(deltaLnTemp), 5*abs(deltaLnN), ...
                      max(abs(deltaLnNJ(gBig)))]);
        lambda = min([1, lambda1]);
 
        % Apply
        deltaTempSafe = max(-10, min(10, lambda*deltaLnTemp));
        gasUpdateSafe = max(-50, min(50, lambda*deltaLnNJ(gasIdx)));
        
        tempNew = temp * exp(deltaTempSafe);
        temp = max(150, min(20000, tempNew)); 
        
        nj(gasIdx) = max(nj(gasIdx) .* exp(gasUpdateSafe), 1e-25);
        
        if numC > 0
            nj(condensedActive) = max(nj(condensedActive) + lambda * deltaNumC, 0); % Additive, can't be negative
        end
        
    end
 
    % Derived properties
    [~, hRT, sR, gRT] = nasa9(A_NASA, B_NASA, T_HI, temp);
    nTotal = sum(nj(gasIdx));
    sMuR = sR;
    sMuR(gasIdx) = sR(gasIdx) - log(nj(gasIdx)/nTotal) - logP;
 
    mbar = sum(nj(gasIdx) .* MOLAR_MASSES(gasIdx)) / nTotal;
    Propellant.h = R * temp * sum(nj .* hRT);
    Propellant.sR = sum(nj .* sMuR);
    Propellant.gRT = gRT;
    Propellant.mbar = mbar;
    Propellant.rho = p * mbar / (R * temp);
    Propellant.temp = temp;
    Propellant.cp = sum(nj .* cpR) * R;
    Propellant.rSpecific = sum(nj) * R;
    Propellant.gamma = Propellant.cp / (Propellant.cp - Propellant.rSpecific);
    Propellant.numCondensed = sum(nj(condensedActive));
    Propellant.massFraction = nj .* MOLAR_MASSES;
    Propellant.massPercentage = Propellant.massFraction * 100;
 
    % Viscosity / conductivity over the gaseous species only
    if needMu == 1
        numG = numel(gasIdx);
        mu = zeros(NUM_SPECIES, 1);
        k = zeros(NUM_SPECIES, 1);
        
        lnTemp = log(temp);
        invTemp = 1 / temp;
        invTemp2 = 1 / (temp^2);
        
        % 1. Evaluate NASA polynomials for Viscosity and Conductivity
        for idx = 1:numG
            i = gasIdx(idx);
            
            % Find the correct Viscosity interval for temperature T
            v_int = find(temp <= V_RANGES_MAX(i, :), 1, 'first');
            if isempty(v_int), v_int = size(V_RANGES_MAX, 2); end % Fallback to highest
            
            rowV = (i - 1) * size(V_RANGES_MAX, 2) + v_int;
            
            % Calculate viscosity
            aV = V_COEFFS(rowV, 1);
            bV = V_COEFFS(rowV, 2);
            cC = V_COEFFS(rowV, 3);
            dV = V_COEFFS(rowV, 4);
            
            mu(i) = exp(aV * lnTemp + bV * invTemp + cC * invTemp2 + dV) * 1e-7; % Convert microPoise to Pa*s
            
            % Find the correct Conductivity interval for temperature T
            c_int = find(temp <= C_RANGES_MAX(i, :), 1, 'first');
            if isempty(c_int), c_int = size(C_RANGES_MAX, 2); end
            
            row_C = (i - 1) * size(C_RANGES_MAX, 2) + c_int;
            
            % Calculate conductivity
            aC = C_COEFFS(row_C, 1);
            bC = C_COEFFS(row_C, 2);
            cC = C_COEFFS(row_C, 3);
            dC = C_COEFFS(row_C, 4);
            
            k(i) = exp(aC * lnTemp + bC * invTemp + cC * invTemp2 + dC) * 1e-4; % W / (m K)
        end
        
        % Wilke interaction matrix (Phi) for mixture calculations
        phi = zeros(numG);
        for a = 1:numG
            for b = 1:numG
                i = gasIdx(a); j = gasIdx(b);
                numer = (1 + sqrt(mu(i)/mu(j))*(MOLAR_MASSES(j)/MOLAR_MASSES(i))^0.25)^2;
                denom = sqrt(8)*sqrt(1 + MOLAR_MASSES(i)/MOLAR_MASSES(j));
                phi(a,b) = numer/denom;
            end
        end
        
        % Calculate final mixture properties
        xg = nj(gasIdx) ./ nTotal; % Mole fractions
        Propellant.mu = sum((xg .* mu(gasIdx)) ./ (phi * xg));
        Propellant.k = sum((xg .* k(gasIdx))  ./ (phi * xg));

    end
    
end