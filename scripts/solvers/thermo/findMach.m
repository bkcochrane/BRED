function [mach, a] = findMach(p, s0R, h0, PropData, temp0, nj0)

    [temp, nj, Propellant] = equilibriumSolve('SP', s0R, p, PropData, temp0, nj0, 0);

    v = sqrt(max(2*(h0 - Propellant.h), 0));
    
    A_IJ = PropData.A_IJ;
    NUM_ELEM = size(A_IJ, 1);
    
    % Separate gaseous components
    if isfield(PropData, 'CONDENSED_LIST')
        CANDIDATES = PropData.CONDENSED_LIST(:).';
    else
        CANDIDATES = [];
    end
    candidates = CANDIDATES(nj(CANDIDATES) > 1e-10);
    gasIdx = setdiff(1:PropData.NUM_SPECIES, CANDIDATES);
    numCandidates = numel(candidates);
    
    nG = sum(nj(gasIdx));
    
    % Build Jacobian matrix
    iNum = NUM_ELEM + numCandidates + 1;
    m = zeros(iNum, iNum);
    m(1:NUM_ELEM, 1:NUM_ELEM) = (A_IJ(:, gasIdx) .* nj(gasIdx).') * A_IJ(:, gasIdx).';
    m(1:NUM_ELEM, 1:NUM_ELEM) = m(1:NUM_ELEM, 1:NUM_ELEM) + eye(NUM_ELEM) * 1e-15;

    for q = 1:numCandidates
        m(1:NUM_ELEM, NUM_ELEM + q) = A_IJ(:, candidates(q));
        m(NUM_ELEM + q, 1:NUM_ELEM) = A_IJ(:, candidates(q)).';
    end

    m(1:NUM_ELEM, iNum) = A_IJ(:, gasIdx) * nj(gasIdx);
    m(iNum, 1:NUM_ELEM) = m(1:NUM_ELEM, iNum).';
    m(iNum, iNum) = 0;
    
    % --------------------------------------------------------------------
    % Derivative with respect to ln(P)
    bP = zeros(iNum, 1);
    bP(1:NUM_ELEM) = A_IJ(:, gasIdx) * nj(gasIdx);
    bP(iNum) = nG;

    xP = m \ bP;
    dLnNG_dLnP = xP(iNum);
    dLnV_dLnP_T = dLnNG_dLnP - 1;
    % --------------------------------------------------------------------
    % Derivative with respect to ln(T) at constant P
    [cpR, hRT, ~, ~] = nasa9(PropData.A_NASA, PropData.B_NASA, PropData.T_HI, temp);
    
    bTemp = zeros(iNum, 1);
    bTemp(1:NUM_ELEM) = -A_IJ(:, gasIdx) * (nj(gasIdx) .* hRT(gasIdx));
    for q = 1:numCandidates
        bTemp(NUM_ELEM + q) = -hRT(candidates(q));
    end
    bTemp(iNum) = -sum(nj(gasIdx) .* hRT(gasIdx));

    xTemp = m \ bTemp;
    piTemp = xTemp(1:NUM_ELEM);
    dnC_dLnT_P = xTemp(NUM_ELEM + 1 : NUM_ELEM + numCandidates); % Condensed derivatives
    dLnNg_dLnT = xTemp(iNum);
    dLnV_dLnT_P = dLnNg_dLnT + 1;
    % --------------------------------------------------------------------

    dlnNj_dlnT = dLnNg_dLnT + A_IJ(:, gasIdx).' * piTemp + hRT(gasIdx);

    R = PropData.R;

    % Equilibrium Cp
    cpFrozen = sum(nj .* cpR) * R;
    cpGas = sum(nj(gasIdx) .* hRT(gasIdx) .* dlnNj_dlnT) * R; 
    cpCondensed = sum(hRT(candidates) .* dnC_dLnT_P) * R;
    cpEquilibrium = cpFrozen + cpGas + cpCondensed;
    
    % Equilibrium Cv
    cvEquilibrium = cpEquilibrium + (nG * R) * (dLnV_dLnT_P^2) / dLnV_dLnP_T;
    
    % Isentropic Exponent (gamma)
    gammaEquilibrium = cpEquilibrium / cvEquilibrium;
    
    % Analytical Speed of Sound
    % a^2 = - gamma / dlnV_dlnP_T * (P/rho) --> where (P/rho) = N_g * R * T
    a2 = -(gammaEquilibrium / dLnV_dLnP_T) * nG * R * temp;
    a = sqrt(a2);
    
    mach = v / a;

end