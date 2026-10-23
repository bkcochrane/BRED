function Ievlev = calculateIevlev(Engine, PropData)

    % Ievlev boundary-layer method from Vasiliev & Kudryavtsev
    % Osnovy teorii i rascheta ZhRD, vol 2, 4th ed, 1993. Equations referenced below.

    x = Engine.Stream.x;
    r = Engine.Stream.r;
    tempChamber = Engine.Chamber.temp;
    rSpecific = Engine.Stream.rSpecific;
    prandtl = Engine.Stream.prandtl;
    rho = Engine.Stream.rho;
    rhoChamber = rho(1);
    mu = Engine.Stream.mu;
    muChamber = mu(1);
    diameterThroat = Engine.Throat.diameter;
    v = Engine.Stream.v;
    tempAW = Engine.Wall.tempAW;
    tempWG = Engine.Wall.tempWG;

    numStations = numel(x);

    N_EXP = 0.15; % eq 11.65 power-law exponent
    expS = 1 / (1 - N_EXP); % 1.2, eq 11.83
    EPS_K = 1.0; % p0c / pk, eq 11.87

    kAdiabatic = Engine.Chamber.cp / (Engine.Chamber.cp - PropData.R / Engine.Chamber.mbar);
    rSpecificChamber = PropData.R / Engine.Chamber.mbar;

    % Reference Reynolds number, eq 11.27
    wMaxRef = sqrt(2 * kAdiabatic * rSpecificChamber * tempChamber / (kAdiabatic - 1));
    reynolds0 = rhoChamber * wMaxRef * diameterThroat / muChamber;

    % Coefficient a, eq 11.74
    aCoeff = 2^1.18 * (2/(kAdiabatic + 1))^(1/(kAdiabatic - 1)) * sqrt(2*kAdiabatic/(kAdiabatic + 1));

    nj1500 = Engine.Chamber.nj1500;
    rSpecific1500K = Engine.Chamber.rSpecific1500K;

    % J_0g for the undissociated gas: 1500 K composition at the chamber temperature
    [~, hRTchamber] = nasa9(PropData.A_NASA, PropData.B_NASA, PropData.T_HI, tempChamber);
    h0Undissoc = PropData.R * tempChamber * sum(nj1500 .* hRTchamber);

    % Relative coordinate and diameter, normalised by the throat diameter
    xBar = x / diameterThroat;
    dBar = 2 * r / diameterThroat;

    betaVec = zeros(numStations, 1);
    beta2 = zeros(numStations, 1);
    tempWBar = zeros(numStations, 1);
    pRatio = zeros(numStations, 1);
    zT = zeros(numStations, 1);
    z = zeros(numStations, 1);

    lengthCombustion = findCombustionLength(Engine);
    iStart = find(x >= x(1) + lengthCombustion, 1, 'first');
    if isempty(iStart)
        iStart = 2;
        fprintf('[Warning] [Ievlev] Combustion length exceeds chamber length; integrating from station 2.');
    end
    iStart = max(iStart, 2);

    for i = 1:numStations
        % beta = w / w_max, eq 11.27. w_max is a single stagnation reference.
        betaVec(i) = min(v(i) / wMaxRef, 0.99);
        beta2(i) = betaVec(i)^2;

        % T_e = R * T0inf / R_1500
        tempE = rSpecific(i) * tempChamber / rSpecific1500K;
        tempWBar(i) = min(max(tempWG(i) / tempE, 0.05), 0.95);

        pRatio(i) = Engine.Stream.p(i) / Engine.Chamber.P;
    end

    % zT from eq 11.83
    for s = 1:numStations
        if s > iStart
            idx = iStart:s;
            den1 = (1 - beta2(idx) ./ (2 * (1 + tempWBar(idx)))).^1.52;
            den2 = (1 - 9 * beta2(idx) ./ (4 * (3 + tempWBar(idx)))).^0.18;

            integrand = pRatio(idx) .* dBar(idx).^expS .* betaVec(idx) ./ (den1 .* den2);
            integral = trapz(xBar(idx), integrand);
        else
            integral = 0;
        end

        zT(s) = 1.2 * 2^1.88 * reynolds0 ...
                / (dBar(s)^expS * (1 + tempWBar(s))^1.52 * (3 + tempWBar(s))^0.18) ...
                * integral;
    end

    % alphaT from eq 11.63, heat flux from eq 11.72
    hG = zeros(numStations, 1);
    qVec = zeros(numStations, 1);

    for i = 1:numStations
        % zT/z, eq 11.76
        denomzTz = max(1 - tempWBar(i) + 0.1 * beta2(i), 1e-3);
        numerzTz = 1 - beta2(i) + beta2(i) * (1 - 0.08696 * (1 - beta2(i)) / denomzTz);
        zTz = (1.769 * max(numerzTz / denomzTz, 1e-3))^0.54;

        z(i) = abs(zT(i) / zTz);

        % eq 11.61, the recovery-to-stagnation enthalpy ratio
        baseNum = max(1 - 0.21 * ((1 - prandtl(i)) / prandtl(i)^(4/3)) ...
                      * beta2(i) / max(1 - tempWBar(i), 1e-3), 1e-3);

        numerHG = (1 / zTz)^(0.089 * prandtl(i)^-0.56) * baseNum^0.9225;
        denomHG = (307.8 + 54.8 * log10(prandtl(i) / 19.5)^2) * prandtl(i)^0.45 * z(i)^0.08 - 650;

        if denomHG <= 0
            qVec(i) = NaN;
            hG(i) = NaN;
            continue;
        end

        alphaT = numerHG / denomHG;

        % J_st: products at the wall temperature with the 1500 K composition
        [~, hRT] = nasa9(PropData.A_NASA, PropData.B_NASA, PropData.T_HI, tempWG(i));
        hWall = PropData.R * tempWG(i) * sum(nj1500 .* hRT);
        deltaH = h0Undissoc - hWall;

        % phi1, eq 11.70
        phi1 = (1 - beta2(i)) / ( (1 - beta2(i)/(2*(1 + tempWBar(i))))^0.82 ...
                                * (1 - 9*beta2(i)/(4*(3 + tempWBar(i))))^0.18 );

        % eq 11.72
        qVec(i) = aCoeff * alphaT * phi1 * Engine.Chamber.P * EPS_K * deltaH ...
                  / ( dBar(i)^2 * (rSpecificChamber * tempChamber)^0.5 ...
                      * (1 + tempWBar(i))^0.82 * (3 + tempWBar(i))^0.18 );

        hG(i) = qVec(i) / max(tempAW(i) - tempWG(i), 200);
        
    end

    iValid = find(isfinite(qVec) & (1:numStations).' > iStart, 1, 'first');
    if isempty(iValid)
        error('calculateIevlev: no valid station downstream of the combustion zone');
    end
    iAnchor = find(x >= x(iValid) + 2 * diameterThroat, 1, 'first');
    if isempty(iAnchor), iAnchor = numStations; end
    qVec(1:iAnchor) = qVec(iAnchor);
    hG(1:iAnchor) = hG(iAnchor);

    assert(all(isfinite(qVec)), 'calculateIevlev: %d non-finite flux values', sum(~isfinite(qVec)));

    Ievlev.q = qVec;
    Ievlev.hG = hG;
    Ievlev.zT = zT;
    Ievlev.beta2 = beta2;
    Ievlev.rSpecific1500K = rSpecific1500K;

    % Better tauW calculation for estimating nozzle efficiency eqs 11.39 & 11.25
    alpha = 0.377 .* z.^(-0.224) + 3.966e-4;
    tauW = alpha .* rho .* wMaxRef.^2;
    
    Ievlev.tauW = tauW;
    
end