function ModifiedBartz = calculateModifiedBartz(Engine)

    % Taken from Di Matteo et al. (2010)
    % SEMI-EMPIRICAL HEAT TRANSFER CORRELATIONS IN COMBUSTION CHAMBERS FOR TRANSIENT SYSTEM MODELLING
    %
    % Applies scaling factors to the orginal Bartz equations to better match experimental data
    
    tempCore = Engine.Stream.temp;
    tempWG = Engine.Wall.tempWG;
    tempRef = 0.5 .* (tempWG + tempCore);
    muRef = Engine.Stream.mu .* (tempRef ./ tempCore).^0.7;
    cpRef = Engine.Stream.cp;
    kRef = Engine.Stream.k .* (tempRef ./ tempCore).^0.7;
    rho = Engine.Stream.rho;
    v = Engine.Stream.v;

    mdot = Engine.Specs.mdot;

    xCoords = Engine.Stream.x;
    r = Engine.Stream.r;
    area = pi .* r.^2;
    diameterThroat = Engine.Throat.diameter;
    radiusCurvature = 0.382 * diameterThroat / 2;
    [~, iThroat] = min(r);
    x = xCoords - xCoords(1);
    xThroat = x(iThroat);
    

    tempAW = Engine.Wall.tempAW;
    
    % ---------------------
    % Original Bartz, eq 10
    % ---------------------

    stantonBartz = 0.026 .* (muRef.^0.2 ./ cpRef.^0.6) .* (kRef ./ muRef).^0.6 ...
                    .* mdot.^(-0.2) .* area.^0.1 .* ((pi / 4) / (radiusCurvature * diameterThroat)).^0.1;
    
    % -------------------------
    % Correction factors, eq 11
    % -------------------------

    % Tuning coefficients
    a = 0;
    b = 0;
    
    % Temperature correction factor
    Kt = (tempAW ./ tempRef).^a;

    % Boundary layer growth correction factor
    Kx = max(0.01, (x ./ xThroat)).^b;
    
    % ---------------------
    % Modified Bartz, eq 12
    %----------------------

    stantonModified = stantonBartz .* Kt .* Kx;
    
    % -------------------------------
    % Combustion length factor, eq 14
    % -------------------------------

    lengthCombustion = findCombustionLength(Engine);
    combustionFactor = min(1, 0.25 .* atan(7 .* (x ./ lengthCombustion - 0.63)) + 0.7);
    stantonCombustion = stantonModified .* combustionFactor;
    
    % -------------------------------
    % Flow acceleration factor, eq 15
    % -------------------------------
    drdx = gradient(r, x);

    Kacc = sqrt(max(1 - abs(drdx), 0));
    stantonFinal = stantonCombustion .* Kacc;

    % Convert to heat transfer coefficient, eq 9
    hG = stantonFinal .* rho .* v .* cpRef;

    ModifiedBartz.hG = hG;
    
end