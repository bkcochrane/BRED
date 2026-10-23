function Pavli = calculatePavli(Engine)

    % Taken from Di Matteo et al. (2010)
    % SEMI-EMPIRICAL HEAT TRANSFER CORRELATIONS IN COMBUSTION CHAMBERS FOR TRANSIENT SYSTEM MODELLING
    %
    % The Pavli correlation is a simpler Stanton-type alternative to Bartz.
    % The paper notes it under-predicts where Bartz over-predicts, as the two bracket the experimental data.

    tempCore = Engine.Stream.temp;
    tempWG = Engine.Wall.tempWG;
    tempRef = 0.5 .* (tempWG + tempCore);
    muRef = Engine.Stream.mu .* (tempRef ./ tempCore).^0.7;
    cpRef = Engine.Stream.cp;
    kRef = Engine.Stream.k .* (tempRef ./ tempCore).^0.7;
    rho = Engine.Stream.rho;
    v = Engine.Stream.v;

    xCoords = Engine.Stream.x;
    r = Engine.Stream.r;
    [~, iThroat] = min(r);
    x = xCoords - xCoords(1);
    xThroat = x(iThroat);

    tempAW = Engine.Wall.tempAW;

    diameterLocal = 2 .* r;
    reynolds = rho .* v .* diameterLocal ./ muRef;
    prandtl = muRef .* cpRef ./ kRef;

    stantonPavli = 0.023 .* reynolds.^(-0.2) .* prandtl.^(-0.6);

    % -------------------------
    % Correction factors, eq 13
    % -------------------------

    % Tuning coefficients
    e = 0;
    f = 0;

    % Temperature correction factor
    Kt = (tempAW ./ tempRef).^e;

    % Boundary layer growth correction factor
    Kx = max(0.01, (x ./ xThroat)).^f;

    stantonModified = stantonPavli .* Kt .* Kx;

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

    Pavli.hG = hG;

end