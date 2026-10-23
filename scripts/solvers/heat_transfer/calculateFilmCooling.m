function Film = calculateFilmCooling(Engine, PropData, filmCoolingFraction)

    % Film cooling method from Vasiliev & Kudryavtsev
    % Osnovy teorii i rascheta ZhRD, vol 2, 4th ed, 1993. Equations referenced below.

    % The Russian subscripts have been westernised (Russian z looks like a 3)

    x = Engine.Stream.x;

    alpha_k = Engine.Wall.hG;
    T_r0 = Engine.Wall.tempAW;

    T_H = PropData.Fuel.T_SUPPLY;
    T_S = PropData.Fuel.TEMP_BOIL;
    c_X = PropData.Fuel.CP;
    Q_S = PropData.Fuel.QS;
    mu_X = PropData.Fuel.MU;
    T_cp = 0.5 * (T_H + T_S);

    mdot_3 = Engine.Fuel.mdot * filmCoolingFraction;
    mbar_3 = mdot_3 / Engine.Specs.mdot;
    mbar_cT = 0.15; % Near wall layer flow as a fraction of the chamber flow, yet to be linked to injector

    D = Engine.Chamber.diameter;

    Re_3 = mdot_3 / (pi * D * mu_X);

    % Film stability coefficient, eta, approximation from figure 12.24
    eta = 0.48 + 0.52 * sech(0.004 * Re_3 ^ 0.81);
    
    % Evapouration length of the liquid film, eq 12.38
    l_X = eta * (mdot_3 / (pi * D)) * (c_X * (T_S - T_H) / (alpha_k(1) * (T_r0(1) - T_cp)) + Q_S / (alpha_k(1) * (T_r0(1) - T_S)));

    % Mixing completeness, eq 12.46 - 12.48
    H_cT = Engine.Chamber.diameter * 0.01; % Distance from the last ring of injectors, not yet linked
    xbar = max(0, (x - x(1) - l_X) ./ H_cT);
    K = 0.1e-2; % Turbulent mixing constant, between 0.05e-2 and 0.2e-2
    M = K * mbar_cT / mbar_3; % Mixing coefficient
    
    % Turbulent mixing efficiency coefficient eq 12.45
    xi = 1 - exp(-M .* xbar.^2);
    
    % Limit of the mixing ratio
    k_m0 = 1.5; % Depends on the injector, not yet implemented, but 1.5 is a good approximation
    % Eq 12.44
    k_mpacy = k_m0 .* xi ./ (1 + (1 + k_m0) * mbar_3 / mbar_cT);

    if ~isfield(Engine.Cooling, 'Film') || ~isfield(Engine.Cooling.Film, 'STable')
        tabulateS(Engine.Chamber.P, Engine.Specs.OF_RATIO_TARGET)
    end

    SRatio = findSRatio(k_mpacy, k_m0);
    SRatio = max(SRatio, 0);

    Film.SRatio = SRatio;
    Film.lengthFilm = l_X;
    Film.k_mpacy = k_mpacy;
    Film.STable = Engine.Cooling.Film.STable;

    function tabulateS(P_CHAMBER, OF_RATIO_TARGET)

        % Tabulate the S function (eq 11.94) against mixture ratio, for the
        % curtain cooling flux rescaling (eq 12.49).
        %
        % S depends only on the propellant pair and the mixture ratio, so this is
        % built once per propellant combination and cached. The wall temperature
        % terms in 11.94 cancel when the ratio S/S0 is taken, so they are
        % omitted here
    
        OF_MIN = 0.1;
        OF_MAX = OF_RATIO_TARGET * 2;
        NUM_POINTS = 30;
        TEMP_WALL_REF = 1000;

        ofVector = flipud(logspace(log10(OF_MIN), log10(OF_MAX), NUM_POINTS).');
        sVector = zeros(NUM_POINTS, 1);

        tempSeed = 3000;
        njSeed = repmat(0.1 / PropData.NUM_SPECIES, PropData.NUM_SPECIES, 1);
        if isfield(PropData, 'CONDENSED_LIST')
            njSeed(PropData.CONDENSED_LIST) = 0;
        end

        for i = 1:NUM_POINTS
            [hReactants, bCirc] = chemistrySolve(ofVector(i), PropData);
            PropDataLocal = PropData;
            PropDataLocal.bCirc = bCirc;

            [tempGas, nj, Propellant] = equilibriumSolve('HP', hReactants, P_CHAMBER, ...
                                                          PropDataLocal, tempSeed, njSeed, 1);

            % J_ct: products at the reference wall temperature, composition
            % frozen at the gas state
            [~, hRT] = nasa9(PropData.A_NASA, PropData.B_NASA, PropData.T_HI, TEMP_WALL_REF);
            hWall = PropData.R * TEMP_WALL_REF * sum(nj .* hRT);

            % mu referenced to 1000 K, eq 11.18 scaling
            muRef = Propellant.mu * (TEMP_WALL_REF / tempGas)^0.7;

            % Reduced group: the (1+Tbar) and (3+Tbar) terms cancel in S/S0
            sVector(i) = (Propellant.h - hWall) * muRef^0.15 ...
                       / ( Propellant.rSpecific^0.425 * tempGas^0.320 );

            if sVector(i) > 0 && isfinite(sVector(i))
                tempSeed = tempGas; % only seed from a solve that worked
                njSeed = nj;
            else
                sVector(i) = NaN; % trimmed below
            end
        end

        % Back to ascending, and drop the points where the chemistry failed
        ofVector = flipud(ofVector);
        sVector = flipud(sVector);
        valid = isfinite(sVector);
        ofVector = ofVector(valid);
        sVector = sVector(valid);

        Engine.Cooling.Film.STable = [ofVector, sVector];
    
    end
    
    function SRatio = findSRatio(k_mpacy, k_m0)
        % S/S0 by interpolation. Held flat below the tabulated minimum
        % where, according to the book, the chemistry is not trustworthy.
        STable = Engine.Cooling.Film.STable;
        ofVector = STable(:,1);
        sVector = STable(:,2);
        k_mpacy = max(k_mpacy, ofVector(1));
        sLocal = interp1(ofVector, sVector, k_mpacy, 'pchip');
        sRef = interp1(ofVector, sVector, k_m0, 'pchip');
        SRatio = sLocal ./ sRef;
    end

end