function Engine = calculateStreamwise(Engine, PropData, region)
    
    switch region
        case 'diverging', x = Engine.Geometry.xDiverging(:); r = Engine.Geometry.rDiverging(:);
        case 'full', x = Engine.Geometry.xCoordinates(:); r = Engine.Geometry.rCoordinates(:);
        otherwise, error('region must be ''diverging'' or ''full''');
    end
    
    [~, PropData.bCirc] = chemistrySolve(Engine.Specs.ofRatio, PropData);
    h0 = Engine.Chamber.h0;
    s0R = Engine.Chamber.s0R;

    % Sort ascending and drop the duplicate throat point
    [x, iu] = unique(x);
    r = r(iu);
    numX = numel(x);

    rhos = zeros(numX,1);
    vs = rhos;
    mus = rhos;
    ps = rhos;
    temps = rhos;
    machs = rhos;
    ks = rhos;
    cps = rhos;
    rSpecifics = rhos;
    gammas = rhos;
    massFractions = zeros(numX, PropData.NUM_SPECIES);
    massPercentages = massFractions;

    areaScale = Engine.Efficiency.combustion * Engine.Efficiency.cdThroat;
    
    % Throat - reference point for both marches
    [~, iThroat] = min(abs(x));
    effectiveThroatArea = (pi * areaScale) * (r(iThroat)^2);
    [tempThroat, njThroat] = store(iThroat, Engine.Throat.p, Engine.Throat.temp, Engine.Throat.nj, false);
    
    % March from throat
    march(iThroat+1 : 1 : numX, 'supersonic', tempThroat, njThroat);
    march(iThroat-1 : -1 : 1, 'subsonic', tempThroat, njThroat);
    
    Engine.Stream.x = x;
    Engine.Stream.r = r;
    Engine.Stream.rho = rhos;
    Engine.Stream.v = vs;
    Engine.Stream.mu = mus;
    Engine.Stream.p = ps;
    Engine.Stream.temp = temps;
    Engine.Stream.mach = machs;
    Engine.Stream.k = ks;
    Engine.Stream.cp = cps;
    Engine.Stream.rSpecific = rSpecifics;
    Engine.Stream.gamma = gammas;
    Engine.Stream.massFraction = massFractions;
    Engine.Stream.massPercentage = massPercentages;
    
    Engine = applySeparation(Engine);
    
    % ------------------------------------------------------------------------
    % March function
    % ------------------------------------------------------------------------
    function march(idx, branch, temp, nj)
        
        pPrevious = Engine.Throat.p;
        mdot = Engine.Specs.mdot;
        
        areaLocalAll = (pi * areaScale) .* (r(idx).^2);
        
        isSupersonic = strcmp(branch, 'supersonic');
        freezeThreshold = 1.3 * effectiveThroatArea;
        
        njFrozen = nj;
        isFrozen = false;
        
        for k = 1:numel(idx)
            i = idx(k);
            areaLocal = areaLocalAll(k);
            
            if ~isFrozen && isSupersonic && (areaLocal > freezeThreshold)
                isFrozen = true;
            end
            
            if isFrozen
                p = solvePressure(s0R, h0, Engine, PropData, temp, njFrozen, mdot, areaLocal, pPrevious, branch, true);
                [temp, njFrozen] = store(i, p, temp, njFrozen, true);
            else
                p = solvePressure(s0R, h0, Engine, PropData, temp, nj, mdot, areaLocal, pPrevious, branch, false);
                [temp, nj] = store(i, p, temp, nj, false);
                njFrozen = nj;
            end
            
            pPrevious = p;
        end
        
    end

    % ------------------------------------------------------------------------
    % Solve and update, then seed value to next point
    % ------------------------------------------------------------------------
    function [temp, nj] = store(i, p, tempGuess, njGuess, isFrozen)
        
        [temp, nj, Propellant] = equilibriumSolve('SP', s0R, p, PropData, tempGuess, njGuess, 1, isFrozen);
        
        rhos(i) = Propellant.rho;
        vs(i) = sqrt(max(2*(h0 - Propellant.h), 0));
        mus(i) = Propellant.mu;
        ps(i) = p;
        temps(i) = Propellant.temp;
        ks(i) = Propellant.k;
        cps(i) = Propellant.cp;
        rSpecifics(i) = Propellant.rSpecific;
        gammas(i) = Propellant.gamma;
        machs(i) = vs(i) / sqrt(Propellant.gamma * Propellant.rSpecific * Propellant.temp);
        massFractions(i, :) = Propellant.massFraction(:)';
        massPercentages(i, :) = Propellant.massPercentage(:)';

    end

end