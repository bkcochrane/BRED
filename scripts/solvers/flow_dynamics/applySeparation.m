function Engine = applySeparation(Engine)

    p = Engine.Stream.p;
    mach = Engine.Stream.mach;
    x = Engine.Stream.x;
    r = Engine.Stream.r;
    P_AMBIENT = Engine.Specs.P_AMBIENT;
    pWall = p;
    area = pi * r.^2;

    pAttached = p;
    pPostSep = nan(numel(x),1);
    Engine.Stream.separates = false;
    Engine.Stream.xSep = NaN;

    if P_AMBIENT > 0
        div = find(x > 0);
        for k = 1:numel(div)
            i = div(k);
            if mach(i) <= 1, continue; end

            % Schmucker separation criterion
            pSep = (1.88 * mach(i) - 1)^(-0.64) * P_AMBIENT;

            if p(i) <= pSep
                
                Engine.Stream.separates = true;
                Engine.Stream.xSep = x(i);
                Engine.Stream.machSep = mach(i);
                Engine.Stream.pSep = p(i);
                
                pAttached(i+1:end) = NaN;
                % Initial jump in pressure across shock
                pJump = 0.9 * P_AMBIENT;
                
                % Subsonic diffusion
                areaSep = area(i);

                for j = i:numel(x)

                    frac = (area(j) - areaSep) / (area(end) - areaSep);

                    pPostSep(j) = pJump + (P_AMBIENT - pJump) * frac;

                end
                
                pWall(i) = NaN;
                Engine.Stream.pShockLow = p(i);
                Engine.Stream.PShockHi = pJump;

                break;

            end
        end
    end

    Engine.Stream.pAttached  = pAttached;
    Engine.Stream.pPostSep = pPostSep;
    
end