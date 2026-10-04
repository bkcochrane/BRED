function Engine = applySeparation(Engine)

    P = Engine.STREAM_P;
    M = Engine.STREAM_M;
    X = Engine.STREAM_X;
    R = Engine.STREAM_R;
    P_AMBIENT = Engine.P_AMBIENT;
    P_WALL = P;
    AREA = pi * R.^2;

    P_PRE = P;
    P_POST = nan(numel(X),1);
    Engine.SEPARATES = false;
    Engine.X_SEP = NaN;

    if P_AMBIENT > 0
        div = find(X > 0);
        for k = 1:numel(div)
            i = div(k);
            if M(i) <= 1, continue; end

            % Schmucker separation criterion
            P_sep = (1.88*M(i) - 1)^(-0.64) * P_AMBIENT;

            if P(i) <= P_sep
                
                Engine.SEPARATES = true;
                Engine.X_SEP = X(i);
                Engine.M_SEP = M(i);
                Engine.P_SEP = P(i);
                
                P_PRE(i+1:end) = NaN;
                % Initial jump in pressure across shock
                P_JUMP = 0.9 * P_AMBIENT;
                
                % Subsonic diffusion
                AREA_sep = AREA(i);

                for j = i:numel(X)

                    frac = (AREA(j) - AREA_sep) / (AREA(end) - AREA_sep);

                    P_POST(j) = P_JUMP + (P_AMBIENT - P_JUMP) * frac;

                end
                
                P_WALL(i) = NaN;
                Engine.P_SHOCK_LO = P(i);
                Engine.P_SHOCK_HI = P_JUMP;

                break;

            end
        end
    end

    Engine.STREAM_P_PRE  = P_PRE;
    Engine.STREAM_P_POST = P_POST;
    
end