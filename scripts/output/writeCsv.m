function writeCsv(Engine)

    % Gas side
    X_HOT = Engine.STREAM_X(:) * 1000; % mm
    R_HOT = Engine.STREAM_R(:) * 1000;
    [X_HOT, iu] = unique(X_HOT);
    R_HOT = R_HOT(iu);
    s0 = [0; cumsum(hypot(diff(X_HOT), diff(R_HOT)))];
    if exist('csaps','file')
        X_HOT = csaps(s0, X_HOT, 0.999, s0);
        R_HOT = csaps(s0, R_HOT, 0.999, s0);
    end
    X_HOT = X_HOT(:);  R_HOT = R_HOT(:);

    % Channel dimensions
    t = Engine.WALL_THICKNESS * 1000; % mm
    Hs = Engine.CHANNEL_HEIGHT_STATION(iu) * 1000; % mm
    Ws = Engine.CHANNEL_WIDTH_STATION(iu)  * 1000; % mm

    % Wall normal
    dX = gradient(X_HOT);  dR = gradient(R_HOT);  L = hypot(dX, dR);
    nX = -dR ./ L;
    nR = dX ./ L;

    % Cold side (channel floor) and outer (channel roof), both normal offsets
    X_FLOOR = X_HOT + t .* nX;
    R_FLOOR = R_HOT + t .* nR;
    X_ROOF = X_HOT + (t + Hs) .* nX;
    R_ROOF = R_HOT + (t + Hs) .* nR;
    
    % Helix
    if isfield(Engine,'SPIRAL_RATE'), RATE = Engine.SPIRAL_RATE / 1000; else, RATE = 0; end
    R_MID = 0.5 * (R_FLOOR + R_ROOF);
    THETA = RATE * (X_HOT - min(X_HOT));

    % Resample for smooth curve
    NPTS = 300;
    s = [0; cumsum(hypot(diff(X_HOT), diff(R_HOT)))];
    sf = linspace(0, s(end), NPTS).';
    [s_unique, unique_idx] = unique(s);
    rs = @(v) interp1(s_unique, v(unique_idx), sf, 'pchip');

    X_HOT = rs(X_HOT);
    R_HOT = rs(R_HOT);
    X_FLOOR = rs(X_FLOOR);
    R_FLOOR = rs(R_FLOOR);
    R_MID = rs(R_MID);
    X_ROOF = rs(X_ROOF);
    R_ROOF = rs(R_ROOF);
    THETA = rs(THETA);
    Z = zeros(NPTS, 1);

    writematrix([X_HOT, R_HOT, Z], "CAD_curves/nozzle_hot_side.csv");
    writematrix([X_FLOOR, R_FLOOR, Z], "CAD_curves/nozzle_channel_floor.csv");
    writematrix([X_ROOF, R_ROOF, Z], "CAD_curves/nozzle_channel_roof.csv");
    
    % Helix for sweep path
    writematrix([X_HOT, 2 * R_MID(1).*cos(THETA), 2 * R_MID(1).*sin(THETA)], "CAD_curves/nozzle_helix.csv");
    writematrix([X_HOT, Z, Z], "CAD_curves/nozzle_centerline.csv")
    
    % Cross section for sweep
    W = min(Ws);
    SECTION = [
        X_FLOOR(1), 2 * R_MID(1), -W/2;
        X_FLOOR(1), 0, -W/2;
        X_FLOOR(1), 0, W/2;
        X_FLOOR(1), 2 * R_MID(1), W/2;
        X_FLOOR(1), 2 * R_MID(1), -W/2;
        ];

    writematrix(SECTION, "CAD_curves/nozzle_channel_section.csv");

    fprintf('\nGeometry csv files written!\n');

end