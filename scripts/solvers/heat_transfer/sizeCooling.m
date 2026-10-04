function Engine = sizeCooling(Engine)

    R_STATIONS = Engine.STREAM_R;
    WIDTH_MIN  = Engine.CHANNEL_WIDTH_MIN; % m
    LAND_WIDTH = Engine.LAND_WIDTH; % m
    R_THROAT   = min(R_STATIONS); % m
    R_MAX      = max(R_STATIONS); % m

    SPIRAL = isfield(Engine,'N_STARTS') && Engine.N_STARTS > 0;

    if SPIRAL
        NUM_CHANNELS = Engine.N_STARTS;
        X_SPAN = max(Engine.STREAM_X) - min(Engine.STREAM_X);
        Engine.SPIRAL_RATE = (Engine.N_TURNS * 2 * pi) / X_SPAN; % rad/m
    else
        % Straight
        NUM_CHANNELS = floor(2 * pi * R_THROAT / (WIDTH_MIN + LAND_WIDTH));
        Engine.SPIRAL_RATE = 0;
    end

    Engine.NUM_CHANNELS = NUM_CHANNELS;

    HEIGHT_THROAT = Engine.CHANNEL_HEIGHT_THROAT;
    HEIGHT_MAX = Engine.CHANNEL_HEIGHT_MAX;
    HEIGHT_LOCAL = HEIGHT_THROAT + ((R_STATIONS - R_THROAT)/(R_MAX - R_THROAT)) .* (HEIGHT_MAX - HEIGHT_THROAT);
    
    WIDTH_LOCAL = WIDTH_MIN * ones(size(R_STATIONS));

    Engine.CHANNEL_WIDTH_STATION  = WIDTH_LOCAL;
    Engine.CHANNEL_HEIGHT_STATION = HEIGHT_LOCAL;

    % To account for the longer length of the helix
    Engine.HELIX_FACTOR = sqrt(1 + (R_STATIONS * Engine.SPIRAL_RATE).^2);

end