function Engine = sizeCooling(Engine)

    rStations = Engine.Stream.r;
    WIDTH_MIN = Engine.Cooling.CHANNEL_WIDTH_MIN; % m
    LAND_WIDTH = Engine.Cooling.LAND_WIDTH; % m
    rThoat = min(rStations); % m
    rMax = max(rStations); % m

    numChannels = floor(2 * pi * rThoat / (WIDTH_MIN + LAND_WIDTH));

    Engine.Cooling.numChannels = numChannels;

    HEIGHT_MIN = Engine.Cooling.CHANNEL_HEIGHT_MIN;
    HEIGHT_MAX = Engine.Cooling.CHANNEL_HEIGHT_MAX;
    heightLocal = HEIGHT_MIN + ((rStations - rThoat)/(rMax - rThoat)) .* (HEIGHT_MAX - HEIGHT_MIN);
    
    widthLocal = WIDTH_MIN * ones(size(rStations));

    Engine.Cooling.widthChannelStation  = widthLocal;
    Engine.Cooling.heightChannelStation = heightLocal;

end