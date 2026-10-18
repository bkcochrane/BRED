function Engine = coolantPressureDrop(Engine)

    r = Engine.Stream.r;
    x = Engine.Stream.x;
    width = Engine.Cooling.widthChannelStation;
    height = Engine.Cooling.heightChannelStation;
    numChannels = Engine.Cooling.numChannels;
    RHO = Engine.Cooling.Coolant.RHO;
    MU = Engine.Cooling.Coolant.MU;
    EPSILON = Engine.Cooling.CHANNEL_ROUGHNESS;

    ds = hypot(diff(x), diff(r));

    % Per-station
    areaChannel = width .* height;
    diameterHydraulic  = 2 .* width .* height ./ (width + height); % m
    v = Engine.Cooling.Coolant.mdot ./ (RHO .* areaChannel .* numChannels); % m/s
    reynolds = RHO .* v .* diameterHydraulic ./ MU;

    % Serghides's solution to the Colebrook–White equation, surface roughness
    a = -2 * log10(EPSILON ./ (diameterHydraulic * 3.7) + 12 ./ reynolds);
    b = -2 * log10(EPSILON ./ (diameterHydraulic * 3.7) + 2.51 .* a ./ reynolds);
    c = -2 * log10(EPSILON ./ (diameterHydraulic * 3.7) + 2.51 .* b ./ reynolds);
    fInvSqrt = a - (b - a).^2 ./ (c - 2 .* b + a);
    f = 1 ./ (fInvSqrt).^2;

    % Darcy-Weisbach dP per segment
    fSegment = 0.5*(f(1:end-1) + f(2:end));
    diameterHydraulicSegment = 0.5 * (diameterHydraulic(1:end-1) + diameterHydraulic(2:end));
    vSegment = 0.5*(v(1:end-1) + v(2:end));
    dPSegment = fSegment .* (ds ./ diameterHydraulicSegment) .* 0.5 .* RHO .* vSegment.^2;

    Engine.Cooling.Coolant.dP = sum(dPSegment); % Pa
    Engine.Cooling.lengthChannel = sum(ds); % m
    Engine.Cooling.Coolant.vMax = max(v); % m/s

end