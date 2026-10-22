function h = coolantH(Engine, widthChannel, heightChannel, tempLocal)

    LAND_WIDTH = Engine.Cooling.LAND_WIDTH;
    numChannels = Engine.Cooling.numChannels;

    areaChannel = widthChannel * heightChannel;
    diameterHydraulic  = 2 * widthChannel * heightChannel / (widthChannel + heightChannel);

    Coolant = Engine.Cooling.Coolant;

    v = Coolant.mdot / (Coolant.RHO * areaChannel * numChannels);
    mu = Coolant.MU * (tempLocal / 300)^-1.5;
    reynolds = Coolant.RHO * v * diameterHydraulic / mu;
    prandtl = mu * Coolant.CP / Coolant.K;
    
    % Sieder-Tate
    nu = 0.027 * reynolds^0.8 * prandtl^(1/3);
    hBare = nu * Coolant.K / diameterHydraulic;

    % Fin efficiency of the lands
    m = sqrt(2 * hBare / (Engine.Wall.K * LAND_WIDTH));
    eta = tanh(m * heightChannel) / (m * heightChannel);
    areaEffectiveRatio = (widthChannel + 2 * eta * heightChannel) / (widthChannel + LAND_WIDTH);

    h = hBare * areaEffectiveRatio;

end