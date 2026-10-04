function Engine = coolantPressureDrop(Engine)

    R = Engine.STREAM_R;
    X = Engine.STREAM_X;
    W = Engine.CHANNEL_WIDTH_STATION;
    H = Engine.CHANNEL_HEIGHT_STATION;
    N = Engine.NUM_CHANNELS;
    HELIX = Engine.HELIX_FACTOR;
    RHO = Engine.RHO_COOLANT;
    MU  = Engine.MU_COOLANT;
    EPSILON = Engine.CHANNEL_ROUGHNESS;

    % Axial slant length per segment, lengthened by the helix factor
    ds_AXIAL = hypot(diff(X), diff(R));
    ds = ds_AXIAL .* 0.5.*(HELIX(1:end-1) + HELIX(2:end)); % True helical segment length

    % Per-station
    AREA_CHANNEL = W .* H;
    DIAMETER_HYDRAULIC  = 2 .* W .* H ./ (W + H); % m
    V = Engine.MDOT_COOLANT ./ (RHO .* AREA_CHANNEL .* N); % m/s
    Re = RHO .* V .* DIAMETER_HYDRAULIC ./ MU;

    % Serghides's solution to the Colebrook–White equation, surface roughness
    A = -2 * log10(EPSILON ./ (DIAMETER_HYDRAULIC * 3.7) + 12 ./ Re);
    B = -2 * log10(EPSILON ./ (DIAMETER_HYDRAULIC * 3.7) + 2.51 .* A ./ Re);
    C = -2 * log10(EPSILON ./ (DIAMETER_HYDRAULIC * 3.7) + 2.51 .* B ./ Re);
    f_inv_sqrt = A - (B - A).^2 ./ (C - 2 .* B + A);
    f = 1 ./ (f_inv_sqrt).^2;

    % Darcy-Weisbach dP per segment
    f_SEG = 0.5*(f(1:end-1) + f(2:end));
    DIAMETER_HYDRAULIC_SEG = 0.5 * (DIAMETER_HYDRAULIC(1:end-1) + DIAMETER_HYDRAULIC(2:end));
    V_SEG = 0.5*(V(1:end-1) + V(2:end));
    dP_SEG = f_SEG .* (ds ./ DIAMETER_HYDRAULIC_SEG) .* 0.5 .* RHO .* V_SEG.^2;

    Engine.COOLANT_DP = sum(dP_SEG); % Pa
    Engine.CHANNEL_LENGTH = sum(ds); % m
    Engine.COOLANT_V_MAX  = max(V); % m/s

end