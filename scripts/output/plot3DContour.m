function plot3DContour(Engine)

    X = Engine.Geometry.X_COORDINATES(:) * 1e3; % mm
    R = Engine.Geometry.R_COORDINATES(:) * 1e3; % mm
    [X, iu] = unique(X);
    R = R(iu);
    
    % Thermal + channel geometry mapped onto the contour stations
    X_S = Engine.Stream.X * 1e3;
    T_WG = interp1(X_S, Engine.Wall.T_WG, X, 'linear', 'extrap');
    T_COOL = interp1(X_S, Engine.Wall.T_COOLANT, X, 'linear', 'extrap');
    WIDTH_STATION  = interp1(X_S, Engine.Cooling.CHANNEL_WIDTH_STATION  * 1e3, X, 'linear', 'extrap');
    HEIGHT_STATION = interp1(X_S, Engine.Cooling.CHANNEL_HEIGHT_STATION * 1e3, X, 'linear', 'extrap');
    
    WALL_THICKNESS = Engine.Wall.THICKNESS * 1e3;
    MAP_GAS  = hot();
    MAP_COOL = nebula();
    
    % Gas side wall
    NTHETA = 60;
    THETA = linspace(0, 2 * pi, NTHETA);
    [X_G, TH_G] = meshgrid(X, THETA);
    R_G = repmat(R.', NTHETA, 1);
    Y_G = R_G .* cos(TH_G);
    Z_G = R_G .* sin(TH_G);
    C_G = repmat(T_WG.', NTHETA, 1);
    TC_GAS = colormapToTruecolor(MAP_GAS, C_G);
    
    % Spiral cooling channels setup
    N_STARTS = Engine.Cooling.NUM_CHANNELS;
    if isfield(Engine.Cooling, 'SPRIAL_RATE')
        SPIRAL_RATE = Engine.Cooling.SPIRAL_RATE / 1e3;
    else
        SPIRAL_RATE = 0;
    end
    R_FLOOR = R + WALL_THICKNESS;
    R_ROOF  = R + WALL_THICKNESS + HEIGHT_STATION;
    THETA_STARTS = linspace(0, 2*pi, N_STARTS + 1);
    THETA_STARTS(end) = [];                           

    % --------------------------------------------------------------------
    % Wall and channels
    % --------------------------------------------------------------------
    fig_1 = figure('Units','centimeters','Position',[5 5 22 12], 'Color', 'w');
    tcl_1 = tiledlayout(fig_1, 1, 2, 'TileSpacing', 'compact');
    
    % Tile 1: Wall Only
    ax1 = nexttile(tcl_1);
    hold(ax1, 'on');
    
    surf(ax1, Y_G, Z_G, -X_G, TC_GAS, 'FaceColor','interp', 'EdgeColor','none');
    
    axis(ax1, 'equal'); axis(ax1, 'off');
    view(ax1, 35, 20);
    lighting(ax1, 'gouraud'); camlight(ax1, 'headlight'); material(ax1, 'dull');
    title(ax1, 'Gas-side wall', 'FontWeight', 'bold');

    % Tile 2: Channels Only
    ax2 = nexttile(tcl_1);
    hold(ax2, 'on');
    
    for c = 1:N_STARTS
        THETA_C = THETA_STARTS(c) + SPIRAL_RATE * (X - min(X));
        dTHETA = (WIDTH_STATION/2) ./ max(R_FLOOR, 1);
        TH_LO = THETA_C - dTHETA;
        TH_HI = THETA_C + dTHETA;
        
        surf3strip(ax2, X, R_ROOF,  TH_LO, TH_HI, T_COOL, MAP_COOL, true);
        surf3strip(ax2, X, R_FLOOR, TH_LO, TH_HI, T_COOL, MAP_COOL, true);
        surf3wall (ax2, X, R_FLOOR, R_ROOF, TH_LO, T_COOL, MAP_COOL, true);
        surf3wall (ax2, X, R_FLOOR, R_ROOF, TH_HI, T_COOL, MAP_COOL, true);
    end
    
    axis(ax2, 'equal'); axis(ax2, 'off');
    view(ax2, 35, 20);
    lighting(ax2, 'gouraud'); camlight(ax2, 'headlight'); material(ax2, 'dull');
    title(ax2, 'Cooling channels', 'FontWeight', 'bold');

    colormap(ax1, MAP_GAS);
    clim(ax1, [min(T_WG) max(T_WG)]);
    cb_gas = colorbar(ax1);
    cb_gas.Layout.Tile = 'south';
    cb_gas.Label.String = 'Gas-side wall T (K)';
    cb_gas.Label.FontWeight = 'normal';
    
    colormap(ax2, MAP_COOL);
    clim(ax2, [min(T_COOL) max(T_COOL)]);
    cb_cool = colorbar(ax2);
    cb_cool.Layout.Tile = 'south';
    cb_cool.Label.String = 'Coolant T (K)';
    cb_cool.Label.FontWeight = 'normal';

    set(gcf, 'WindowStyle', 'docked');
    
    % --------------------------------------------------------------------
    % Full assembly
    % --------------------------------------------------------------------
    fig_2 = figure('Units','centimeters','Position',[10 8 16 12], 'Color', 'w');
    tcl_2 = tiledlayout(fig_2, 1, 1, 'TileSpacing', 'compact');
    
    ax3 = nexttile(tcl_2);
    hold(ax3, 'on');
    
    surf(ax3, X_G, Y_G, Z_G, TC_GAS, 'FaceColor','interp', 'EdgeColor','none');
    
    for c = 1:N_STARTS
        THETA_C = THETA_STARTS(c) + SPIRAL_RATE * (X - min(X));
        dTHETA = (WIDTH_STATION/2) ./ max(R_FLOOR, 1);
        TH_LO = THETA_C - dTHETA;
        TH_HI = THETA_C + dTHETA;

        surf3strip(ax3, X, R_ROOF,  TH_LO, TH_HI, T_COOL, MAP_COOL, false);
        surf3strip(ax3, X, R_FLOOR, TH_LO, TH_HI, T_COOL, MAP_COOL, false);
        surf3wall (ax3, X, R_FLOOR, R_ROOF, TH_LO, T_COOL, MAP_COOL, false);
        surf3wall (ax3, X, R_FLOOR, R_ROOF, TH_HI, T_COOL, MAP_COOL, false);
    end
    
    axis(ax3, 'equal'); axis(ax3, 'off');
    view(ax3, 35, 10);
    lighting(ax3, 'gouraud'); camlight(ax3, 'headlight'); material(ax3, 'dull');
    title(ax3, 'Full Assembly', 'FontWeight', 'bold');

    ax_hidden = axes(tcl_2, 'Visible','off', 'HandleVisibility','off');
    colormap(ax3, MAP_GAS);
    colormap(ax_hidden, MAP_COOL);
    
    clim(ax3, [min(T_WG) max(T_WG)]);
    clim(ax_hidden, [min(T_COOL) max(T_COOL)]);
    
    cb_gas3 = colorbar(ax3);
    cb_gas3.Layout.Tile = 'south';
    cb_gas3.Label.String = 'Gas-side wall T (K)';
    cb_gas3.Label.FontWeight = 'normal';
    
    cb_cool3 = colorbar(ax_hidden);
    cb_cool3.Layout.Tile = 'south';
    cb_cool3.Label.String = 'Coolant T (K)';
    cb_cool3.Label.FontWeight = 'normal';

    set(gcf, 'WindowStyle', 'docked');

    % ------------------------------------------------------------------------
    % Helper: Swept strip at constant radius
    % ------------------------------------------------------------------------
    function surf3strip(AX, X, R_FACE, TH_LO, TH_HI, T, MAP, VERT)
        R_M  = [R_FACE, R_FACE];
        TH_M = [TH_LO, TH_HI];
        TC = colormapToTruecolor(MAP, [T, T]);
        
        X_cyl = [X, X];
        Y_cyl = R_M.*cos(TH_M);
        Z_cyl = R_M.*sin(TH_M);
        
        if VERT
            surf(AX, Y_cyl, Z_cyl, -X_cyl, TC, 'FaceColor','interp', 'EdgeColor','none');
        else
            surf(AX, X_cyl, Y_cyl, Z_cyl, TC, 'FaceColor','interp', 'EdgeColor','none');
        end
    end
    
    % ------------------------------------------------------------------------
    % Helper: Swept wall at constant theta
    % ------------------------------------------------------------------------
    function surf3wall(AX, X, R_LO, R_HI, TH, T, MAP, VERT)
        R_M  = [R_LO, R_HI];
        TH_M = [TH, TH];
        TC = colormapToTruecolor(MAP, [T, T]);
        
        X_cyl = [X, X];
        Y_cyl = R_M.*cos(TH_M);
        Z_cyl = R_M.*sin(TH_M);
        
        if VERT
            surf(AX, Y_cyl, Z_cyl, -X_cyl, TC, 'FaceColor','interp', 'EdgeColor','none');
        else
            surf(AX, X_cyl, Y_cyl, Z_cyl, TC, 'FaceColor','interp', 'EdgeColor','none');
        end
    end
    
    % ------------------------------------------------------------------------
    % Helper: Convert Z/C data to truecolor RGB via a colormap
    % ------------------------------------------------------------------------
    function TC = colormapToTruecolor(MAP, C_DATA)
        IDX = round(rescale(C_DATA, 1, height(MAP)));
        TC = reshape(MAP(IDX,:), [size(C_DATA), 3]);
    end

end