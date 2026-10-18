function plotThrottleCurve(Engine, GasData)

    NominalEngine = createEngine(Engine, GasData);

    NominalEngine.LOCK_GEOMETRY = true;
    AREA_THROAT = NominalEngine.AREA_THROAT;
    AREA_EXIT = NominalEngine.AREA_EXIT;

    [H_REACTANTS, GasData.B_I_CIRC] = chemistrySolve(Engine.OF_RATIO, GasData);

    % Throttle sweep: 100% down to 20% chamber pressure
    THROTTLE_PCT = linspace(0.2, 1, 20);
    throttles = numel(THROTTLE_PCT);

    P_AMBIENT = linspace(0, 101325, 20);
    pambients = numel(P_AMBIENT);

    THROTTLE_THRUST = zeros(pambients,throttles);
    THROTTLE_MDOT = zeros(pambients,throttles);
    THROTTLE_ISP = zeros(pambients,throttles);
    SEPARATED = false(pambients,throttles);
    
    

    for i = 1:pambients

        for j = 1:throttles

            fprintf('\nSolving engine at Pa = %.0f, Throttle = %.2f', ...
                    P_AMBIENT(i), THROTTLE_PCT(j));

            CurrentEngine = NominalEngine;
            
            CurrentEngine.P_AMBIENT = P_AMBIENT(i);
            CurrentEngine.P_CHAMBER = CurrentEngine.P_CHAMBER * THROTTLE_PCT(j);
    
            % Chamber + throat at the throttled pressure (per 1 kg/s)
            [T_CHAMBER, Nj_CHAMBER, PropellantChamber] = solveChamber(CurrentEngine.P_CHAMBER, GasData.NUM_SPECIES, H_REACTANTS, GasData);
            H0 = PropellantChamber.H;  S0_R = PropellantChamber.S_R;
            [~, ~, ~, ~, AREA_THROAT_PER_MDOT] = solveThroat(CurrentEngine.P_CHAMBER, H0, S0_R, T_CHAMBER, Nj_CHAMBER, GasData, 1.0);
    
            % Actual mdot from the fixed physical throat area
            MDOT = AREA_THROAT / AREA_THROAT_PER_MDOT;
            
            CurrentEngine.MDOT = MDOT;
            CurrentEngine = solveRealEngine(CurrentEngine, GasData);
    
            THROTTLE_MDOT(i, j) = CurrentEngine.MDOT;
            THROTTLE_THRUST(i, j) = CurrentEngine.THRUST;
            THROTTLE_ISP(i, j) = CurrentEngine.THRUST / (CurrentEngine.MDOT * 9.80665);
    
            CurrentEngine = applySeparation(CurrentEngine);
            SEPARATED(i, j) = CurrentEngine.SEPARATES;
            
        end
    end
    
    [P_Ambient_Mesh, Thrust_Mesh] = ndgrid(P_AMBIENT, THROTTLE_PCT);
    
    Thrust_Map = griddedInterpolant(P_Ambient_Mesh, Thrust_Mesh, THROTTLE_THRUST, 'spline');
    Mdot_Map = griddedInterpolant(P_Ambient_Mesh, Thrust_Mesh, THROTTLE_MDOT, 'spline');
    Isp_Map = griddedInterpolant(P_Ambient_Mesh, Thrust_Mesh, THROTTLE_ISP, 'spline');
    
    SEPARATED = double(SEPARATED);
    Separation_Map = griddedInterpolant(P_Ambient_Mesh, Thrust_Mesh, SEPARATED, 'nearest');

    THROTTLE_PCT = THROTTLE_PCT * 100; % Convert to percentage for the plot axes
    [Pa_PlotMesh, Th_PlotMesh] = ndgrid(P_AMBIENT, THROTTLE_PCT); % Scaled mesh for plotting
    
    % 2. Calculate the 3D Separation Boundary Line
    % Scans each ambient pressure to find the exact throttle where separation begins
    sep_Throttle = NaN(pambients, 1);
    sep_Thrust = NaN(pambients, 1);
    sep_Mdot = NaN(pambients, 1);
    
    for k = 1:pambients
        idx = find(SEPARATED(k, :), 1, 'first');
        if ~isempty(idx)
            sep_Throttle(k) = THROTTLE_PCT(idx);
            sep_Thrust(k) = THROTTLE_THRUST(k, idx) / 1000;
            sep_Mdot(k) = THROTTLE_MDOT(k, idx);
        end
    end

    % 3. Generate the 3D Surfaces
    figure('Units','centimeters','Position',[5 5 30 14], 'Color','w');
    
    % --- LEFT PANEL: THRUST SURFACE ---
    subplot(1, 2, 1);
    surf(Th_PlotMesh, Pa_PlotMesh, THROTTLE_THRUST / 1000, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    hold on;
    % Overlay the red separation boundary line
    sep_bound_thrust = plot3(sep_Throttle, P_AMBIENT, sep_Thrust, 'r--', 'LineWidth', 3, 'DisplayName', 'Separation Boundary');
    
    xlabel('Throttle (%)', 'FontWeight','bold');
    ylabel('Ambient Pressure (Pa)', 'FontWeight','bold');
    zlabel('Thrust (kN)', 'FontWeight','bold');
    title(sprintf('Thrust Map (fixed A_t = %.2f cm^2, A_e = %.2f cm^2)', AREA_THROAT*1e4, AREA_EXIT*1e4));
    
    grid on; box on; set(gca, 'YDir', 'reverse'); % Reversing Y puts sea level in front
    colormap parula; view([-45 30]);
    legend(sep_bound_thrust, 'Separation boundary', 'Location', 'best');
    
    % --- RIGHT PANEL: ISP SURFACE ---
    subplot(1, 2, 2);
    surf(Th_PlotMesh, Pa_PlotMesh, THROTTLE_MDOT, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    hold on;
    % Overlay the red separation boundary line
    sep_bound_mdot = plot3(sep_Throttle, P_AMBIENT, sep_Mdot, 'r--', 'LineWidth', 3, 'DisplayName', 'Separation Boundary');
    
    xlabel('Throttle (%)', 'FontWeight','bold');
    ylabel('Ambient Pressure (Pa)', 'FontWeight','bold');
    zlabel('Specific Impulse (s)', 'FontWeight','bold');
    title('Specific Impulse Map');
    
    grid on; box on; set(gca, 'YDir', 'reverse');
    colormap parula; view([-45 30]);
    legend(sep_bound_mdot, 'Separation boundary', 'Location', 'best');
    
end
