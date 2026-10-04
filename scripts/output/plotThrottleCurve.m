function plotThrottleCurve(Engine, GAS_DATA)

    nominal = createEngine(Engine, GAS_DATA);
    AREA_THROAT = nominal.AREA_THROAT;
    AREA_EXIT = nominal.AREA_EXIT;
    M_EXIT = nominal.STREAM_M(end); % Exit Mach constant at fixed Area Ratio

    [H_REACTANTS, GAS_DATA.B_I_CIRC] = chemistrySolve(Engine.OF_RATIO);

    % Throttle sweep: 100% down to 20% chamber pressure
    THROTTLE_PCT = linspace(1.0, 0.20, 40);
    N = numel(THROTTLE_PCT);

    THROTTLE_THRUST = zeros(1,N);
    THROTTLE_ISP    = zeros(1,N);
    THROTTLE_MDOT   = zeros(1,N);
    SEPARATED = false(1,N);

    for i = 1:N
        Pc = nominal.P_CHAMBER * THROTTLE_PCT(i);

        % Chamber + throat at the throttled pressure (per 1 kg/s)
        [Tc, Nj_c, prop_c] = solveChamber(Pc, GAS_DATA.NUM_SPECIES, H_REACTANTS, GAS_DATA);
        H0 = prop_c.H;  S0_R = prop_c.S_R;
        [Pt, Tt, Njt, ~, At_per_mdot] = solveThroat(Pc, H0, S0_R, Tc, Nj_c, GAS_DATA, 1.0);

        % Actual mdot from the FIXED physical throat area
        MDOT = AREA_THROAT / At_per_mdot;

        % Exit pressure at the fixed physical area ratio
        AREA_EXIT_PER_MDOT = AREA_EXIT / MDOT;
        local = nominal; local.P_THROAT = Pt; local.P_CHAMBER = Pc;
        P_EXIT = solvePressure(S0_R, H0, local, GAS_DATA, Tt, Njt, 1.0, AREA_EXIT_PER_MDOT, Pt, 'supersonic');

        [~,~,prop_e] = findEquilibrium('SP', S0_R, P_EXIT, GAS_DATA, Tt, Njt, 0);
        V_EXIT = sqrt(max(2*(H0 - prop_e.H), 0)) * Engine.ETA_NOZZLE;

        THROTTLE_MDOT(i) = MDOT;
        THROTTLE_THRUST(i) = MDOT * V_EXIT + (P_EXIT - Engine.P_AMBIENT) * AREA_EXIT;
        THROTTLE_ISP(i) = THROTTLE_THRUST(i) / (MDOT * 9.80665);

        % Schmucker separation check
        P_sep = (1.88*M_EXIT - 1)^(-0.64) * Engine.P_AMBIENT;
        SEPARATED(i) = P_EXIT < P_sep;
    end

    MDOT_PCT = (THROTTLE_MDOT / nominal.MDOT) * 100;
    SEP_idx = find(SEPARATED, 1, 'first');

    figure('Units','centimeters','Position',[5 5 20 12], 'Color','w');
    neb = nebula(3);

    yyaxis left
    plot(MDOT_PCT, THROTTLE_THRUST/1000, '-', 'LineWidth',2.5, 'Color',neb(1,:));
    ylabel('Thrust (kN)', 'Color',neb(1,:), 'FontWeight','bold');
    ax = gca; ax.YColor = neb(1,:);

    yyaxis right
    hold on
    if isempty(SEP_idx)
        plot(MDOT_PCT, THROTTLE_ISP, '-', 'LineWidth',2.5, 'Color',neb(2,:));
    else
        plot(MDOT_PCT(1:SEP_idx),     THROTTLE_ISP(1:SEP_idx),     '-',  'LineWidth',2.5, 'Color',neb(2,:));
        plot(MDOT_PCT(SEP_idx:end),   THROTTLE_ISP(SEP_idx:end),   '--', 'LineWidth',2.5, 'Color',neb(2,:));
        xline(MDOT_PCT(SEP_idx), ':', 'Color','r', 'LineWidth',2, ...
            'Label','Separation', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment', 'bottom');
    end
    ylabel('Specific Impulse (s)', 'Color',neb(2,:), 'FontWeight','bold');
    ax.YColor = neb(2,:);

    xlabel('Throttle Setting (% of Nominal Mass Flow)', 'FontWeight','bold');
    ax.XColor = 'k';
    title(sprintf('Throttle Curve (fixed A_t = %.2f cm^2, A_e = %.2f cm^2)', ...
        AREA_THROAT*1e4, AREA_EXIT*1e4));
    xlim([min(MDOT_PCT) 100]); grid on; box on;
    set(gcf, 'WindowStyle', 'docked');
    
end
