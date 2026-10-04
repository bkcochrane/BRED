function plotResults(Engine)

    [X_raw, iu] = unique(Engine.STREAM_X);
    R_raw=Engine.STREAM_R(iu); M_raw=Engine.STREAM_M(iu);
    P_raw=Engine.STREAM_P(iu); T_raw=Engine.STREAM_T(iu);
    T_WG=Engine.WALL_TWG(iu); T_WI=Engine.WALL_TWI(iu); T_COOL=Engine.WALL_TCOOLANT(iu);

    X_fine = linspace(min(X_raw), max(X_raw), 100).';
    M = interp1(X_raw,M_raw,X_fine,'pchip');
    P = interp1(X_raw,P_raw,X_fine,'pchip'); % Pa
    T = interp1(X_raw,T_raw,X_fine,'pchip'); % K
    P_PRE = Engine.STREAM_P_PRE(iu) / Engine.P_CHAMBER; % Pre separation
    P_POST = Engine.STREAM_P_POST(iu) / Engine.P_CHAMBER; % Post separation
    XF = X_fine * 1e3; XR = X_raw * 1e3; RR = R_raw * 1e3; % mm
    x_min = min(XR); x_max=max(XR); r_max = max(RR); XLIMS=[x_min x_max];

    fig_w=18; fig_h=24;
    fig_engine_results = figure('Units','centimeters','Position',[2 2 fig_w fig_h],'Color','w'); clf;
    colors = nebula(3);

    subText = sprintf(['Thrust = %.1f kN, O/F = %.1f, P_c = %.0f bar\n', ...
        'Area Ratio = %.1f, ISP = %.0f s, C* = %.0f m/s'], ...
        Engine.TARGET_THRUST / 1000, Engine.OF_RATIO, Engine.P_CHAMBER / 1e5, ...
        Engine.AREA_RATIO, Engine.ISP, Engine.C_STAR);
    annotation('textbox',[0.1 0.945 0.8 0.05],'String', ...
        {'\bf\fontsize{14}Nozzle Design',['\rm\fontsize{9}' subText]}, ...
        'EdgeColor','none','HorizontalAlignment','center');

    LEFT=0.12; TOP=0.88; BOTTOM=0.07; GAP=0.06;
    WIDTH=0.7;
    H_CONT_MAX=0.24;

    panel_w_cm = fig_w * WIDTH;
    panel_h_cm = fig_h * H_CONT_MAX;
    data_w = x_max - x_min;
    data_h = 2 * r_max * 1.1;
    panel_aspect = panel_h_cm / panel_w_cm;
    data_aspect  = data_h / data_w;

    if data_aspect > panel_aspect
        H_CONT  = H_CONT_MAX;
        WIDTH_C = (panel_h_cm / data_aspect) / fig_w;
    else
        WIDTH_C = WIDTH;
        H_CONT  = (panel_w_cm * data_aspect) / fig_h;
    end

    % ax1: Contour
    left_c = LEFT + (WIDTH - WIDTH_C)/2;
    ax1 = axes('Position', [left_c, TOP-H_CONT, WIDTH_C, H_CONT]);
    plot(ax1, XR, RR, 'k','LineWidth',1.5); hold(ax1,'on');
    plot(ax1, XR, -RR, 'k','LineWidth',1.5);
    xlim(ax1, XLIMS); ylim(ax1, [-r_max r_max]*1.1);
    ax1.DataAspectRatio = [1 1 1];
    ax1.PositionConstraint = 'innerposition';
    ylabel(ax1,'r (mm)'); grid(ax1,'on'); box(ax1,'on'); ax1.XTickLabel=[];
    title(ax1,'Nozzle contour');

    Y_CONT  = TOP - H_CONT;
    avail   = (Y_CONT - GAP) - BOTTOM;
    H_PANEL = (avail - GAP)/2;
    y_flow  = Y_CONT - GAP - H_PANEL;
    y_wall  = BOTTOM;
    PANEL_LEFT  = left_c;
    PANEL_WIDTH = WIDTH_C; 

    % ax2: Mach / P ---
    ax2 = axes('Position',[PANEL_LEFT, y_flow, PANEL_WIDTH, H_PANEL]);
    hold on
    yyaxis(ax2,'left');
    hM=plot(ax2,XF,M,'LineWidth',1.5,'Color',colors(1,:));
    ylabel(ax2,'Mach','Color',colors(1,:)); ax2.YColor=colors(1,:);
    yline(ax2,1,':','M = 1','Color',colors(1,:));
    yyaxis(ax2,'right');
    hPP = plot(ax2, XR, P_PRE, '-', 'LineWidth',1.5, 'Color',colors(2,:));
    plot(ax2, XR, P_POST, '-', 'LineWidth',1.5, 'Color',colors(2,:));
    if Engine.SEPARATES
        hP = plot(ax2, XF, P / Engine.P_CHAMBER, '--', 'LineWidth',1.0, 'Color',colors(2,:));
        xs = Engine.X_SEP*1e3;
        [~, is] = min(abs(X_raw - Engine.X_SEP));
        ylo = Engine.STREAM_P_PRE(is) / Engine.P_CHAMBER;
        yhi = Engine.STREAM_P_POST(is) / Engine.P_CHAMBER;
        plot(ax2, [xs xs], [ylo yhi], ':', 'Color', colors(2,:), 'LineWidth',1.2);
        
    end
    ylabel(ax2,'P/P_C','Color',colors(2,:)); ax2.YColor=colors(2,:);
    grid(ax2,'on'); box(ax2,'on'); ax2.XTickLabel=[]; xlim(ax2,XLIMS);
    title(ax2,'Flow properties');
    
    % --- ax3: T
    d = 0.09;
    ax3 = axes('Position',[PANEL_LEFT, y_flow, PANEL_WIDTH + d, H_PANEL]);
    hT=plot(ax3,XF,T,'LineWidth',1.5,'Color',colors(3,:));
    ax3.Color='none'; ax3.YAxisLocation='right'; ax3.XColor='none';
    ax3.YColor=colors(3,:); box(ax3,'off'); ax3.XTick=[];
    ylabel(ax3,'T (K)','Color',colors(3,:));
    xlim(ax3,XLIMS); ax3.YLim=[0 max(T)*1.1];
    legend_items = [hM hPP hT];
    legend_entry = {'Mach Number', 'Isentropic Wall Pressure', 'Temperature'};
    if Engine.SEPARATES
        legend_items = [hM hPP hP hT];
        legend_entry = {'Mach Number','Wall Pressure', 'Isentropic Wall Pressure', 'Temperature'};
    end
    legend(ax3, legend_items, legend_entry,'Location','west');
    x_span = XLIMS(2) - XLIMS(1);
    xmax_ext = XLIMS(1) + x_span * (PANEL_WIDTH + d) / PANEL_WIDTH;
    xlim(ax3, [XLIMS(1), xmax_ext]); ylim([min(T) * 0.995, max(T) * 1.005])

    % ax4: Wall temps ---
    ax4 = axes('Position',[PANEL_LEFT, y_wall, PANEL_WIDTH, H_PANEL]);
    hold(ax4,'on');
    plot(ax4,XR,T_WG,'LineWidth',1.5,'Color',colors(1,:));
    plot(ax4,XR,T_WI,'LineWidth',1.5,'Color',colors(2,:));
    plot(ax4,XR,T_COOL,'LineWidth',1.5,'Color',colors(3,:));
    ylabel(ax4,'T (K)'); xlabel(ax4,'x (mm)');
    grid(ax4,'on'); box(ax4,'on');
    ylim(ax4,[0 max(T_WG)*1.1]); xlim(ax4,XLIMS);
    legend(ax4,{'Gas-side wall','Coolant-side wall','Coolant'},'Location','southwest');
    title(ax4,'Wall temperatures');

    axs = findall(gcf, 'Type', 'axes');
    set(axs, 'FontSize', 10);
    for a = axs'
        a.Title.FontSize = 12;
    end
    
    print(fig_engine_results, 'results/engine_results.pdf', '-dpdf');

    fig_engine_results.WindowStyle = 'docked';

end