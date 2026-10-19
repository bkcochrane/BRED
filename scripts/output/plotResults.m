function plotResults(Engine)

    [xRaw, iu] = unique(Engine.Stream.x);
    rRaw = Engine.Stream.r(iu);
    machRaw = Engine.Stream.mach(iu);
    pRaw = Engine.Stream.p(iu);
    tempRaw = Engine.Stream.temp(iu);
    tempWGBartz = Engine.Wall.Bartz.tempWG(iu);
    tempWIBartz = Engine.Wall.Bartz.tempWI(iu);
    tempCoolantBartz = Engine.Wall.Bartz.tempCoolant(iu);
    tempWGIevlev = Engine.Wall.Ievlev.tempWG(iu);
    tempWIIevlev = Engine.Wall.Ievlev.tempWI(iu);
    tempCoolantIevlev = Engine.Wall.Ievlev.tempCoolant(iu);

    xFine = linspace(min(xRaw), max(xRaw), 100).';
    mach = interp1(xRaw,machRaw,xFine,'pchip');
    p = interp1(xRaw,pRaw,xFine,'pchip'); % Pa
    temp = interp1(xRaw,tempRaw,xFine,'pchip'); % K
    pPreSep = Engine.Stream.pAttached(iu) / Engine.Chamber.P; % Pre separation
    pPostSep = Engine.Stream.pPostSep(iu) / Engine.Chamber.P; % Post separation
    xFine = xFine * 1e3; XR = xRaw * 1e3; RR = rRaw * 1e3; % mm
    xMin = min(XR); x_max=max(XR); r_max = max(RR); XLIMS=[xMin x_max];

    figW=18; figH=24;
    figEngineResults = figure('Units','centimeters','Position',[2 2 figW figH],'Color','w'); clf;

    colors = [0, 0.5434, 0.9833;
              0.6068, 0.1418, 0.6805;
              1.0, 0, 0.3180];

    subText = sprintf(['Thrust = %.1f kN, O/F = %.1f, P_c = %.0f bar\n', ...
        'Area Ratio = %.1f, ISP = %.0f s, C* = %.0f m/s'], ...
        Engine.Specs.TARGET_THRUST / 1000, Engine.Specs.ofRatio, Engine.Chamber.P / 1e5, ...
        Engine.Specs.areaRatio, Engine.Specs.Isp, Engine.Specs.cStar);
    annotation('textbox',[0.1 0.945 0.8 0.05],'String', ...
        {'\bf\fontsize{14}Nozzle Design',['\rm\fontsize{9}' subText]}, ...
        'EdgeColor','none','HorizontalAlignment','center');

    left=0.12; top=0.88; bottom=0.07; gap=0.06;
    width=0.7;
    hContMax=0.24;

    panelWcm = figW * width;
    panelHcm = figH * hContMax;
    dataW = x_max - xMin;
    dataH = 2 * r_max * 1.1;
    panelAspect = panelHcm / panelWcm;
    dataAspect  = dataH / dataW;

    if dataAspect > panelAspect
        hCont  = hContMax;
        widthC = (panelHcm / dataAspect) / figW;
    else
        widthC = width;
        hCont  = (panelWcm * dataAspect) / figH;
    end

    % ax1: Contour
    left_c = left + (width - widthC)/2;
    ax1 = axes('Position', [left_c, top-hCont, widthC, hCont]);
    plot(ax1, XR, RR, 'k','LineWidth',1.5); hold(ax1,'on');
    plot(ax1, XR, -RR, 'k','LineWidth',1.5);
    xlim(ax1, XLIMS); ylim(ax1, [-r_max r_max]*1.1);
    ax1.DataAspectRatio = [1 1 1];
    ax1.PositionConstraint = 'innerposition';
    ylabel(ax1,'r (mm)'); grid(ax1,'on'); box(ax1,'on'); ax1.XTickLabel=[];
    title(ax1,'Nozzle contour');

    yCont = top - hCont;
    avail = (yCont - gap) - bottom;
    hPanel = (avail - gap)/2;
    yFlow = yCont - gap - hPanel;
    yWall = bottom;
    panelLeft = left_c;
    panelWidth = widthC; 

    % ax2: Mach / P ---
    ax2 = axes('Position',[panelLeft, yFlow, panelWidth, hPanel]);
    hold on
    yyaxis(ax2,'left');
    hM=plot(ax2,xFine,mach,'LineWidth',1.5,'Color',colors(1,:));
    ylabel(ax2,'Mach','Color',colors(1,:)); ax2.YColor=colors(1,:);
    yline(ax2,1,':','M = 1','Color',colors(1,:));
    yyaxis(ax2,'right');
    hPP = plot(ax2, XR, pPreSep, '-', 'LineWidth',1.5, 'Color',colors(2,:));
    plot(ax2, XR, pPostSep, '-', 'LineWidth',1.5, 'Color',colors(2,:));
    if Engine.Stream.separates
        hP = plot(ax2, xFine, p / Engine.Chamber.P, '--', 'LineWidth',1.0, 'Color',colors(2,:));
        xs = Engine.Stream.xSep*1e3;
        [~, is] = min(abs(xRaw - Engine.Stream.xSep));
        ylo = Engine.Stream.pAttached(is) / Engine.Chamber.P;
        yhi = Engine.Stream.pPostSep(is) / Engine.Chamber.P;
        plot(ax2, [xs xs], [ylo yhi], ':', 'Color', colors(2,:), 'LineWidth',1.2);
        
    end
    ylabel(ax2,'P/P_C','Color',colors(2,:)); ax2.YColor=colors(2,:);
    grid(ax2,'on'); box(ax2,'on'); ax2.XTickLabel=[]; xlim(ax2,XLIMS);
    title(ax2,'Flow properties');
    
    % --- ax3: T
    d = 0.09;
    ax3 = axes('Position',[panelLeft, yFlow, panelWidth + d, hPanel]);
    hT=plot(ax3,xFine,temp,'LineWidth',1.5,'Color',colors(3,:));
    ax3.Color='none'; ax3.YAxisLocation='right'; ax3.XColor='none';
    ax3.YColor=colors(3,:); box(ax3,'off'); ax3.XTick=[];
    ylabel(ax3,'T (K)','Color',colors(3,:));
    xlim(ax3,XLIMS); ax3.YLim=[0 max(temp)*1.1];
    legend_items = [hM hPP hT];
    legend_entry = {'Mach Number', 'Isentropic Wall Pressure', 'Temperature'};
    if Engine.Stream.separates
        legend_items = [hM hPP hP hT];
        legend_entry = {'Mach Number','Wall Pressure', 'Isentropic Wall Pressure', 'Temperature'};
    end
    legend(ax3, legend_items, legend_entry,'Location','west');
    x_span = XLIMS(2) - XLIMS(1);
    xmax_ext = XLIMS(1) + x_span * (panelWidth + d) / panelWidth;
    xlim(ax3, [XLIMS(1), xmax_ext]); ylim([min(temp) * 0.995, max(temp) * 1.005])

    % ax4: Wall temps ---
    ax4 = axes('Position',[panelLeft, yWall, panelWidth, hPanel]);
    p1 = plot(ax4, XR, tempWGBartz, 'LineWidth', 1.5, 'Color', colors(1,:));
    hold(ax4, 'on');
    p2 = plot(ax4, XR, tempWIBartz, 'LineWidth', 1.5, 'Color', colors(2,:));
    p3 = plot(ax4, XR, tempCoolantBartz, 'LineWidth', 1.5, 'Color', colors(3,:));
    
    plot(ax4, XR, tempWGIevlev, '--', 'LineWidth', 1.5, 'Color', colors(1,:));
    plot(ax4, XR, tempWIIevlev, '--', 'LineWidth', 1.5, 'Color', colors(2,:));
    plot(ax4, XR, tempCoolantIevlev, '--', 'LineWidth', 1.5, 'Color', colors(3,:));
    
    ylabel(ax4,'T (K)'); xlabel(ax4,'x (mm)');
    grid(ax4,'on'); box(ax4,'on');
    ylim(ax4,[0 max(max(tempWGBartz), max(tempWGIevlev))*1.1]); 
    xlim(ax4,XLIMS);
    title(ax4,'Wall temperatures');

    legend(ax4, [p1, p2, p3], { ...
        'Gas-side wall', ...
        'Coolant-side wall', ...
        'Coolant' }, ...
        'Location', 'southwest');

    axOverlay = axes('Position', get(ax4, 'Position'), ...
               'Color', 'none', ...
               'XAxisLocation', 'top', ...
               'YAxisLocation', 'right', ...
               'XColor', 'none', ...
               'YColor', 'none');
             
    hold(axOverlay, 'on');
    xlim(axOverlay, xlim(ax4)); 
    ylim(axOverlay, ylim(ax4));

    d1 = plot(axOverlay, NaN, NaN, '-k', 'LineWidth', 1.5);
    d2 = plot(axOverlay, NaN, NaN, '--k', 'LineWidth', 1.5);

    legend(axOverlay, [d1, d2], { ...
        'Bartz', ...
        'Ievlev' }, ...
        'Location', 'northwest');
        
    set(axOverlay, 'HitTest', 'off', 'PickableParts', 'none');

    % ---

    axs = findall(gcf, 'Type', 'axes');
    set(axs, 'FontSize', 10);
    for a = axs'
        a.Title.FontSize = 12;
    end
    
    if ~exist('results', 'dir')
        mkdir('results');
    end
    
    print(figEngineResults, 'results/engine_results.pdf', '-dpdf');

end