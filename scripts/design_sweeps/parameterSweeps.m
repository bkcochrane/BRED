function parameterSweeps(Engine, GAS_DATA)

    fig_engine_sweeps = figure('Units','centimeters','Position',[2 2 21 29.7], ...
       'PaperType','A4','PaperOrientation','portrait', ...
       'PaperUnits','centimeters','PaperPosition',[0 0 21 29.7]);
    
    BASE_THRUST=1000; BASE_OF=2.4; BASE_PC=20e5; BASE_PR=0.5; BASE_LSTAR=1.3; BASE_CONE=90; BASE_CR=6;
    
    sweeps(1) = struct('field','TARGET_THRUST', 'values',400:200:2000, 'label','Thrust (N)', 'scale',1, 'title','Varying thrust');
    sweeps(2) = struct('field','P_EXIT', 'values',(0.4:0.1:1), 'label','Exit pressure ratio', 'scale',1, 'title','Varying exit pressure ratio');
    Engine.EXIT_MODE = 'pressure';
    sweeps(3) = struct('field','P_CHAMBER', 'values',5:5:40, 'label','Chamber pressure (bar)','scale',1e5, 'title','Varying chamber pressure');
    sweeps(4) = struct('field','L_STAR', 'values',0.8:0.2:1.8, 'label','L* (m)', 'scale',1, 'title','Varying L*');
    sweeps(5) = struct('field','CONTRACTION_RATIO','values',5:1:9, 'label','Contraction ratio', 'scale',1, 'title','Varying contraction ratio');
    % sweeps(6) = struct('field','OF_RATIO', 'values',1.8:0.2:3.2, 'label','O/F ratio', 'scale',1, 'title','Varying O/F ratio');
    % sweeps(7) = struct('field','CONE_PERCENT', 'values',[60 80 90], 'label','Bell (%)', 'scale',1, 'title','Varying bell percent');
    
    tl = tiledlayout(numel(sweeps),1,'TileSpacing','compact','Padding','loose');

    xmax = 0; rmax = 0;
    xmin = 0;
    data = cell(numel(sweeps),1);
    for s = 1:numel(sweeps)
        vals = sweeps(s).values;
        contours = cell(numel(vals),1);
        for i = 1:numel(vals)

            Engine.TARGET_THRUST=BASE_THRUST; Engine.OF_RATIO=BASE_OF; Engine.P_CHAMBER=BASE_PC;
            Engine.P_EXIT=BASE_PR*Engine.P_AMBIENT; Engine.L_STAR=BASE_LSTAR;
            Engine.CONE_PERCENT=BASE_CONE; Engine.CONTRACTION_RATIO=BASE_CR;
            
            if strcmp(sweeps(s).field,'P_EXIT')
                Engine.P_EXIT = vals(i) * Engine.P_AMBIENT;
            else
                Engine.(sweeps(s).field) = vals(i) * sweeps(s).scale;
            end

            Engine = createEngine(Engine, GAS_DATA);
            X = Engine.X_COORDINATES*1e3;  R = Engine.R_COORDINATES*1e3;
            contours{i} = [X; R];
            xmax = max(xmax, max(X));  rmax = max(rmax, max(abs(R)));
            xmin = min(xmin, min(X));
        end

        data{s} = contours;

    end
    XL = [xmin*1.05, xmax*1.05];
    YL = [-rmax*1.05, rmax*1.05];
    
    for s = 1:numel(sweeps)
        nexttile;
        vals = sweeps(s).values;
        colors = nebula(numel(vals));
        hold on
        for i = 1:numel(vals)
            XR = data{s}{i};
            plot(XR(1,:),  XR(2,:), 'Color',colors(i,:), LineWidth=2);
            plot(XR(1,:), -XR(2,:), 'Color',colors(i,:), LineWidth=2);
        end
        axis equal; grid on; box on; xlim(XL); ylim(YL);
        xlabel('x (mm)'); ylabel('r (mm)'); title(sweeps(s).title);
        colormap(gca, nebula);
        if numel(vals) > 1
            d = vals(2)-vals(1);  clim([vals(1)-d/2, vals(end)+d/2]);
        end
        cb = colorbar; cb.Label.String = sweeps(s).label; cb.Ticks = vals;
    end

    subText = sprintf('\nBase values: Thrust = %.1f kN, O/F = %.1f, P_c = %.0f bar\nNPR = %.1f, L* = %.1f m, Cone = %d%%, CR = %d\n', ...
    BASE_THRUST/1000, BASE_OF, BASE_PC/1e5, BASE_PR, BASE_LSTAR, BASE_CONE, BASE_CR);

    title(tl, 'Nozzle Design Sweeps', 'FontWeight', 'bold', 'FontSize', 14);
    subtitle(tl, subText, 'FontWeight', 'normal', 'FontSize', 9);
    
    print(fig_engine_sweeps, 'nozzle_sweeps.pdf', '-dpdf');

end