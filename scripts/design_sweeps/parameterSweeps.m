function parameterSweeps(Engine, PropData)

    fig_engine_sweeps = figure('Units','centimeters','Position',[2 2 21 29.7], ...
       'PaperType','A4','PaperOrientation','portrait', ...
       'PaperUnits','centimeters','PaperPosition',[0 0 21 29.7]);
    
    BASE_THRUST = Engine.TARGET_THRUST; 
    BASE_OF = Engine.OF_RATIO; 
    BASE_PC = Engine.P_CHAMBER; 
    BASE_LSTAR = Engine.L_STAR; 
    BASE_CONE = Engine.CONE_PERCENT; 
    BASE_CR = Engine.CONTRACTION_RATIO;
    
    sweep_mult = 0.5 : 0.25 : 1.5;

    sweeps(1) = struct('field','TARGET_THRUST', 'values', (BASE_THRUST/1000) * sweep_mult, 'label','Thrust (kN)', 'scale',1000, 'title','Varying thrust');
    
    if strcmpi(Engine.EXIT_MODE, 'pressure')
        BASE_PR = Engine.EXIT_PR;
        sweeps(2) = struct('field','P_EXIT', 'values', BASE_PR * sweep_mult, 'label','Exit pressure ratio', 'scale',1, 'title','Varying exit pressure ratio');
        exit_str = sprintf('NPR = %.1f', BASE_PR);
    else
        BASE_AR = Engine.AREA_RATIO_TARGET;
        sweeps(2) = struct('field','AREA_RATIO_TARGET', 'values', BASE_AR * sweep_mult, 'label','Area Ratio', 'scale',1, 'title','Varying Area Ratio');
        exit_str = sprintf('Area Ratio = %.1f', BASE_AR);
    end
    
    sweeps(3) = struct('field','P_CHAMBER', 'values', (BASE_PC/1e5) * sweep_mult, 'label','Chamber pressure (bar)', 'scale',1e5, 'title','Varying chamber pressure');
    
    sweeps(4) = struct('field','L_STAR', 'values', BASE_LSTAR * sweep_mult, 'label','L* (m)', 'scale',1, 'title','Varying L*');
    
    sweeps(5) = struct('field','CONTRACTION_RATIO', 'values', BASE_CR * sweep_mult, 'label','Contraction ratio', 'scale',1, 'title','Varying contraction ratio');
    
    tl = tiledlayout(numel(sweeps),1,'TileSpacing','compact','Padding','loose');

    xmax = 0; rmax = 0;
    xmin = 0;
    data = cell(numel(sweeps),1);
    
    for s = 1:numel(sweeps)
        vals = sweeps(s).values;
        contours = cell(numel(vals),1);
        
        for i = 1:numel(vals)

            fprintf('  Evaluating %s sweep: %d of %d...\n', sweeps(s).field, i, numel(vals));

            Engine.LOCK_GEOMETRY = false;
            
            Engine.TARGET_THRUST = BASE_THRUST; 
            Engine.OF_RATIO = BASE_OF; 
            Engine.P_CHAMBER = BASE_PC;
            Engine.L_STAR = BASE_LSTAR;
            Engine.CONE_PERCENT = BASE_CONE; 
            Engine.CONTRACTION_RATIO = BASE_CR;
            
            if strcmpi(Engine.EXIT_MODE, 'pressure')
                Engine.EXIT_PR = BASE_PR;
                Engine.P_EXIT = BASE_PR * Engine.P_AMBIENT;
            else
                Engine.AREA_RATIO_TARGET = BASE_AR;
            end
            
            if strcmp(sweeps(s).field, 'P_EXIT')
                Engine.EXIT_PR = vals(i) * sweeps(s).scale;
                Engine.P_EXIT = Engine.EXIT_PR * Engine.P_AMBIENT;
            else
                Engine.(sweeps(s).field) = vals(i) * sweeps(s).scale;
            end

            Engine = createEngine(Engine, PropData);
            
            X = Engine.X_COORDINATES*1e3;  
            R = Engine.R_COORDINATES*1e3;
            contours{i} = [X; R];
            
            xmax = max(xmax, max(X));  
            rmax = max(rmax, max(abs(R)));
            xmin = min(xmin, min(X));
            
        end
        fprintf('  Completed %s sweep!\n\n', sweeps(s).field);
        data{s} = contours;
    end
    
    XL = [xmin*1.05, xmax*1.05];
    YL = [-rmax*1.05, rmax*1.05];
    
    for s = 1:numel(sweeps)
        nexttile;
        vals = sweeps(s).values;
        if exist('nebula', 'file')
            colors = nebula(numel(vals));
            cmap = nebula;
        else
            colors = parula(numel(vals));
            cmap = parula;
        end
        hold on
        for i = 1:numel(vals)
            XR = data{s}{i};
            plot(XR(1,:),  XR(2,:), 'Color',colors(i,:), LineWidth=2);
            plot(XR(1,:), -XR(2,:), 'Color',colors(i,:), LineWidth=2);
        end
        axis equal; grid on; box on; xlim(XL); ylim(YL);
        xlabel('x (mm)'); ylabel('r (mm)'); title(sweeps(s).title);
        colormap(gca, cmap);
        
        if numel(vals) > 1
            d = vals(2)-vals(1);  
            clim([vals(1)-d/2, vals(end)+d/2]);
        end
        cb = colorbar; cb.Label.String = sweeps(s).label; cb.Ticks = vals;
    end

    subText = sprintf('\nBase values: Thrust = %.1f kN, O/F = %.1f, P_c = %.0f bar\n%s, L* = %.1f m, Cone = %d%%, CR = %d\n', ...
        BASE_THRUST/1000, BASE_OF, BASE_PC/1e5, exit_str, BASE_LSTAR, BASE_CONE, BASE_CR);

    title(tl, 'Nozzle Design Sweeps', 'FontWeight', 'bold', 'FontSize', 14);
    subtitle(tl, subText, 'FontWeight', 'normal', 'FontSize', 9);
    
    print(fig_engine_sweeps, 'nozzle_sweeps.pdf', '-dpdf');

end