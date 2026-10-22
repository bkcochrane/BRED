function parameterSweeps(Engine, PropData)

    BASE_THRUST = Engine.Specs.TARGET_THRUST; 
    BASE_OF = Engine.Specs.OF_RATIO_TARGET; 
    BASE_PC = Engine.Chamber.P; 
    BASE_LSTAR = Engine.Geometry.LSTAR_TARGET; 
    BASE_CONE = Engine.Geometry.CONE_PERCENT; 
    BASE_CR = Engine.Geometry.CONTRACTION_RATIO_TARGET;
    
    sweepMultiplier = 0.5 : 0.1 : 1.5;

    sweeps(1) = struct('field','Specs.TARGET_THRUST', ...
                       'values', (BASE_THRUST/1000) * sweepMultiplier, ...
                       'label','Thrust (kN)', 'scale',1000, 'title','Varying thrust');
    
    if strcmpi(Engine.Specs.EXIT_MODE, 'pressure')
        BASE_PR = Engine.Specs.EXIT_PR;
        sweeps(2) = struct('field','Specs.EXIT_PR', ...
                           'values', BASE_PR * sweepMultiplier, ...
                           'label','Exit pressure ratio', 'scale',1, 'title','Varying exit pressure ratio');
        exitStr = sprintf('NPR = %.1f', BASE_PR);
    else
        BASE_AR = Engine.Specs.AREA_RATIO_TARGET;
        sweeps(2) = struct('field','Specs.AREA_RATIO_TARGET', ...
                           'values', BASE_AR * sweepMultiplier, ...
                           'label','Area Ratio', 'scale',1, 'title','Varying Area Ratio');
        exitStr = sprintf('Area Ratio = %.1f', BASE_AR);
    end
    
    sweeps(3) = struct('field','Chamber.P', ...
                       'values', (BASE_PC/1e5) * sweepMultiplier, ...
                       'label','Chamber pressure (bar)', 'scale',1e5, 'title','Varying chamber pressure');
    
    sweeps(4) = struct('field','Geometry.LSTAR_TARGET', ...
                       'values', BASE_LSTAR * sweepMultiplier, ...
                       'label','L* (m)', 'scale',1, 'title','Varying L*');
    
    sweeps(5) = struct('field','Geometry.CONTRACTION_RATIO_TARGET', ...
                       'values', BASE_CR * sweepMultiplier, ...
                       'label','Contraction ratio', 'scale',1, 'title','Varying contraction ratio');
    
    xmax = 0;
    rmax = 0;
    xmin = 0;
    data = cell(numel(sweeps),1);
    
    for s = 1:numel(sweeps)
        vals = sweeps(s).values;
        contours = cell(numel(vals),1);
        
        for i = 1:numel(vals)

            fprintf('\n  Evaluating %s sweep: %d of %d...\n', sweeps(s).field, i, numel(vals));

            % Reset base states prior to perturbation
            Engine.Geometry.lock = false;
            Engine.Specs.TARGET_THRUST = BASE_THRUST; 
            Engine.Specs.OF_RATIO_TARGET = BASE_OF; 
            Engine.Chamber.P = BASE_PC;
            Engine.Geometry.LSTAR_TARGET = BASE_LSTAR;
            Engine.Geometry.CONE_PERCENT = BASE_CONE; 
            Engine.Geometry.CONTRACTION_RATIO_TARGET = BASE_CR;
            
            if strcmpi(Engine.Specs.EXIT_MODE, 'pressure')
                Engine.Specs.EXIT_PR = BASE_PR;
                Engine.Exit.p = BASE_PR * Engine.Specs.P_AMBIENT;
            else
                Engine.Specs.AREA_RATIO_TARGET = BASE_AR;
            end
            
            % Assign target sweep value across nested structs
            valToSet = vals(i) * sweeps(s).scale;
            Engine = setNestedField(Engine, sweeps(s).field, valToSet);
            
            % Maintain coupled pressure state if sweeping exit pressure ratio
            if strcmp(sweeps(s).field, 'Specs.EXIT_PR')
                Engine.Exit.p = Engine.Specs.EXIT_PR * Engine.Specs.P_AMBIENT;
            end

            Engine = createIdealEngine(Engine, PropData);
            
            x = Engine.Geometry.xCoordinates * 1e3;  
            r = Engine.Geometry.rCoordinates * 1e3;
            contours{i} = [x; r];
            
            xmax = max(xmax, max(x));  
            rmax = max(rmax, max(abs(r)));
            xmin = min(xmin, min(x));
            
        end
        fprintf('\n  Completed %s sweep!\n\n', sweeps(s).field);
        data{s} = contours;
    end

    fprintf('=====================================================\n\n');
    fprintf('Parameter sweeps completed!\n\n')
    fprintf('=====================================================\n');

    figEngineSweeps = figure('Units','centimeters','Position',[2 2 21 29.7], ...
       'PaperType','A4','PaperOrientation','portrait', ...
       'PaperUnits','centimeters','PaperPosition',[0 0 21 29.7]);
    
    tl = tiledlayout(figEngineSweeps, numel(sweeps), 1, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    xL = [xmin*1.05, xmax*1.05];
    yL = [-rmax*1.05, rmax*1.05];
    
    for s = 1:numel(sweeps)
        nexttile(tl);
        vals = sweeps(s).values;
        nVals = numel(vals);
        
        if exist('nebula', 'file')
            colors = nebula(nVals);
            cmap = nebula(256);
        else
            colors = turbo(nVals);
            cmap = turbo(256);
        end
        
        chamberLengths = cellfun(@(c) max(c(1,:)) - min(c(1,:)), data{s});
        [~, drawOrder] = sort(chamberLengths, 'descend');

        hold on
        for idx = 1:nVals
            i = drawOrder(idx);
            XR = data{s}{i};
            plot(XR(1,:),  XR(2,:), 'Color', colors(i,:), 'LineWidth', 2);
            plot(XR(1,:), -XR(2,:), 'Color', colors(i,:), 'LineWidth', 2);
        end
        hold off
        
        axis equal; grid on; box on; xlim(xL); ylim(yL);
        xlabel('x (mm)'); ylabel('r (mm)'); title(sweeps(s).title);
        colormap(gca, cmap);
        
        if nVals > 1
            d = vals(2) - vals(1);  
            clim([vals(1) - d/2, vals(end) + d/2]);
        end
        cb = colorbar; 
        cb.Label.String = sweeps(s).label; 
        cb.Ticks = linspace(min(vals), max(vals), 5);
    end

    subText = sprintf('\nBase values: Thrust = %.1f kN, O/F = %.1f, P_c = %.0f bar\n%s, L* = %.1f m, Cone = %d%%, CR = %d\n', ...
        BASE_THRUST/1000, BASE_OF, BASE_PC/1e5, exitStr, BASE_LSTAR, BASE_CONE, BASE_CR);

    title(tl, 'Nozzle Design Sweeps', 'FontWeight', 'bold', 'FontSize', 14);
    subtitle(tl, subText, 'FontWeight', 'normal', 'FontSize', 9);
    
    print(figEngineSweeps, 'nozzle_sweeps.pdf', '-dpdf');
end

% Helper function to set fields separated by dots (e.g. 'Specs.TARGET_THRUST')
function S = setNestedField(S, fieldPath, value)
    parts = strsplit(fieldPath, '.');
    if isscalar(parts)
        S.(parts{1}) = value;
    else
        subField = parts{1};
        rest = strjoin(parts(2:end), '.');
        S.(subField) = setNestedField(S.(subField), rest, value);
    end
end