function findOptimalRatio(Engine, PropData)
    
    fprintf('\nFinding optimal O/F ratio...\n');

    % Define the three chamber pressures (Pa)
    pNominal = Engine.Chamber.P;
    pLow = pNominal - 5e5;
    pHigh = pNominal + 5e5;
    
    pValues = [pLow, pNominal, pHigh];
    
    % Setup O/F vectors
    ofRatios = 0.2: 0.2 : 8;
    ofFine = linspace(min(ofRatios), max(ofRatios), 500);
    
    % Preallocate storage for the fine sweeps
    ispFine = zeros(3, length(ofFine));
    cStarFine = zeros(3, length(ofFine));
    
    % Create a base engine
    fprintf('Creating base engine...\n')
    BaseEngine = createEngine(Engine, PropData);
    
    for pIdx = 1:3
        LocalEngine = BaseEngine;
        LocalEngine.Chamber.P = pValues(pIdx);
        Isp = zeros(size(ofRatios));
        cStar = zeros(size(ofRatios));

        fprintf('Solving engine for Pc=%.1f bar...\n', LocalEngine.Chamber.P/1e5)

        for of = 1:length(ofRatios)
            LocalEngine.Specs.ofRatio = ofRatios(of);
            try
                LoopEngine = solveEngine(LocalEngine, PropData);
                Isp(of) = LoopEngine.Specs.Isp;
                cStar(of) = LoopEngine.Specs.cStar;
            catch ME
                fprintf('FAILED at Pc=%.1f bar, O/F=%.3f: %s\n', LocalEngine.Chamber.P/1e5, LocalEngine.Specs.ofRatio, ME.message);
                Isp(of) = NaN; 
                cStar(of) = NaN;
            end
        end

        % Strip NaNs from solver failures
        validIdx = ~isnan(Isp) & ~isnan(cStar);
        
        % Interpolate
        ispFine(pIdx, :) = interp1(ofRatios(validIdx), Isp(validIdx), ofFine, 'spline');
        cStarFine(pIdx, :) = interp1(ofRatios(validIdx), cStar(validIdx), ofFine, 'spline');
        
        % Find the peaks for the nominal pressure line
        if pIdx == 2
            
            % ---
            % Find Peak Isp
            % ---
            SEARCH_RADIUS = 1.0;
            while SEARCH_RADIUS <= 3.0 % Hard limit to prevent infinite loops
                
                % Define bounds
                bounds = find(ofFine >= (Engine.Specs.OF_RATIO_TARGET - SEARCH_RADIUS) & ofFine <= (Engine.Specs.OF_RATIO_TARGET + SEARCH_RADIUS));
                
                % Find max within bounds
                [~, localMax] = max(ispFine(pIdx, bounds));
                IspIdx = bounds(localMax);
                
                % If the max is sitting exactly on the boundary edge, the gradient is not 0 (it is still climbing).
                % If it is completely inside the bounds, it is a true peak where the gradient crosses 0.
                if localMax > 1 && localMax < length(bounds)
                    break;
                else
                    SEARCH_RADIUS = SEARCH_RADIUS + 0.1; % Extend bounds and try again
                end
            end
            ofPeakIsP = ofFine(IspIdx);
            
            
            % ---
            % Find Peak C*
            % ---
            SEARCH_RADIUS = 1.0;
            while SEARCH_RADIUS <= 3.0 
                
                bounds = find(ofFine >= (Engine.Specs.OF_RATIO_TARGET - SEARCH_RADIUS) & ofFine <= (Engine.Specs.OF_RATIO_TARGET + SEARCH_RADIUS));
                
                [~, localMax] = max(cStarFine(pIdx, bounds));
                cStarIdx = bounds(localMax);
                
                if localMax > 1 && localMax < length(bounds)
                    break;
                else
                    SEARCH_RADIUS = SEARCH_RADIUS + 0.1;
                end
            end
            ofPeakCStar = ofFine(cStarIdx);
            
            fprintf('  -> Optimal O/F ratios found: Peak Isp = %.2f | Peak C* = %.2f\n', ofPeakIsP, ofPeakCStar);
        end
        
    end
    
    % ------------------------------------------------------------------------
    % Plotting
    % ------------------------------------------------------------------------
    figure('Units','centimeters','Position',[5 5 16 10]);
    if exist('nebula', 'file')
        COLOR_THEME = nebula(2);
    else
        COLOR_THEME = parula(2);
    end
    
    cmIsp = COLOR_THEME(1,:);
    cmCStar = COLOR_THEME(2,:);

    lineStyles = {'--', '-', ':'};
    
    % Isp
    yyaxis left
    hold on;
    for pIdx = 1:3
        plot(ofFine, ispFine(pIdx, :), lineStyles{pIdx}, 'Color', cmIsp, 'LineWidth', 2.5);
    end
    ylabel('I_{sp} (s)', 'Color', cmIsp);
    ax = gca;
    ax.YColor = cmIsp; 
    
    % Peak marker for nominal Isp
    xline(ofPeakIsP, ':', 'Color', cmIsp, 'LineWidth', 2, ...
        'Label', sprintf('Peak I_{sp}: O/F = %.2f', ofPeakIsP), ...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
    
    % C*
    yyaxis right
    hold on;
    for pIdx = 1:3
        plot(ofFine, cStarFine(pIdx, :), lineStyles{pIdx}, 'Color', cmCStar, 'LineWidth', 2.5);
    end
    ylabel('c^* (m/s)', 'Color', cmCStar);
    ax.YColor = cmCStar;
    
    % Peak marker for nominal C*
    xline(ofPeakCStar, ':', 'Color', cmCStar, 'LineWidth', 2, ...
        'Label', sprintf('Peak c^*: O/F = %.2f', ofPeakCStar), ...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
    
    xlabel('O/F Ratio');
    ax.XColor = 'k';
    title(sprintf('Performance vs O/F Ratio (Nominal P_c = %.1f bar)', pNominal/1e5));
    
    axis tight;
    box on;
    grid on;
    
    h1 = plot(NaN, NaN, '--k', 'LineWidth', 2);
    h2 = plot(NaN, NaN, '-k', 'LineWidth', 2);
    h3 = plot(NaN, NaN, ':k', 'LineWidth', 2);
    legend([h1, h2, h3], {sprintf('%.1f bar', pLow/1e5), ...
                          sprintf('%.1f bar', pNominal/1e5), ...
                          sprintf('%.1f bar', pHigh/1e5)}, ...
                          'Location', 'best');
                          
    set(gcf, 'WindowStyle', 'docked');

end