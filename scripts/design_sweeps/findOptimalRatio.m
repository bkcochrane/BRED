function findOptimalRatio(Engine, GAS_DATA)
    
    % Define the three chamber pressures (Pa)
    P_NOMINAL = Engine.P_CHAMBER;
    P_LOW = P_NOMINAL - 5e5;
    P_HIGH = P_NOMINAL + 5e5;
    
    Pc_values = [P_LOW, P_NOMINAL, P_HIGH];
    line_styles = {'--', '-', ':'}; % Below = dashed, Nominal = solid, Above = dotted
    
    % Setup O/F vectors
    OF_RATIOs = Engine.OF_RATIO - 1.5 : 0.05 : Engine.OF_RATIO + 1.5;
    OF_fine = linspace(min(OF_RATIOs), max(OF_RATIOs), 500);
    
    % Preallocate storage for the fine sweeps (3 rows, one for each pressure)
    ISPs_fine_all = zeros(3, length(OF_fine));
    C_STARs_fine_all = zeros(3, length(OF_fine));
    
    % Create a base engine so we don't permanently overwrite the nominal properties
    base_engine = createEngine(Engine, GAS_DATA);
    
    % Loop through the 3 pressures
    for p_idx = 1:3
        local_engine = base_engine;
        local_engine.P_CHAMBER = Pc_values(p_idx);
        
        ISPs = zeros(size(OF_RATIOs));
        C_STARs = zeros(size(OF_RATIOs));
        
        for of = 1:length(OF_RATIOs)
            local_engine.OF_RATIO = OF_RATIOs(of);
            
            try
                loop_engine = solveEngine(local_engine, GAS_DATA);
                ISPs(of) = loop_engine.ISP;
                C_STARs(of) = loop_engine.C_STAR;
            catch ME
                fprintf('FAILED at Pc=%.1f bar, O/F=%.3f: %s\n', local_engine.P_CHAMBER/1e5, local_engine.OF_RATIO, ME.message);
                ISPs(of) = NaN; 
                C_STARs(of) = NaN;
            end
        end
        
        % Interpolate this specific pressure sweep
        ISPs_fine_all(p_idx, :) = interp1(OF_RATIOs, ISPs, OF_fine, 'spline');
        C_STARs_fine_all(p_idx, :) = interp1(OF_RATIOs, C_STARs, OF_fine, 'spline');
        
        % If this is the nominal run (index 2), extract the peaks
        if p_idx == 2
            [~, idx_isp] = max(ISPs_fine_all(p_idx, :));
            OF_PEAK_ISP = OF_fine(idx_isp);
            
            [~, idx_cstar] = max(C_STARs_fine_all(p_idx, :));
            OF_PEAK_CSTAR = OF_fine(idx_cstar);
        end
    end
    
    % ------------------------------------------------------------------------
    % Plotting
    % ------------------------------------------------------------------------
    figure('Units','centimeters','Position',[5 5 16 10]);
    neb = nebula(2);
    cm_isp = neb(1,:);
    cm_cstar = neb(2,:);
    
    % --- Left Axis (Isp) ---
    yyaxis left
    hold on;
    for p_idx = 1:3
        plot(OF_fine, ISPs_fine_all(p_idx, :), line_styles{p_idx}, 'Color', cm_isp, 'LineWidth', 2.5);
    end
    ylabel('I_{sp} (s)', 'Color', cm_isp);
    ax = gca;
    ax.YColor = cm_isp; 
    
    % Peak marker for nominal Isp
    xline(OF_PEAK_ISP, ':', 'Color', cm_isp, 'LineWidth', 2, ...
        'Label', sprintf('Peak I_{sp}: O/F = %.2f', OF_PEAK_ISP), ...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
    
    % --- Right Axis (C*) ---
    yyaxis right
    hold on;
    for p_idx = 1:3
        plot(OF_fine, C_STARs_fine_all(p_idx, :), line_styles{p_idx}, 'Color', cm_cstar, 'LineWidth', 2.5);
    end
    ylabel('c^* (m/s)', 'Color', cm_cstar);
    ax.YColor = cm_cstar;
    
    % Peak marker for nominal C*
    xline(OF_PEAK_CSTAR, ':', 'Color', cm_cstar, 'LineWidth', 2, ...
        'Label', sprintf('Peak c^*: O/F = %.2f', OF_PEAK_CSTAR), ...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
    
    % --- Formatting ---
    xlabel('O/F Ratio');
    ax.XColor = 'k';
    title(sprintf('Performance vs O/F Ratio (Nominal P_c = %.1f bar)', P_NOMINAL/1e5));
    
    axis tight;
    box on;
    grid on;
    
    % Create hidden "dummy lines" purely to generate a clean legend for the line styles
    h1 = plot(NaN, NaN, '--k', 'LineWidth', 2);
    h2 = plot(NaN, NaN, '-k', 'LineWidth', 2);
    h3 = plot(NaN, NaN, ':k', 'LineWidth', 2);
    legend([h1, h2, h3], {sprintf('%.1f bar', P_LOW/1e5), ...
                          sprintf('%.1f bar', P_NOMINAL/1e5), ...
                          sprintf('%.1f bar', P_HIGH/1e5)}, ...
                          'Location', 'best');
                          
    set(gcf, 'WindowStyle', 'docked');

end