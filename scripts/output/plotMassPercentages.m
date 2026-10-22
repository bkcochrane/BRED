function plotMassPercentages(Engine, PropData)

    x = Engine.Stream.x; 
    massPercentage = Engine.Stream.massPercentage;
    
    if nargin > 1 && isfield(PropData, 'SPECIES_NAMES')
        speciesNames = PropData.SPECIES_NAMES;
    else
        numSpecies = size(massPercentage, 2);
        speciesNames = arrayfun(@(x) sprintf('Species %d', x), 1:numSpecies, 'UniformOutput', false);
    end

    THRESHOLD = 0.5; 
    maxPctAlongEngine = max(massPercentage, [], 1);
    significantIdx = maxPctAlongEngine > THRESHOLD;
    
    plotPct = massPercentage(:, significantIdx);
    plotNames = speciesNames(significantIdx);

    figure('Name', 'Species Composition', 'Color', 'w', 'Position', [100, 100, 800, 500]);
    hold on; grid on;

    plot(x, plotPct, 'LineWidth', 2);

    xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);

    xlabel('x (m)', 'FontWeight', 'bold');
    ylabel('Mass Percentage (%)', 'FontWeight', 'bold');
    title('Exhaust Species Composition');
    box on
    legend(plotNames, 'Location', 'eastoutside', 'Box', 'off');
    
    xlim([min(x) max(x)]);
    ylim([0, 5]);
    
    hold off;
end