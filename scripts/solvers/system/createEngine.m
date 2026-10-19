function [Ideal, Engine] = createEngine(Engine, PropData)
    
    fprintf('Creating engine...\n')

    % ------------------------------------------------------------------------
    
    % Ideal Engine
    
    % ------------------------------------------------------------------------

    Ideal = createIdealEngine(Engine, PropData);

    % ------------------------------------------------------------------------
    
    % Real Engine
    
    % ------------------------------------------------------------------------
    
    Engine = createRealEngine(Ideal, PropData);
    
    % ------------------------------------------------------------------------
     
    % Full contour
    
    % ------------------------------------------------------------------------
    
    fprintf('\nCalculating engine geometry...\n')
    Engine = calculateProfile(Engine);
    fprintf('Engine geometry calculated!\n')

    fprintf('\nCalculating streamwise flow dynamics...\n')
    Engine = calculateStreamwise(Engine, PropData, 'full');
    fprintf('Flow dynamics calculated!\n')
    
    fprintf('\nSolving thermal circuit...\n')
    Engine = solveThermal(Engine, PropData);
    fprintf('Thermal circuit solved!\n')

    fprintf('\nEngine created!\n')

end