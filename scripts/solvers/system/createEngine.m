function Engine = createEngine(Engine, GAS_DATA)

    % ------------------------------------------------------------------------
    
    % Ideal Engine
    
    % ------------------------------------------------------------------------

    Engine.ETA_COMBUSTION = 1.0;
    Engine.ETA_NOZZLE = 1.0;
    Engine.DISCHARGE_COEFF = 1.0;
    
    Ideal = solveEngine(Engine, GAS_DATA);
    Ideal = calculateProfile(Ideal);
    
    printNozzleResults(Ideal, 'IDEAL');
    
    % ------------------------------------------------------------------------
    
    % Real Engine
    
    % ------------------------------------------------------------------------
    
    [Engine] = solveRealEngine(Ideal, GAS_DATA);
    
    printNozzleResults(Engine, 'REAL');
    
    Engine = sizeInjector(Engine);
    printInjector(Engine);
    
    % ------------------------------------------------------------------------
     
    % Full contour
    
    % ------------------------------------------------------------------------
    
    Engine = calculateProfile(Engine);
    Engine = calculateStreamwise(Engine, GAS_DATA, 'full');

end