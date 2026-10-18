function [Ideal, Engine] = createEngine(Engine, PropData)
    
    fprintf('Creating engine...\n')

    % ------------------------------------------------------------------------
    
    % Ideal Engine
    
    % ------------------------------------------------------------------------

    Engine.Efficiency.combustion = 1.0;
    Engine.Efficiency.nozzle = 1.0;
    Engine.Efficiency.cdThroat = 1.0;
    Engine.Geometry.LStar = Engine.Geometry.LSTAR_TARGET;
    Engine.Specs.ofRatio = Engine.Specs.OF_RATIO_TARGET;

    fprintf('\nSolving for Ideal conditions...\n')

    Ideal = solveEngine(Engine, PropData);
    Ideal = calculateProfile(Ideal);

    fprintf('Ideal engine solved!\n')
    
    % ------------------------------------------------------------------------
    
    % Real Engine
    
    % ------------------------------------------------------------------------
    
    fprintf('\nSolving for Real engine...\n')
    Real = Ideal;

    MAX_ITER = 20;
    TOLERANCE = 1e-4;
    oldMdot =  Real.Specs.mdot;
    
    for i = 1:MAX_ITER
        
        Real.Specs.mdot = oldMdot;

        % Calculate real physics
        Real = sizeInjector(Real, PropData);
        Real = solveRealEngine(Real, PropData);
        Real = calculateProfile(Real);
        Real = calculateStreamwise(Real, PropData, 'full');
        if Engine.Cooling.REGEN_COOLING == true
            Real = solveThermal(Real, PropData);
        end

        mdot = Real.Specs.mdot;
        
        % Check relative convergence
        newMdot = mdot;
        
        relErrMdot = abs(newMdot - oldMdot) / max(newMdot, 1e-12);

        fprintf('Iteration %d       oldMdot: %.4f | newMdot: %.4f | Error: %.2e\n', i, oldMdot, newMdot, relErrMdot);
        
        if relErrMdot < TOLERANCE
            fprintf('Real engine converged in %d iterations.\n', i);
            break;
        end
        
        if i == MAX_ITER
            warning('Real engine solver hit the maximum number of iterations without fully converging.');
        end
        
        % Apply relaxation factor if solver not converging
        if i > 5
            alpha = 0.7;
        else
            alpha = 1.0;
        end

        mdot = oldMdot + alpha * (newMdot - oldMdot);

        oldMdot = mdot;

    end

    Engine = Real;

    if Engine.Cooling.REGEN_COOLING == true && Engine.Cooling.Coolant.tempMargin < 0

        error('Engine cannot be regeneratively cooled: Increase film cooling or switch to an external cooling circuit.')

    end

    Engine.Geometry.LOCK = true;
    
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