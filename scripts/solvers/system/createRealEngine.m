function Real = createRealEngine(Ideal, PropData)

    fprintf('\nSolving for Real engine...\n')
    
    Real = Ideal;

    MAX_ITER = 30;
    TOLERANCE = 1e-4;
    oldMdot =  Real.Specs.mdot;
    oldEta = NaN;
    
    for i = 1:MAX_ITER
        
        Real.Specs.mdot = oldMdot;

        % Calculate real physics
        Real = sizeInjector(Real, PropData);
        Real = solveRealEngine(Real, PropData);
        Real = calculateProfile(Real);
        Real = calculateStreamwise(Real, PropData, 'full');
        if Real.Cooling.REGEN_COOLING == true
            Real = solveThermal(Real, PropData);
        end

        % Check relative convergence
        newMdot = Real.Specs.mdot;
        newEta = Real.Efficiency.combustion;
        
        relErrMdot = abs(newMdot - oldMdot) / max(newMdot, 1e-12);
        if isnan(oldEta)
            relErrEta = Inf;
        else
            relErrEta = abs(newEta - oldEta) / max(newEta, 1e-12);
        end

        fprintf('Iteration %d       oldMdot: %.4f | newMdot: %.4f | Error: %.2e\n', i, oldMdot, newMdot, relErrMdot);
        
        if relErrMdot < TOLERANCE && relErrEta < TOLERANCE
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

        oldMdot = oldMdot + alpha * (newMdot - oldMdot);
        oldEta = newEta;
    end

    if Real.Cooling.REGEN_COOLING == true && Real.Cooling.COOLANT_LIMITED

        fprintf('[Warning] Engine cannot be regeneratively cooled: Increase film cooling or switch to an external cooling circuit.')

    end

    Real.Geometry.lock = true;

end