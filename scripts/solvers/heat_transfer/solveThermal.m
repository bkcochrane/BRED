function Engine = solveThermal(Engine, PropData)

    filmPercent = Engine.Cooling.FILM_COOLING_FRACTION;
    fuelMdot = Engine.Fuel.mdot;

    if filmPercent > 0

        Engine.Cooling.Film.mdot = fuelMdot * filmPercent;

    end
    
    if Engine.Cooling.REGEN_COOLING == true

        Engine.Cooling.Coolant = PropData.Fuel;
        Engine.Cooling.Coolant.mdot = fuelMdot * (1 - filmPercent);
        Engine.Cooling.Coolant.T_SUPPLY = PropData.Fuel.T_SUPPLY;

    else

        Engine.Cooling.Coolant.P_REF = [1, 2, 5, 10, 20, 50, 100, 150, 220.6]; 
        Engine.Cooling.Coolant.T_REF = [373.1, 393.4, 425.0, 453.2, 485.6, 537.1, 584.3, 615.4, 647.1]; 
        
        Engine.Cooling.Coolant.tSat = @(P) interp1(P_REF, T_REF, P/1e5, 'pchip', 'extrap');
        
        tempSupply = 298; % K
        Engine.Cooling.Coolant.T_SUPPLY = tempSupply;
        
        Engine.Cooling.Coolant.RHO = 1000.0; % kg/m3
        Engine.Cooling.Coolant.MU = 8.9e-4; % Pa*s
        Engine.Cooling.Coolant.CP = 4184.0; % J/(kg*K)
        Engine.Cooling.Coolant.K = 0.606; % W/(m*K)

    end
    
    Engine = sizeCoolingFlow(Engine, PropData, 'Bartz');

    EngingBartz  = calculateHeatFlux(Engine, PropData, 'Bartz');
    EngineIevlev = calculateHeatFlux(Engine, PropData, 'Ievlev');

    Engine.Wall = EngineIevlev.Wall;
    Engine.Wall.Bartz = EngingBartz.Wall;
    Engine.Wall.Ievlev = EngineIevlev.Wall;

    Engine = coolantPressureDrop(Engine);

end