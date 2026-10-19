function plotThrottleCurve(Engine, PropData)

   fprintf('\nCalculating throttle curves...\n')

    pchambers = linspace(Engine.Chamber.P * 0.1, Engine.Chamber.P, 100);
    
    thrustVacVec = zeros(size(pchambers));
    mdotFuelVec = zeros(size(pchambers));
    mdotOxVec = zeros(size(pchambers));
    
    MIN_OF = 0.5;
    DESIGN_OF = Engine.Specs.OF_RATIO_TARGET;
    
    for i = 1:length(pchambers)
        
        BaseEngine = Engine;
        BaseEngine.Specs.OF_RATIO_TARGET = DESIGN_OF;
        BaseEngine.Chamber.P = pchambers(i);
        
        fprintf('Thermal Balancing Pc = %.1f bar...\n', BaseEngine.Chamber.P/1e5);
        
        BaseEngine = solveEngine(BaseEngine, PropData);
        
        if Engine.Cooling.REGEN_COOLING == true

            BaseEngine.Coolant = PropData.Fuel;
            BaseEngine.Coolant.MDOT = BaseEngine.MDOT_FUEL;
        
            BaseEngine = solveThermal(BaseEngine, PropData);
            
            msgLength = 0;
            while BaseEngine.Coolant.MDOT > BaseEngine.MDOT_FUEL
                
                if BaseEngine.OF_RATIO <= MIN_OF
                    warning('off', 'all');
                    break;
                end
                
                BaseEngine.OF_RATIO = BaseEngine.OF_RATIO - 0.02;
                
                BaseEngine = solveEngine(BaseEngine, PropData);
                
                BaseEngine = solveThermal(BaseEngine, PropData);
                
                fprintf(repmat('\b', 1, msgLength));
                msg = sprintf('  Cooling demand too high -> Lowering O/F to %.2f', BaseEngine.OF_RATIO);
                fprintf('%s', msg);
                msgLength = length(msg);
                
            end
            
            if msgLength > 0
                fprintf('\n');
            end
            
            if BaseEngine.SPECS.OF_RATIO <= MIN_OF
                thrustVacVec(i) = NaN;
                mdotFuelVec(i) = NaN;
                mdotOxVec(i) = NaN;
            else
                CurrentEngine = BaseEngine;
                CurrentEngine.P_AMBIENT = 0;
                
                CurrentEngine = solveEngine(CurrentEngine, PropData);
                
                thrustVacVec(i) = CurrentEngine.THRUST;
                mdotFuelVec(i) = CurrentEngine.MDOT_FUEL;
                mdotOxVec(i) = CurrentEngine.MDOT_OX;
            end

        end
        
    end
    
    valid = ~isnan(thrustVacVec);
    P_CHAMBER_valid = pchambers(valid)';
    THRUST_valid = thrustVacVec(valid)';
    MDOT_FUEL_valid = mdotFuelVec(valid)';
    MDOT_OX_valid = mdotOxVec(valid)';
    
    figure('Name', 'Engine Performance Maps', 'Position', [100, 100, 1200, 400]);
    
    subplot(1,3,1);
    plot(P_CHAMBER_valid/1e5, THRUST_valid/1000, 'LineWidth', 2);
    grid on;
    xlabel('Chamber Pressure (bar)');
    ylabel('Vacuum Thrust (kN)');
    title('Vacuum Thrust Map');
    
    subplot(1,3,2);
    plot(P_CHAMBER_valid/1e5, MDOT_FUEL_valid, 'LineWidth', 2);
    grid on;
    xlabel('Chamber Pressure (bar)');
    ylabel('Fuel Mass Flow (kg/s)');
    title('Fuel Flow Map');
    
    subplot(1,3,3);
    plot(P_CHAMBER_valid/1e5, MDOT_OX_valid, 'LineWidth', 2);
    grid on;
    xlabel('Chamber Pressure (bar)');
    ylabel('Ox Mass Flow (kg/s)');
    title('Oxidizer Flow Map');
    
end
