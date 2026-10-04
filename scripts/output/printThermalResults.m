function printThermalResults(Engine)
    
    fprintf('\n');
    fprintf(' Cooling\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f K\n', 'Max gas-side wall T', max(Engine.WALL_TWG));
    fprintf(' %-25s %10.2f K\n', 'Max coolant-side wall T', max(Engine.WALL_TWI));
    fprintf('\n');
    fprintf(' %-25s %10.2f kg/s\n', 'mDot Coolant', Engine.MDOT_COOLANT);
    fprintf(' %-25s %10.2f K\n', 'Coolant outlet T', max(Engine.WALL_TCOOLANT));
    fprintf('\n');
    fprintf(' %-25s %10.2f m\n',   'Channel length',    Engine.CHANNEL_LENGTH);
    fprintf(' %-25s %10.2f m/s\n', 'Peak coolant V',     Engine.COOLANT_V_MAX);
    fprintf(' %-25s %10.2f bar\n', 'Coolant dP',         Engine.COOLANT_DP/1e5);
    if isfield(Engine,'SPIRAL_RATE') && Engine.SPIRAL_RATE > 0
        fprintf(' %-25s %10d x %.1f turns\n', 'Spiral channels', Engine.NUM_CHANNELS, Engine.N_TURNS);
    end
    fprintf('=====================================================\n');

end