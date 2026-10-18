function printThermalResults(Engine)
    
    fprintf('\n');
    fprintf(' Cooling\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f K\n', 'Max gas-side wall T', max(Engine.Wall.tempWG));
    fprintf(' %-25s %10.2f K\n', 'Max coolant-side wall T', max(Engine.Wall.tempWI));
    fprintf('\n');
    fprintf(' %-25s %10.2f kg/s\n', 'mDot Coolant', Engine.Cooling.Coolant.mdot);
    fprintf(' %-25s %10.2f K\n', 'Coolant outlet T', max(Engine.Wall.tempCoolant));
    fprintf('\n');
    fprintf(' %-25s %10d\n',   'Num. of channels',    Engine.Cooling.numChannels);
    fprintf(' %-25s %10.2f m\n',   'Channel length',    Engine.Cooling.lengthChannel);
    fprintf(' %-25s %10.2f m/s\n', 'Peak coolant V',     Engine.Cooling.Coolant.vMax);
    fprintf(' %-25s %10.2f bar\n', 'Coolant dP',         Engine.Cooling.Coolant.dP/1e5);
    if isfield(Engine,'SPIRAL_RATE') && Engine.SPIRAL_RATE > 0
        fprintf(' %-25s %10d x %.1f turns\n', 'Spiral channels', Engine.Cooling.numChannels, Engine.Cooling.numTurns);
    end
    fprintf('=====================================================\n');

end