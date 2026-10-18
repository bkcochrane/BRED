function printNozzleResults(Engine, name)

    fprintf('\n');
    fprintf('=====================================================\n');
    fprintf(' %s ENGINE\n', name);
    fprintf('=====================================================\n');
    fprintf('\n');
    fprintf(' Station     P (bar)    T (K)\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-10s %8.3g %9.0f\n', 'Chamber', Engine.Chamber.P / 1e5, Engine.Chamber.temp);
    fprintf(' %-10s %8.2f %9.0f     Mach = %.3f\n', 'Throat', Engine.Throat.p / 1e5, Engine.Throat.temp, Engine.Throat.mach);
    fprintf(' %-10s %8.2f %9.0f    vExit = %.1f m/s\n', 'Exit', Engine.Exit.p / 1e5, Engine.Exit.temp, Engine.Exit.v);
    fprintf('=====================================================\n');
    fprintf('\n');
    fprintf(' Performance\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f\n', 'O/F ratio', Engine.Specs.ofRatio);
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f kg/s\n', 'mDot total', Engine.Specs.mdot);
    fprintf(' %-25s %10.2f kg/s\n', 'mDot oxidiser', Engine.Ox.mdot);
    fprintf(' %-25s %10.2f kg/s\n', 'mDot fuel', Engine.Fuel.mdot);
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f kN\n', 'Thrust', Engine.Specs.thrust / 1e3);
    fprintf(' %-25s %10.1f s\n', 'Isp', Engine.Specs.Isp);
    fprintf(' %-25s %10.2f m/s\n', 'C*', Engine.Specs.cStar);
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f m\n', 'L*', Engine.Geometry.LStar);
    fprintf(' %-25s %10.2f mm^3\n', 'Chamber volume', Engine.Chamber.volume * 1e9);
    fprintf(' %-25s %10.2f mm\n', 'Chamber diameter', Engine.Chamber.diameter * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Chamber length', Engine.Geometry.lengthChamber * 1000);
    fprintf('\n');
    fprintf(' %-25s %10.2f\n', 'Area ratio', Engine.Specs.areaRatio);
    fprintf(' %-25s %10.2f mm\n', 'Throat diameter', Engine.Throat.diameter * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Exit diameter', Engine.Exit.diameter * 1000);
    fprintf('\n');
    fprintf(' %-25s %10.2f mm\n', 'Nozzle length', Engine.Geometry.lengthNozzle * 1000);
    fprintf('\n');
    fprintf(' %-25s %10.2f mm\n', 'Engine length', Engine.Geometry.lengthEngine * 1000);
    fprintf('=====================================================\n');

end