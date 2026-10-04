function printNozzleResults(Engine, NAME)

    fprintf('\n');
    fprintf('=====================================================\n');
    fprintf(' %s ENGINE\n', NAME);
    fprintf('=====================================================\n');
    fprintf('\n');
    fprintf(' Station     P (bar)    T (K)\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-10s %8.3g %9.0f\n', 'Chamber', Engine.P_CHAMBER / 1e5, Engine.T_CHAMBER);
    fprintf(' %-10s %8.2f %9.0f      Mach = %.3f\n', 'Throat', Engine.P_THROAT / 1e5, Engine.T_THROAT, Engine.MACH_THROAT);
    fprintf(' %-10s %8.2f %9.0f    V_EXIT = %.1f m/s\n', 'Exit', Engine.P_EXIT / 1e5, Engine.T_EXIT, Engine.V_EXIT);
    fprintf('=====================================================\n');
    fprintf('\n');
    fprintf(' Performance\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f\n', 'O/F ratio', Engine.OF_RATIO);
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f kg/s\n', 'mDot total', Engine.MDOT);
    fprintf(' %-25s %10.2f kg/s\n', 'mDot oxidiser', Engine.MDOT_OX);
    fprintf(' %-25s %10.2f kg/s\n', 'mDot fuel', Engine.MDOT_FUEL);
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f kN\n', 'Thrust', Engine.THRUST/1e3);
    fprintf(' %-25s %10.1f s\n', 'Isp', Engine.ISP);
    fprintf(' %-25s %10.2f m/s\n', 'C*', Engine.C_STAR);
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10.2f m\n', 'L*', Engine.L_STAR);
    fprintf(' %-25s %10.2f mm^2\n', 'Chamber volume', Engine.VOLUME_CHAMBER * 1000000);
    fprintf(' %-25s %10.2f mm\n', 'Chamber diameter', Engine.DIAMETER_CHAMBER * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Chamber length', Engine.LENGTH_CHAMBER * 1000);
    fprintf('\n');
    fprintf(' %-25s %10.2f\n', 'Area ratio', Engine.AREA_RATIO);
    fprintf(' %-25s %10.2f mm\n', 'Throat diameter', Engine.DIAMETER_THROAT * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Exit diameter', Engine.DIAMETER_EXIT * 1000);
    fprintf('\n');
    fprintf(' %-25s %10.2f mm\n', 'Nozzle length', Engine.NOZZLE_LENGTH * 1000);
    fprintf('=====================================================\n');

end