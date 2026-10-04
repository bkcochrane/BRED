function printInjector(Engine)

    fprintf('\n');
    fprintf(' Injector\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10d\n', 'Triplet sets', Engine.INJECTOR_N_PORT_SETS);
    fprintf(' %-25s %10d\n', 'Fuel ports', Engine.INJECTOR_N_PORTS_FUEL);
    fprintf(' %-25s %10d\n', 'Oxidiser ports', Engine.INJECTOR_N_PORTS_OX);
    fprintf('\n');
    fprintf(' %-25s %10.2f mm\n', 'Fuel port diameter', Engine.INJECTOR_DIAMETER_FUEL * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Ox port diameter', Engine.INJECTOR_DIAMETER_OX * 1000);
    fprintf(' %-25s %10.1f deg\n','Impingement angle', rad2deg(Engine.INJECTOR_IMPINGEMENT_ANGLE));
    fprintf('=====================================================\n');

end