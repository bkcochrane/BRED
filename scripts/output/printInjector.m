function printInjector(Engine)

    fprintf('\n');
    fprintf(' Injector\n');
    fprintf('-----------------------------------------------------\n');
    fprintf(' %-25s %10d\n', 'Impinging sets', Engine.Injector.numPortSets);
    fprintf(' %-25s %10d\n', 'Fuel ports', Engine.Injector.numPortsFuel);
    fprintf(' %-25s %10d\n', 'Oxidiser ports', Engine.Injector.numPortsOx);
    fprintf('\n');
    fprintf(' %-25s %10.2f mm\n', 'Fuel port diameter', Engine.Injector.diameterPortsFuel * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Ox port diameter', Engine.Injector.diameterPortsOx * 1000);
    fprintf(' %-25s %10.1f deg\n','Impingement angle', rad2deg(Engine.Injector.impingementAngle));
    fprintf('\n');
    fprintf(' %-25s %10.2f mm\n', 'Inner ring diameter', Engine.Injector.innerRingDiameter * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Middle ring diameter', Engine.Injector.middleRingDiameter * 1000);
    fprintf(' %-25s %10.2f mm\n', 'Outer ring diameter', Engine.Injector.outerRingDiameter * 1000);
    fprintf('=====================================================\n');

end