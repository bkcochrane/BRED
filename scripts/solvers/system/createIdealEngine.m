function Ideal = createIdealEngine(Ideal, PropData)

    Ideal.Efficiency.combustion = 1.0;
    Ideal.Efficiency.nozzle = 1.0;
    Ideal.Efficiency.cdThroat = 1.0;
    Ideal.Geometry.LStar = Ideal.Geometry.LSTAR_TARGET;
    Ideal.Specs.ofRatio = Ideal.Specs.OF_RATIO_TARGET;

    fprintf('\nSolving for Ideal conditions...\n')

    Ideal = solveEngine(Ideal, PropData);
    Ideal = calculateProfile(Ideal);

    fprintf('Ideal engine solved!\n')

end