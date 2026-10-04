function Engine = solveThermal(Engine, GAS_DATA)
    
    Engine = calculateHeatFlux(Engine, GAS_DATA);
    Engine = sizeCoolingFlow(Engine, GAS_DATA);
    Engine = coolantPressureDrop(Engine);
    printThermalResults(Engine);

end