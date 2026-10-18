function Engine = solveThermal(Engine, PropData)
    
    Engine = calculateHeatFlux(Engine, PropData);
    Engine = sizeCoolingFlow(Engine, PropData);
    Engine = coolantPressureDrop(Engine);
    
end