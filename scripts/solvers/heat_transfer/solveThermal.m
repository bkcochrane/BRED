function Engine = solveThermal(Engine, PropData)
    
    Engine = calculateHeatFlux(Engine, PropData, 'Bartz');
    Engine = sizeCoolingFlow(Engine, PropData);
    Engine = coolantPressureDrop(Engine);
    Engine.Wall.Bartz = Engine.Wall;

    Engine = calculateHeatFlux(Engine, PropData, 'Ievlev');
    Engine.Wall.Ievlev = Engine.Wall;

end