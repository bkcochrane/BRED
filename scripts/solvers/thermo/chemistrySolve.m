function [hReactants, bCirc] = chemistrySolve(OF_RATIO, PropData)

    Fuel = PropData.Fuel;
    Ox = PropData.Ox;

    massFuel = 1 / (1 + OF_RATIO);
    massOx   = OF_RATIO / (1 + OF_RATIO);
    
    molesFuel = massFuel / Fuel.MOLECULAR_MASS;
    molesOx   = massOx / Ox.MOLECULAR_MASS;
    
    hReactants = (molesFuel * Fuel.HEAT_OF_FORMATION * 1e3) + (molesOx * Ox.HEAT_OF_FORMATION * 1e3);
    
    NUM_ELEM = numel(PropData.ELEMENTS);
    bCirc = zeros(NUM_ELEM, 1);
    
    for e = 1:NUM_ELEM
        sym = PropData.ELEMENTS{e};
        
        f_atoms = 0;
        if isfield(Fuel.Atoms, sym)
            f_atoms = Fuel.Atoms.(sym);
        end
        
        o_atoms = 0;
        if isfield(Ox.Atoms, sym)
            o_atoms = Ox.Atoms.(sym);
        end
        
        bCirc(e) = (f_atoms * molesFuel) + (o_atoms * molesOx);
    end
    
end