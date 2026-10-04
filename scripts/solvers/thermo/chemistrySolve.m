function [H_REACTANTS, B_I_CIRC] = chemistrySolve(OF_RATIO, GasData)

    MASS_FUEL = 1 / (1 + OF_RATIO);
    MASS_OX   = OF_RATIO / (1 + OF_RATIO);
    
    Ox = struct('atoms', struct());
    switch upper(GasData.OXIDIZER)
        case 'LOX'
            Ox.M = 31.9988;
            Ox.Hf = -12979.0e3;  % Liquid at 90 K
            Ox.atoms.O = 2;
        case 'HTP' 
            % 90% H2O2 / 10% H2O by mass approximation
            Ox.M = 32.65; 
            Ox.Hf = -197300e3;   % Liquid at 298.15 K
            Ox.atoms.H = 2.12; 
            Ox.atoms.O = 1.94;
        otherwise
            error('Unknown oxidizer: %s', GasData.OXIDIZER);
    end
    
    Fuel = struct('atoms', struct());
    switch upper(GasData.FUEL)
        case 'RP-1'
            Fuel.M = 13.976; % Empirical CH_1.95
            Fuel.Hf = -24717.7e3; % Liquid at 298.15 K
            Fuel.atoms.C = 1;
            Fuel.atoms.H = 1.95;
        case {'METHANE', 'CH4'}
            Fuel.M = 16.043;
            Fuel.Hf = -89360e3;  % Liquid at 111 K
            Fuel.atoms.C = 1;
            Fuel.atoms.H = 4;
        case 'LH2'
            Fuel.M = 2.016;
            Fuel.Hf = -9012e3;   % Liquid at 20 K
            Fuel.atoms.H = 2;
        otherwise
            error('Unknown fuel: %s', GasData.FUEL);
    end
    
    MOLES_FUEL = MASS_FUEL / Fuel.M;
    MOLES_OX   = MASS_OX / Ox.M;
    
    H_REACTANTS = (MOLES_FUEL * Fuel.Hf) + (MOLES_OX * Ox.Hf);
    
    NUM_ELEM = numel(GasData.ELEMENTS);
    B_I_CIRC = zeros(NUM_ELEM, 1);
    
    for e = 1:NUM_ELEM
        sym = GasData.ELEMENTS{e};
        
        f_atoms = 0;
        if isfield(Fuel.atoms, sym)
            f_atoms = Fuel.atoms.(sym);
        end
        
        o_atoms = 0;
        if isfield(Ox.atoms, sym)
            o_atoms = Ox.atoms.(sym);
        end
        
        B_I_CIRC(e) = (f_atoms * MOLES_FUEL) + (o_atoms * MOLES_OX);
    end
    
end