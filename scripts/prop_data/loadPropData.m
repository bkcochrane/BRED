function PropData = loadPropData(FUEL_NAME, OX_NAME)

    if nargin < 1, FUEL_NAME = 'RP-1'; end
    if nargin < 2, OX_NAME = 'LOX'; end

    Fuel = struct();
    Fuel.Atoms = struct();
    
    switch upper(FUEL_NAME)
        case {'RP-1', 'RP1'}
            Fuel.MOLECULAR_MASS = 13.976; % CH_1.95
            Fuel.HEAT_OF_FORMATION = -24717.7e3; % Liquid at 298.15 K (J/mol)
            Fuel.Atoms.C = 1;
            Fuel.Atoms.H = 1.95;

            Fuel.P_REF = [1 2 5 10 15 20 22.6]; 
            Fuel.T_REF = [450 480 525 565 595 620 645]; 
            Fuel.tSat = @(p) interp1(P, T, p/1e5, 'pchip', 'extrap');
            Fuel.T_SUPPLY = 298; % K

            Fuel.RHO = 820.0; % kg/m3
            Fuel.MU = 1.48e-3; % Pa*s
            Fuel.CP = 1880.0; % J/(kg*K)
            Fuel.K = 0.145; % W/(m*K)
            Fuel.SURFACE_TENSION = 0.024; % N/m

            Fuel.QS = 246e3; % J/kg
            Fuel.TEMP_BOIL = 490; % K
            Fuel.calcK = @(t) 0.002 + 0.0001 * t; % Vapor thermal conductivity

        case {'METHANE', 'CH4', 'LCH4'}
            Fuel.MOLECULAR_MASS = 16.043;
            Fuel.HEAT_OF_FORMATION = -89360; % Liquid at 111 K (J/mol)
            Fuel.Atoms.C = 1;
            Fuel.Atoms.H = 4;

            Fuel.P_REF = [1 5 10 15 20 25 30 35 40 45.99];
            Fuel.T_REF = [111.6 131.3 143.4 152.6 160.4 167.3 173.6 179.4 184.9 190.56];
            Fuel.tSat = @(p) interp1(P, T, p/1e5, 'pchip', 'extrap');
            Fuel.T_SUPPLY = 95; % K

            Fuel.RHO = 445.8; % kg/m3
            Fuel.MU = 1.95e-4; % Pa*s
            Fuel.CP = 3350.0; % J/(kg*K)
            Fuel.K = 0.218; % W/(m*K)
            Fuel.SURFACE_TENSION = 0.014; % N/m
            
            Fuel.QS = 510e3; % J/kg
            Fuel.TEMP_BOIL = 111; % K
            Fuel.calcK = @(t) 0.010 + 0.0001 * t;

        case {'LH2', 'HYDROGEN', 'H2'}
            Fuel.MOLECULAR_MASS = 2.016;
            Fuel.HEAT_OF_FORMATION = -9012; % Liquid at 20 K (J/mol) ->
            Fuel.Atoms.H = 2;

            Fuel.P_REF = [1 2 3 4 5 6 8 10 12.96];
            Fuel.T_REF = [20.3 22.8 24.5 25.9 27.0 28.0 29.8 31.3 32.9];
            Fuel.tSat = @(p) interp1(P, T, p/1e5, 'pchip', 'extrap');
            Fuel.T_SUPPLY = 20; % K

            Fuel.RHO = 70.85; % kg/m3
            Fuel.MU = 1.35e-5; % Pa*s
            Fuel.CP = 9630.0; % J/(kg*K)
            Fuel.K = 0.104; % W/(m*K)
            Fuel.SURFACE_TENSION = 0.002; % N/m
            
            Fuel.QS = 446e3; % J/kg
            Fuel.TEMP_BOIL = 20; % K
            Fuel.calcK = @(t) 0.050 + 0.0005 * t;

        case {'IPA', 'ISOPROPANOL', 'C3H8O'}
            Fuel.MOLECULAR_MASS = 60.096;
            Fuel.HEAT_OF_FORMATION = -318100; % Liquid at 298.15 K (J/mol)
            Fuel.Atoms.C = 3;
            Fuel.Atoms.H = 8;
            Fuel.Atoms.O = 1;

            Fuel.P_REF = [1, 2, 5, 10, 15, 20, 25, 30, 40, 47.6];
            Fuel.T_REF = [355.4, 373.5, 403.5, 430.8, 449.2, 463.8, 476.0, 486.6, 504.0, 508.3];
            Fuel.tSat = @(p) interp1(P, T, p/1e5, 'pchip', 'extrap');
            Fuel.T_SUPPLY = 293; % K

            Fuel.RHO = 785.0; % kg/m3
            Fuel.MU = 1.30e-3; % Pa*s
            Fuel.CP = 2680.0; % J/(kg*K)
            Fuel.K = 0.137; % W/(m*K)
            Fuel.SURFACE_TENSION = 0.0217; % N/m
    
            Fuel.QS = 663e3; % J/kg
            Fuel.TEMP_BOIL = 355.4; % K
            Fuel.calcK = @(t) max(0.005, -0.012 + 0.00011 * t);

        otherwise
            error('loadGasData:UnknownFuel', 'Unknown fuel: %s', FUEL_NAME);
    end
    
    Ox = struct();
    Ox.Atoms = struct(); 
    
    switch upper(OX_NAME)
        case 'LOX'
            Ox.MOLECULAR_MASS = 31.9988;
            Ox.HEAT_OF_FORMATION = -12979e3; % Liquid at 95 K (J/mol)
            Ox.Atoms.O = 2;

            Ox.RHO = 1141.0; % kg/m3
            Ox.SURFACE_TENSION = 0.013; % N/m
            
            % Vaporization Properties
            Ox.QS = 213e3; % J/kg
            Ox.TEMP_BOIL = 153; % K (Elevated pressure reference)
            Ox.calcK = @(t) 0.0015 + 0.00008 * t;

        case 'HTP' 
            % 90% H2O2 / 10% H2O by mass approximation
            Ox.MOLECULAR_MASS = 32.65; 
            Ox.HEAT_OF_FORMATION = -197300; % Liquid at 298 K (J/mol)
            Ox.Atoms.H = 2.12; 
            Ox.Atoms.O = 1.94;
            
            Ox.RHO = 1400.0; % kg/m3
            Ox.SURFACE_TENSION = 0.070; % N/m
            
            Ox.QS = 1530e3; % J/kg
            Ox.TEMP_BOIL = 414; % K
            Ox.calcK = @(t) 0.015 + 0.00006 * t;

        case {'N2O', 'NOX', 'NITROUS'}
            Ox.MOLECULAR_MASS = 44.013; 
            Ox.HEAT_OF_FORMATION = 65000; % Liquid N2O (Gas is +82050e3 minus ~17050e3 latent heat) (J/mol)
            Ox.Atoms.N = 2; 
            Ox.Atoms.O = 1;
            
            % Liquid properties at ~275 K (2 deg C, ~32 bar)
            % Note: If running warm at 293 K (20 C), RHO drops to ~785 and SURFACE_TENSION to ~0.0018
            Ox.RHO = 900.0; % kg/m3 
            Ox.SURFACE_TENSION = 0.0050; % N/m
            
            % Droplet vaporization properties
            Ox.QS = 250e3; % J/kg (376e3 at 1 bar, ~250e3 at 20 bar chamber pressure)
            Ox.TEMP_BOIL = 256.0; % K (Saturation temp at 20 bar; 184.7 K at 1 bar)
            Ox.calcK = @(t) 0.002 + 0.00006 * t; % Vapour thermal conductivity W/(m*K)
            
        otherwise
            error('loadPropData:UnknownOxidizer', 'Unknown oxidizer: %s', OX_NAME);
    end

    PROPELLANTS = [OX_NAME, '/', FUEL_NAME];
    
    switch PROPELLANTS
        case {'LOX/RP-1', 'LOX/METHANE', 'LOX/CH4', 'LOX/IPA'}
            % Hydrocarbon combustion species
            GAS_SPECIES = {'CO2','CO','H2O','H2','O2','OH','O','H', 'CH4','C'};
            CONDENSED_SPECIES = {'C_gr_'};
            ELEMENT_LIST = {'C','H','O'};
            
        case 'LOX/LH2'
            % Hydrolox combustion species (No Carbon)
            GAS_SPECIES = {'H2O','H2','O2','OH','O','H'};
            CONDENSED_SPECIES = {};
            ELEMENT_LIST = {'H','O'};

        case {'NOX/IPA', 'N2O/IPA', 'IPA/NOX', 'IPA/N2O'}
            % Hydrocarbon / Alcohol + Nitrous Oxide combustion species
            GAS_SPECIES = {'CO2','CO','H2O','H2','O2','OH','O','H','N2','NO','N'};
            CONDENSED_SPECIES = {'C_gr_'};
            ELEMENT_LIST = {'C','H','O','N'};
            
        otherwise
            error('loadPropData:UnknownPropCombo', ...
                'Species list for combination %s is not defined in the library.', PROPELLANTS);
    end
    
    SPECIES_LIST = [GAS_SPECIES, CONDENSED_SPECIES];
    NUM_SPECIES = numel(SPECIES_LIST);
    NUM_GAS = numel(GAS_SPECIES);
    NUM_ELEM = numel(ELEMENT_LIST);
    CONDENSED_IDX = (NUM_GAS+1) : NUM_SPECIES;
    
    % Load databases
    scriptDir = fileparts(mfilename('fullpath'));
    thermoMatFile = fullfile(scriptDir, 'thermo.mat');
    transMatFile  = fullfile(scriptDir, 'trans.mat');
    
    % Thermo database
    thermoInpFile = fullfile(scriptDir, 'thermo.inp');
    if ~isfile(thermoMatFile)
        fprintf('"%s" not found. Compiling from thermo.inp...\n', thermoMatFile);
        if ~isfile(thermoInpFile)
            error('Missing %s! Cannot compile thermodynamic database.', thermoInpFile);
        end
        compileThermo(thermoInpFile, thermoMatFile);
    end
    
    % Transport database
    transInpFile = fullfile(scriptDir, 'trans.inp');
    if ~isfile(transMatFile)
        fprintf('"%s" not found. Compiling from trans.inp...\n', transMatFile);
        if ~isfile(transInpFile)
            error('Missing %s! Cannot compile transport database.', transInpFile);
        end
        compileTransport(transInpFile, transMatFile);
    end
    
    ThermoData = load(thermoMatFile);
    TransData  = load(transMatFile);

    % Extract data
    
    MAX_T_INT = 3;
    A_IJ = zeros(NUM_ELEM, NUM_SPECIES);
    MOLAR_MASSES = zeros(NUM_SPECIES, 1);
    A_NASA = zeros(NUM_SPECIES * MAX_T_INT, 7);
    B_NASA = zeros(NUM_SPECIES * MAX_T_INT, 2);
    T_HI = inf(NUM_SPECIES, MAX_T_INT);
    
    MAX_V_INT = 3;
    V_COEFFS = zeros(NUM_SPECIES * MAX_V_INT, 4);
    V_RANGES_MAX = inf(NUM_SPECIES, MAX_V_INT);
    C_COEFFS = zeros(NUM_SPECIES * MAX_V_INT, 4);
    C_RANGES_MAX = inf(NUM_SPECIES, MAX_V_INT);

    for k = 1:NUM_SPECIES

        SPECIES = SPECIES_LIST{k};

        % Thermo
        
        if ~isfield(ThermoData, SPECIES)
            error('Species %s is missing from thermo.mat', SPECIES);
        end
        SpeciesThermo = ThermoData.(SPECIES);
        
        MOLAR_MASSES(k) = SpeciesThermo.MOLAR_MASS;
        
        for e = 1:NUM_ELEM
            sym = ELEMENT_LIST{e};
            if isfield(SpeciesThermo.ATOMS, sym)
                A_IJ(e, k) = SpeciesThermo.ATOMS.(sym);
            end
        end
        
        nT = size(SpeciesThermo.RANGES, 1);
        rT = (k-1) * MAX_T_INT + (1:nT);
        A_NASA(rT, :) = SpeciesThermo.A;
        B_NASA(rT, :) = SpeciesThermo.B;
        T_HI(k, 1:nT) = SpeciesThermo.RANGES(:, 2).';
        
        % Transport

        if ismember(k, CONDENSED_IDX) || ~isfield(TransData, SPECIES)
            continue;
        end
        
        SpeciesTrans = TransData.(SPECIES);

        nV = size(SpeciesTrans.V_RANGES, 1);
        rV = (k-1) * MAX_V_INT + (1:nV);
        V_COEFFS(rV, :) = SpeciesTrans.V_COEFFS;
        V_RANGES_MAX(k, 1:nV) = SpeciesTrans.V_RANGES(:, 2).';
        
        nC = size(SpeciesTrans.C_RANGES, 1);
        rC = (k-1) * MAX_V_INT + (1:nC);
        C_COEFFS(rC, :) = SpeciesTrans.C_COEFFS;
        C_RANGES_MAX(k, 1:nC) = SpeciesTrans.C_RANGES(:, 2).';
    end

    PropData = struct(...
        'FUEL', FUEL_NAME, 'OXIDISER', OX_NAME, ...
        'Fuel', Fuel, 'Ox', Ox, ...
        'SPECIES_NAMES', {SPECIES_LIST}, ...
        'ELEMENTS', {ELEMENT_LIST}, ...
        'A_IJ', A_IJ, 'bCirc', 0, ...
        'A_NASA', A_NASA, 'B_NASA', B_NASA, 'T_HI', T_HI, ...
        'M', MOLAR_MASSES, 'R', 8314.4626, ...
        'NUM_SPECIES', NUM_SPECIES, 'CONDENSED_LIST', CONDENSED_IDX, ...
        'V_COEFFS', V_COEFFS, 'V_RANGES_MAX', V_RANGES_MAX, ...
        'C_COEFFS', C_COEFFS, 'C_RANGES_MAX', C_RANGES_MAX);

end