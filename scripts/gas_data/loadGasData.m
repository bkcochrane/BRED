function GasData = loadGasData(fuelName, oxName)

    if nargin < 1, fuelName = 'RP-1'; end
    if nargin < 2, oxName = 'LOX'; end
    
    propCombo = upper([oxName, '/', fuelName]);
    
    switch propCombo
        case 'LOX/RP-1'
            GAS_SPECIES = {'CO2','CO','H2O','H2','O2','OH','O','H'};
            CONDENSED_SPECIES = {'C_gr_'};
            ELEMENT_LIST = {'C','H','O'};
        case {'LOX/METHANE', 'LOX/CH4'}
            GAS_SPECIES = {'CO2','CO','H2O','H2','O2','OH','O','H'};
            CONDENSED_SPECIES = {'C_gr_'};
            ELEMENT_LIST = {'C','H','O'};
        case 'LOX/LH2'
            GAS_SPECIES = {'H2O','H2','O2','OH','O','H'};
            CONDENSED_SPECIES = {};
            ELEMENT_LIST = {'H','O'};
        otherwise
            error('loadGasData:UnknownPropellant', ...
                'Propellant combination %s is not defined in the library.', propCombo);
    end
    
    SPECIES_LIST = [GAS_SPECIES, CONDENSED_SPECIES];
    NUM_SPECIES = numel(SPECIES_LIST);
    NUM_GAS = numel(GAS_SPECIES);
    NUM_ELEM = numel(ELEMENT_LIST);
    CONDENSED_IDX = (NUM_GAS+1) : NUM_SPECIES;
    
    % Load databases
    script_dir = fileparts(mfilename('fullpath'));
    thermoMatFile = fullfile(script_dir, 'thermo.mat');
    transMatFile = fullfile(script_dir, 'trans.mat');
    
    % Thermo database
    if ~isfile(thermoMatFile)
        fprintf('"%s" not found. Compiling from thermo.inp...\n', thermoMatFile);
        if ~isfile('thermo.inp')
            error('Missing thermo.inp! Cannot compile the thermodynamic database.');
        end
        compileThermo('thermo.inp', thermoMatFile);
    end
    
    % Transport database
    if ~isfile(transMatFile)
        fprintf('"%s" not found. Compiling from trans.inp...\n', transMatFile);
        if ~isfile('trans.inp')
            error('Missing trans.inp! Cannot compile the transport database.');
        end
        compileTransport('trans.inp', transMatFile);
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
            error('Species %s is missing from thermo_complete_database.mat', SPECIES);
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

    GasData = struct(...
        'FUEL', fuelName, 'OXIDIZER', oxName, ...
        'SPECIES_NAMES', {SPECIES_LIST}, ...
        'ELEMENTS', {ELEMENT_LIST}, ...
        'A_IJ', A_IJ, 'B_I_CIRC', 0, ...
        'A_NASA', A_NASA, 'B_NASA', B_NASA, 'T_HI', T_HI, ...
        'M', MOLAR_MASSES, 'R', 8314.4626, ...
        'NUM_SPECIES', NUM_SPECIES, 'CONDENSED_LIST', CONDENSED_IDX, ...
        'V_COEFFS', V_COEFFS, 'V_RANGES_MAX', V_RANGES_MAX, ...
        'C_COEFFS', C_COEFFS, 'C_RANGES_MAX', C_RANGES_MAX);

end