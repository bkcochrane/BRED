function compileThermo(FILENAME, OUTPUT_MAT_FILE)
    if nargin < 1 || isempty(FILENAME)
        FILENAME = 'thermo.inp';
    end
    if nargin < 2 || isempty(OUTPUT_MAT_FILE)
        script_dir = fileparts(mfilename('fullpath'));
        OUTPUT_MAT_FILE = fullfile(script_dir, 'thermo.mat');
    end
    lines = readLines(FILENAME);
    ThermoData = struct();
    SPECIES_NAMES = {};
    ELEMENT_LIST = {};
    all_atoms = {};
    i = 1;
    while i <= numel(lines) - 2
        L = pad80(lines{i});
        trimmedL = strtrim(L);
        
        % Record headers begin in column 1. Skip indented lines, comments, and END cards.
        if isempty(trimmedL) || L(1) == ' ' || L(1) == '!' || L(1) == '#' || ...
           startsWith(upper(trimmedL), 'END') || strcmpi(trimmedL, 'thermo')
            i = i + 1;
            continue;
        end
        rawName = strtrim(L(1:18));
        if isempty(rawName)
            i = i + 1;
            continue;
        end
        H = pad80(lines{i+1});
        nInt = str2double(H(1:2));
        
        % Check if nInt is valid (0 = single-temp reactant like H2(L), 1..5 = standard)
        if isnan(nInt) || nInt < 0 || nInt > 5
            i = i + 1;
            continue;
        end
        
        % Guard against non-species lines starting in col 1 (e.g. temperature header line)
        if fld(H(53:65)) <= 0
            i = i + 1;
            continue;
        end

        % Convert to a valid MATLAB variable name (e.g. C(gr) -> C_gr_)
        validName = matlab.lang.makeValidName(rawName);
        
        % Ensure no overwriting if thermo.inp has multiple phases with the same name
        counter = 1;
        baseName = validName;
        while isfield(ThermoData, validName)
            counter = counter + 1;
            validName = sprintf('%s_%d', baseName, counter);
        end
        
        SPECIES_NAMES{end+1} = validName;
        S = struct();
        S.RAW_NAME = rawName;
        S.MOLAR_MASS = fld(H(53:65));
        S.HEAT_OF_FORMATION = fld(H(66:80));
        
        % Find and record atoms (preserving exact original element order & negative ions)
        local_atoms = struct();
        for p = 0:4
            c = 11 + 8*p;
            sym = strtrim(H(c:c+1));
            
            % Only skip empty or purely numeric strings like '0' or '00'
            if isempty(sym) || ~all(isstrprop(sym, 'alpha'))
                continue;
            end
            
            if ~ismember(sym, ELEMENT_LIST)
                ELEMENT_LIST{end+1} = sym;
            end
            
            val = fld(H(c+2:c+7));
            local_atoms.(sym) = val;
        end
        all_atoms{end+1} = local_atoms;
        S.ATOMS = local_atoms;
        
        % Get a and b coefficients, based on the standard Thermo.inp file formatting
        r = i + 2;
        if nInt == 0
            % Single-point cryogenic liquid reactant (e.g. H2(L), HNO3(L))
            TL = pad80(lines{r});
            T_ref = fld(TL(1:11));
            S.RANGES = [T_ref, T_ref];
            S.A = zeros(1, 7);
            % Store H_form / R in b1 so standard NASA polynomial evaluation yields exact enthalpy
            S.B = [S.HEAT_OF_FORMATION / 8.314462618, 0];
            r = r + 1;
        else
            S.RANGES = zeros(nInt, 2);
            S.A = zeros(nInt, 7);
            S.B = zeros(nInt, 2);
            for t = 1:nInt
                TL = pad80(lines{r});
                C1 = pad80(lines{r+1});
                C2 = pad80(lines{r+2});
                S.RANGES(t,:) = [fld(TL(1:11)),  fld(TL(12:22))];
                S.A(t,1:5) = [fld(C1(1:16)),  fld(C1(17:32)), fld(C1(33:48)), ...
                                 fld(C1(49:64)), fld(C1(65:80))];
                S.A(t,6:7) = [fld(C2(1:16)),  fld(C2(17:32))];
                S.B(t,:) = [fld(C2(49:64)), fld(C2(65:80))];
                r = r + 3;
            end
        end
        ThermoData.(validName) = S;
        i = r;
    end
    
    % Create A_ij matrix on all species and elements
    NUM_SPECIES = numel(SPECIES_NAMES);
    NUM_ELEMENTS = numel(ELEMENT_LIST);
    A_IJ = zeros(NUM_ELEMENTS, NUM_SPECIES);
    
    for k = 1:NUM_SPECIES
        local_atoms = all_atoms{k};
        fields = fieldnames(local_atoms);
        for f = 1:numel(fields)
            sym = fields{f};
            idx = strcmp(ELEMENT_LIST, sym);
            A_IJ(idx, k) = local_atoms.(sym);
        end
    end
    ExportData = ThermoData; 
    ExportData.A_IJ_MATRIX_FULL = A_IJ;
    ExportData.SPECIES_NAMES_FULL = SPECIES_NAMES;
    ExportData.ELEMENT_LIST_FULL = ELEMENT_LIST;
    
    save(OUTPUT_MAT_FILE, '-struct', 'ExportData');
    fprintf('\nSUCCESS! Parsed %d species containing %d unique elements.\n', NUM_SPECIES, NUM_ELEMENTS);
    fprintf('Saved all data to: %s\n', OUTPUT_MAT_FILE);
end

function lines = readLines(filename)
    fid = fopen(filename, 'r');
    if fid < 0, error('Cannot open %s', filename); end
    cleanup = onCleanup(@() fclose(fid));
    raw   = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = raw{1};
    keep = true(numel(lines), 1);
    for i = 1:numel(lines)
        s = lines{i};
        if isempty(strtrim(s)) || s(1) == '!' || s(1) == '#', keep(i) = false; end
    end
    lines = lines(keep);
end
function s = pad80(s)
    if numel(s) < 80, s = [s, blanks(80 - numel(s))]; end
end
function v = fld(s)
    s = strtrim(s);
    if isempty(s), v = 0; return; end
    v = str2double(strrep(strrep(s, 'D', 'E'), 'd', 'e'));
end