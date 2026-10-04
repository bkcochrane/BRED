function compileTrans(FILENAME, OUTPUT_MAT_FILE)

    if nargin < 1 || isempty(FILENAME)
        FILENAME = 'trans.inp';
    end
    if nargin < 2 || isempty(OUTPUT_MAT_FILE)
        script_dir = fileparts(mfilename('fullpath'));
        OUTPUT_MAT_FILE = fullfile(script_dir, 'trans.mat');
    end

    lines = readLines(FILENAME);

    TransData = struct();
    SPECIES_NAMES = {};

    i = 1;
    while i <= numel(lines)
        L = pad80(lines{i});

        % Skip comments, blank lines, data lines (which start with ' V' or ' C'), or the 'transport property coefficients' header block.
        if isempty(strtrim(L)) || L(1) == '!' || L(1) == '#' || L(1) == ' ' || strncmpi(L, 'transport', 9)
            i = i + 1;
            continue;
        end

        % Species names and binary pair names exist in columns 1 through 34.
        rawName = strtrim(L(1:34));
        
        % Check for number of temperature ranges
        tokens = regexp(L, 'V(\d+)C(\d+)', 'tokens');
        if isempty(tokens)
            i = i + 1;
            continue; 
        end
        
        NUM_V = str2double(tokens{1}{1});
        NUM_C = str2double(tokens{1}{2});

        % Convert to a valid MATLAB variable name
        validName = matlab.lang.makeValidName(regexprep(rawName, '\s+', '_'));

        % Ensure no overwriting if the file has duplicate names
        counter = 1;
        baseName = validName;
        while isfield(TransData, validName)
            counter = counter + 1;
            validName = sprintf('%s_%d', baseName, counter);
        end

        SPECIES_NAMES{end+1} = validName;

        S = struct();
        S.RAW_NAME = rawName;
        S.V_RANGES = zeros(NUM_V, 2);
        S.V_COEFFS = zeros(NUM_V, 4);
        S.C_RANGES = zeros(NUM_C, 2);
        S.C_COEFFS = zeros(NUM_C, 4);

        r = i + 1;
        
        % Read Viscosity coefficient lines
        for t = 1:NUM_V
            VL = pad80(lines{r});
            S.V_RANGES(t, :) = [fld(VL(3:10)), fld(VL(11:18))];
            S.V_COEFFS(t, :) = [fld(VL(19:35)), fld(VL(36:50)), fld(VL(51:65)), fld(VL(66:80))];
            r = r + 1;
        end

        % Read Conductivity coefficient lines
        for t = 1:NUM_C
            CL = pad80(lines{r});
            S.C_RANGES(t, :) = [fld(CL(3:10)), fld(CL(11:18))];
            S.C_COEFFS(t, :) = [fld(CL(19:35)), fld(CL(36:50)), fld(CL(51:65)), fld(CL(66:80))];
            r = r + 1;
        end

        TransData.(validName) = S;
        i = r; 
    end

    ExportData = TransData; 
    ExportData.SPECIES_NAMES_FULL = SPECIES_NAMES;

    save(OUTPUT_MAT_FILE, '-struct', 'ExportData');
    fprintf('\nSUCCESS! Parsed %d species and binary pairs.\n', numel(SPECIES_NAMES));
    fprintf('Saved all data to: %s\n', OUTPUT_MAT_FILE);

end

% ------------------------------------------------------------------------
% Helper Functions
% ------------------------------------------------------------------------
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
    
    % Standardize Fortran D and E notations
    s = strrep(strrep(s, 'D', 'E'), 'd', 'e');
    
    % NASA CEA format uses 'E 00' instead of 'E+00' for positive exponents.
    s = regexprep(s, '([Ee])\s+', '$1+');
    
    v = str2double(s);

end