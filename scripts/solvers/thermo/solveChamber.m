function [T_CHAMBER, N_J_CHAMBER, Propellant_Chamber] = solveChamber(P_CHAMBER, NUM_SPECIES, H_REACTANTS, GAS_DATA)
    
    % Inital guesses
    T_guess = 3800; % K
    N_guess = repmat(0.1/NUM_SPECIES, NUM_SPECIES, 1);

    if isfield(GAS_DATA,'CONDENSED_LIST')
        N_guess(GAS_DATA.CONDENSED_LIST) = 0;
    end

    [T, N_J, Propellant_Chamber] = equilibriumSolve('HP', H_REACTANTS, P_CHAMBER, GAS_DATA, T_guess, N_guess, 0);
    
    T_CHAMBER = T;
    N_J_CHAMBER = N_J;

end