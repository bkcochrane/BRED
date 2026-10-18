function [tempChamber, njChamber, PropellantChamber] = solveChamber(P_CHAMBER, NUM_SPECIES, hReactants, PropData)
    
    % Inital guesses
    guessTemp = 3800; % K
    guessN = repmat(0.1/NUM_SPECIES, NUM_SPECIES, 1);

    if isfield(PropData,'CONDENSED_LIST')
        guessN(PropData.CONDENSED_LIST) = 0;
    end

    [temp, nj, PropellantChamber] = equilibriumSolve('HP', hReactants, P_CHAMBER, PropData, guessTemp, guessN, 1);
    
    tempChamber = temp;
    njChamber = nj;

end