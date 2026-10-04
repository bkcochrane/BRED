function P = solvePressure(S0_R, H0, Engine, GAS_DATA, T, N_J, MDOT, AREA, P_REF, branch, IS_FROZEN)

    G_TARGET = MDOT / AREA;
    flux_err = @(P) findMassFlux(P, S0_R, H0, GAS_DATA, T, N_J, IS_FROZEN) - G_TARGET;

    if strcmp(branch,'supersonic')
        bracket = [10, Engine.P_THROAT]; % P_REF = P_THROAT, supersonic root below it
    else
        bracket = [Engine.P_THROAT, Engine.P_CHAMBER];  % subsonic root above throat
    end

    try
        P = fzero(flux_err, bracket);
    catch
        P = fminbnd(@(P) abs(flux_err(P)), bracket(1), bracket(2));
    end

end