function [G, T, N_J] = findMassFlux(P, S0_R, H0, GAS_DATA, T, N_J, IS_FROZEN)

    if nargin < 7;IS_FROZEN = false;end

    [T, N_J, Propellant] = equilibriumSolve('SP', S0_R, P, GAS_DATA, T, N_J, 0, IS_FROZEN);

    V  = sqrt(max(2 * (H0 - Propellant.H), 0));
    G  = Propellant.RHO * V;

end