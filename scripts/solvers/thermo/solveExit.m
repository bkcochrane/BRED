% ------------------------------------------------------------------------
% Exit solve - maximise mass flux - G = rho * v (P)
% ------------------------------------------------------------------------
function [P_EXIT, T_EXIT, N_J_EXIT, V_EXIT, AREA_EXIT] = solveExit(H0, S0_R, P_EXIT, GAS_DATA, T_THROAT, N_J_THROAT, MDOT)
    
    [T_EXIT, N_J_EXIT, Propellant_Exit] = equilibriumSolve('SP', S0_R, P_EXIT, GAS_DATA, T_THROAT, N_J_THROAT, 0);
    V_EXIT = sqrt(2 * (H0 - Propellant_Exit.H));
    AREA_EXIT = MDOT / (Propellant_Exit.RHO * V_EXIT);

end