% -----------------------------------------------------------------------
% Throat solve - maximise mass flux - G = rho * v (P)
% ------------------------------------------------------------------------
function [P_THROAT, T_THROAT, N_J_THROAT, V_T, AREA_THROAT, C_STAR, MACH_T] = solveThroat(P_CHAMBER, H0, S0_R, T_CHAMBER, N_J_CHAMBER, GAS_DATA, MDOT)
        
    % Pressure sweep to find throat
    P_SWEEP = P_CHAMBER * (0.9 : -0.01 : 0.2);
    G  = zeros(size(P_SWEEP));
    
    T = T_CHAMBER;
    N_J = N_J_CHAMBER;
    
    for i = 1:numel(P_SWEEP)
        [G(i), T, N_J] = findMassFlux(P_SWEEP(i), S0_R, H0, GAS_DATA, T, N_J);
    end
    
    % Find location of peak on the coarse sweep
    [~,ip] = max(G);
    lo = P_SWEEP(min(ip+1,end)); % Lower pressure bound
    hi = P_SWEEP(max(ip-1,1));   % Higher pressure bound

    % Refine peak location by maximizing mass flux directly
    flux_err = @(P) -findMassFlux(P, S0_R, H0, GAS_DATA, T_CHAMBER, N_J_CHAMBER);
    P_THROAT = fminbnd(flux_err, lo, hi);
    
    % Get propellant properties at throat
    [T_THROAT, N_J_THROAT, Propellant_Throat] = equilibriumSolve('SP', S0_R, P_THROAT, GAS_DATA, T_CHAMBER, N_J_CHAMBER, 0);
    
    V_T = sqrt(max(2 * (H0 - Propellant_Throat.H), 0));
    AREA_THROAT = MDOT / (Propellant_Throat.RHO * V_T);
    C_STAR = P_CHAMBER * AREA_THROAT / MDOT;
    
    [MACH_T, ~] = findMach(P_THROAT, S0_R, H0, GAS_DATA, T_CHAMBER, N_J_CHAMBER);

end