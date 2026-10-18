% -----------------------------------------------------------------------
% Throat solve - maximise mass flux - G = rho * v (P)
% ------------------------------------------------------------------------
function [pThroat, tempThroat, njThroat, vThroat, areaThroat, cStar, machThroat] = solveThroat(P_CHAMBER, h0, s0R, tempChamber, njChamber, PropData, mdot)
        
    % Pressure sweep to find throat
    pSweep = P_CHAMBER * (0.9 : -0.01 : 0.2);
    massFlux = zeros(size(pSweep));
    
    temp = tempChamber;
    nj = njChamber;
    
    for i = 1:numel(pSweep)
        [massFlux(i), temp, nj] = findMassFlux(pSweep(i), s0R, h0, PropData, temp, nj);
    end
    
    % Find location of peak on the coarse sweep
    [~,ip] = max(massFlux);
    lo = pSweep(min(ip+1,end)); % Lower pressure bound
    hi = pSweep(max(ip-1,1));   % Higher pressure bound

    % Refine peak location by maximizing mass flux directly
    negMassFlux = @(p) -findMassFlux(p, s0R, h0, PropData, tempChamber, njChamber);
    pThroat = fminbnd(negMassFlux, lo, hi);
    
    % Get propellant properties at throat
    [tempThroat, njThroat, PropellantThroat] = equilibriumSolve('SP', s0R, pThroat, PropData, tempChamber, njChamber, 0);
    
    vThroat = sqrt(max(2 * (h0 - PropellantThroat.h), 0));
    areaThroat = mdot / (PropellantThroat.rho * vThroat);
    cStar = P_CHAMBER * areaThroat / mdot;
    
    [machThroat, ~] = findMach(pThroat, s0R, h0, PropData, tempChamber, njChamber);

end