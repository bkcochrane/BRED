% ------------------------------------------------------------------------
% Exit solve - maximise mass flux - G = rho * v (P)
% ------------------------------------------------------------------------
function [pExit, tExit, njExit, vExit, areaExit] = solveExit(h0, s0R, pExit, PropData, tempThroat, njThroat, mdot)
    
    [tExit, njExit, PropellantExit] = equilibriumSolve('SP', s0R, pExit, PropData, tempThroat, njThroat, 0);
    vExit = sqrt(2 * (h0 - PropellantExit.h));
    areaExit = mdot / (PropellantExit.rho * vExit);

end