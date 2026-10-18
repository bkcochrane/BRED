function [massFlux, temp, nj] = findMassFlux(p, s0R, h0, PropData, temp, nj, isFrozen)
    
    if nargin < 7; isFrozen = false; end

    [temp, nj, Propellant] = equilibriumSolve('SP', s0R, p, PropData, temp, nj, 0, isFrozen);

    v = sqrt(max(2 * (h0 - Propellant.h), 0));
    massFlux = Propellant.rho * v;

end