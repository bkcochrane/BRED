function lengthCombustion = findCombustionLength(Engine)

    % Length over which combustion is incomplete, taken as the distance a
    % droplet travels while vaporizing. Reuses the Ranz-Marshall vaporization
    % time from solveRealEngine so this is consistent with the L* sizing.
    % Bounded below by the impingement point (nothing burns before the jets
    % meet) and above by the chamber length (combustion must be complete by
    % the throat).

    vDroplet = Engine.Chamber.vDroplet;
    timeVapourise = Engine.Chamber.timeVapourise;
    lengthCombustion = vDroplet * timeVapourise;

    lengthCombustion = max(lengthCombustion, Engine.Injector.lengthImpinge);
    lengthCombustion = min(lengthCombustion, Engine.Geometry.lengthChamber);

end