function Engine = sizeInjector(Engine, PropData)

    RHO_FUEL = PropData.Fuel.RHO; % kg/m3
    RHO_OX = PropData.Ox.RHO; % kg/m3
    Cd = Engine.Injector.Cd;
    t = Engine.Injector.PLATE_THICKNESS; % m
    mdotFuel = Engine.Fuel.mdot; % kg/s
    mdotOx = Engine.Ox.mdot; % kg/s
    DIAMETER_MIN = Engine.Injector.PORT_MIN_DIAMETER; % m

    P_CHAMBER = Engine.Chamber.P; % Pa
    dP_FUEL = Engine.Fuel.P - P_CHAMBER; % Pa
    dP_OX = Engine.Ox.P - P_CHAMBER; % Pa

    areaFuel = mdotFuel / (Cd * sqrt(2 * RHO_FUEL * dP_FUEL)); % m2
    areaOx = mdotOx / (Cd * sqrt(2 * RHO_OX * dP_OX)); % m2

    numOxPerFuel = 2; % Start with a standard triplet (2 ox, 1 fuel)
    MAX_DIAMETER_RATIO = 1.5; % Define "too large" (e.g., Ox max 50% larger than fuel)

    while numOxPerFuel < 6 % Safety cap at 5 ox ports per fuel port
        % Calculate area and diameter ratio for the current number of ports
        areaRatioPorts = (1.0 / numOxPerFuel) * (areaOx / areaFuel);
        diameterRatioPorts = sqrt(areaRatioPorts);
        
        if diameterRatioPorts <= MAX_DIAMETER_RATIO
            break; % The size is acceptable, exit the loop
        end
        
        % If the ox diameter is still too large, add another ox port and try again
        numOxPerFuel = numOxPerFuel + 1; 
    end

    % Apply the calculated diameter ratio
    if diameterRatioPorts > 1
        diameterFuel = DIAMETER_MIN;
        diameterOx = DIAMETER_MIN * diameterRatioPorts;
    else
        diameterFuel = DIAMETER_MIN / diameterRatioPorts;
        diameterOx = DIAMETER_MIN;
    end

    numPortsFuel = floor(areaFuel / (pi / 4 * diameterFuel^2));
    numPortsOx = numOxPerFuel * numPortsFuel;

    vFuel = mdotFuel / (RHO_FUEL * areaFuel); % m/s
    vOx = mdotOx / (RHO_OX * areaOx); % m/s

    momentumRatio = (mdotFuel * vFuel) / (mdotOx * vOx);
    
    minGammaAngle = deg2rad(45);
    maxGammaAngle = deg2rad(60);
    
    if momentumRatio > 1
        gammaAngle = minGammaAngle; 
    else
        gammaAngleIdeal = acos(momentumRatio);
        gammaAngle = max(minGammaAngle, min(maxGammaAngle, gammaAngleIdeal));
    end

    JET_DIAMETERS = 5; % port diameters
    lengthJet = JET_DIAMETERS * diameterOx; % m - along the ox jet
    lengthImpingement = lengthJet * cos(gammaAngle); % m - axial from the face
    spacing = lengthImpingement * tan(gammaAngle);
    
    web = spacing - 0.5*(diameterFuel + diameterOx);
    if web < 0.001
        fprintf('[Warning] Land between fuel and ox ports is %.2f mm - too thin to drill\n', web*1000);
    end

    momentumFuel = mdotFuel * vFuel;
    momentumOx = mdotOx * vOx * cos(gammaAngle);
    ratioMomentum = momentumFuel / momentumOx;

    stiffnessFuel = dP_FUEL / P_CHAMBER;
    stiffnessOx = dP_OX / P_CHAMBER;

    LDFuel = t / diameterFuel;
    LDOx = t / diameterOx / cos(gammaAngle);

    if ratioMomentum > 1
        fprintf('[Warning] Fuel momentum greater than oxidiser momentum - mixing efficiency will drop.\n')
    end
    if ratioMomentum < 0.67
        fprintf('[Warning] Oxidiser momentum far greater than fuel momentum - mixing efficiency will drop.\n')
    end
    if stiffnessFuel < 0.15 || stiffnessOx < 0.15
        fprintf('[Warning] Injector stiffness below 0.15 - combustion instability risk.\n');
    end
    if LDFuel < 3 || LDOx < 3
        fprintf('[Warning] Orifice L/D below 3 - discharge coefficient will be unstable.\n');
    end

    Engine.Injector.impingementAngle = gammaAngle;
    Engine.Injector.lengthImpinge = lengthJet;
    Engine.Injector.momentumRatio = ratioMomentum;
    Engine.Injector.numPortSets = numPortsFuel;
    
    Engine.Injector.numPortsFuel = numPortsFuel;
    Engine.Injector.diameterPortsFuel = diameterFuel;
    Engine.Injector.areaPortsFuel = areaFuel;
    
    Engine.Injector.numPortsOx = numPortsOx;
    Engine.Injector.diameterPortsOx = diameterOx;
    Engine.Injector.areaPortsOx = areaOx;

    Engine.Injector.innerRingDiameter = 2 * spacing + diameterFuel;
    Engine.Injector.middleRingDiameter = 2 * spacing + Engine.Injector.innerRingDiameter;
    Engine.Injector.outerRingDiameter = 2 * spacing + Engine.Injector.middleRingDiameter;

    if Engine.Injector.outerRingDiameter > (Engine.Chamber.diameter - 2 * spacing)
        fprintf('[Warning] Injector outer ring exceeds chamber diameter or is dangerously close to the wall.\n');
    end    

    Engine.Fuel.v = vFuel;
    Engine.Fuel.stiffness = stiffnessFuel;
    Engine.Fuel.dP = dP_FUEL;
    
    Engine.Ox.v = vOx;
    Engine.Ox.stiffness = stiffnessOx;
    Engine.Ox.dP = dP_OX;
    
    % Droplet sizing - Anderson 1995
    rho = [PropData.Fuel.RHO, PropData.Ox.RHO];
    velocity = [vFuel, vOx];
    d0 = [diameterFuel, diameterOx];
    sigma = [PropData.Fuel.SURFACE_TENSION, PropData.Ox.SURFACE_TENSION];
    rhoGas = Engine.Chamber.rho;
    rhoLiquid = [PropData.Fuel.RHO, PropData.Ox.RHO];

    weber = (rho .* velocity.^2 .* d0) ./ sigma;
    fTheta = (1 - cos(gammaAngle)).^2 / sin(gammaAngle).^3;
    s = rhoGas ./ rhoLiquid;

    dD = (2.62 / 12.^(1/3)) .* s.^(-1/6) .* (weber .* fTheta).^(-1/3) .* d0;
    
    Engine.Fuel.d0 = dD(1);
    Engine.Ox.d0 = dD(2);

end