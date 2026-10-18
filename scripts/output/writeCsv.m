function writeCsv(Engine)

    % Gas side
    xHot = Engine.Stream.x(:) * 1000; % mm
    rHot = Engine.Stream.r(:) * 1000;
    [xHot, iu] = unique(xHot);
    rHot = rHot(iu);
    s0 = [0; cumsum(hypot(diff(xHot), diff(rHot)))];
    if exist('csaps','file')
        xHot = csaps(s0, xHot, 0.999, s0);
        rHot = csaps(s0, rHot, 0.999, s0);
    end
    xHot = xHot(:);  rHot = rHot(:);

    % Channel dimensions
    t = Engine.Wall.THICKNESS * 1000; % mm
    heights = Engine.Cooling.heightChannelStation(iu) * 1000; % mm
    widths = Engine.Cooling.widthChannelStation(iu)  * 1000; % mm

    % Wall normal
    dX = gradient(xHot);  dR = gradient(rHot);  L = hypot(dX, dR);
    nX = -dR ./ L;
    nR = dX ./ L;

    % Cold side (channel floor) and outer (channel roof), both normal offsets
    xFloor = xHot + t .* nX;
    rFloor = rHot + t .* nR;
    xRoof = xHot + (t + heights) .* nX;
    rRoof = rHot + (t + heights) .* nR;
    
    % Helix
    if isfield(Engine.Cooling,'SPIRAL_RATE'), RATE = Engine.Cooling.SPIRAL_RATE / 1000; else, RATE = 0; end
    rMid = 0.5 * (rFloor + rRoof);
    theta = RATE * (xHot - min(xHot));

    % Resample for smooth curve
    numPoints = 300;
    s = [0; cumsum(hypot(diff(xHot), diff(rHot)))];
    sf = linspace(0, s(end), numPoints).';
    [sUnique, unique_idx] = unique(s);
    rs = @(v) interp1(sUnique, v(unique_idx), sf, 'pchip');

    xHot = rs(xHot);
    rHot = rs(rHot);
    xFloor = rs(xFloor);
    rFloor = rs(rFloor);
    rMid = rs(rMid);
    xRoof = rs(xRoof);
    rRoof = rs(rRoof);
    theta = rs(theta);
    z = zeros(numPoints, 1);

    writematrix([xHot, rHot, z], "CAD_curves/nozzle_hot_side.csv");
    writematrix([xFloor, rFloor, z], "CAD_curves/nozzle_channel_floor.csv");
    writematrix([xRoof, rRoof, z], "CAD_curves/nozzle_channel_roof.csv");
    
    % Helix for sweep path
    writematrix([xHot, 2 * rMid(1).*cos(theta), 2 * rMid(1).*sin(theta)], "CAD_curves/nozzle_helix.csv");
    writematrix([xHot, z, z], "CAD_curves/nozzle_centerline.csv")
    
    % Cross section for sweep
    w = min(widths);
    section = [
        xFloor(1), 2 * rMid(1), -w/2;
        xFloor(1), 0, -w/2;
        xFloor(1), 0, w/2;
        xFloor(1), 2 * rMid(1), w/2;
        xFloor(1), 2 * rMid(1), -w/2;
        ];

    writematrix(section, "CAD_curves/nozzle_channel_section.csv");

    fprintf('\nGeometry csv files written!\n');

end