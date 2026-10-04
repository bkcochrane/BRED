function T_SAT = waterTSat(P)

    % Steam-table saturation
    P_BARs = [1 5 10 20 30 40 50 70 100 150 200 220.6]; % bar
    T_SATs = [373.1 425.0 453.0 485.5 507.0 523.5 537.1 559.0 584.1 615.3 638.9 647.1]; % K
    T_SAT = interp1(P_BARs, T_SATs, P/1e5, 'pchip');

end