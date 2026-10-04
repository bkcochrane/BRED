clear variables;
close all;
clc;

addpath(genpath(pwd));

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

Engine.TARGET_THRUST = 1000; % N
Engine.EXIT_MODE = 'pressure';
Engine.AREA_RATIO_TARGET = 10;
Engine.EXIT_PR = 0.5;
Engine.P_AMBIENT = 101325; % Pa - Ambient pressure
Engine.P_EXIT = Engine.EXIT_PR * Engine.P_AMBIENT; % Pa
Engine.OF_RATIO = 2.4;
Engine.P_CHAMBER = 20e5; % Pa
Engine.L_STAR = 1.3; % m
Engine.CONTRACTION_RATIO = 10; % AREA_CHAMBER / AREA_THROAT
Engine.CONE_PERCENT = 90;
Engine.BETA = deg2rad(30); % radians - Converging angle
Engine.BLEND_RADIUS_THROAT_FACTOR = 1.5;
Engine.BLEND_RADIUS_CHAMBER_FACTOR = 2;
Engine.D0 = 68e-6;
Engine.ETA_MIXING = 0.98;
Engine.ETA_KINETIC = 0.96;
Engine.DISCHARGE_COEFF = 0.97;
Engine.T_COOLANT_IN = 300; % K
Engine.CP_COOLANT = 4200;
Engine.RHO_COOLANT = 998; % kg/m3
Engine.MU_COOLANT = 1.0e-3; % Pa s
Engine.K_COOLANT = 0.6; % W/(m K)
Engine.P_COOLANT = 10e5; % Pa
Engine.MDOT_COOLANT = 2; % kg/s
Engine.CHANNEL_WIDTH_MIN = 0.0025; % m
Engine.CHANNEL_HEIGHT_THROAT = 0.001; % m
Engine.CHANNEL_HEIGHT_MAX = 0.0045; % m
Engine.CHANNEL_ROUGHNESS = 15e-6; % m
Engine.LAND_WIDTH = 0.001; % m
Engine.WALL_THICKNESS = 0.001; % m
Engine.WALL_K = 100; % W/(m K)
Engine.N_STARTS = 0;
Engine.N_TURNS = 1;
Engine.RHO_FUEL = 810; % kg/m3
Engine.RHO_OX = 1140; % kg/m3
Engine.P_FUEL = 30e5; % Pa
Engine.P_OX = 25e5; % Pa
Engine.INJECTOR_DISCHARGE_COEFFICIENT = 0.8;
Engine.INJECTOR_PLATE_THICKNESS = 0.010; % m
Engine.INJECTOR_PORT_MIN_DIAMETER = 0.0015; % m
Engine.INJECTOR_IMPINGEMENT_ANGLE = deg2rad(45); % radians

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

GasData = loadGasData('RP-1', 'LOX');

% parameterSweep(Engine, GAS_DATA)
% EngineSweeps(Engine, GAS_DATA);

Engine = createEngine(Engine, GasData);
Engine = solveThermal(Engine, GasData);

plotResults(Engine)
% plot3DContour(Engine)
writeCsv(Engine)

% plotThrottleCurve(Engine, GAS_DATA)