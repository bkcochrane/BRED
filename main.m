clear all;
close all;
clc;

addpath(genpath(pwd));

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

Engine = engine_RAPTOR;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

PropData = loadPropData(Engine.Fuel.NAME, Engine.Ox.NAME);

% findOptimalRatio(Engine, PropData);
% parameterSweeps(Engine, PropData)

[Ideal, Engine] = createEngine(Engine, PropData);

printNozzleResults(Ideal, 'IDEAL');
printNozzleResults(Engine, 'REAL');
printInjector(Engine);
printThermalResults(Engine);

plotResults(Engine)
% plotMassPercentages(Engine, PropData)
% plot3DContour(Engine)
% writeCsv(Engine)

% plotThrottleCurve(Engine, PropData)