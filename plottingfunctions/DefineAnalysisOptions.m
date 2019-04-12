function AnalysisOptions = DefineAnalysisOptions
%Define structure with fields corresponding to options for analyses to run
%when loading a sessions and directories where to save/load preprocessed
%files
AnalysisOptions = Tstructure('AnalysisOptions');

% SpeedThreshold
% PlotVar1DMaps.SmthTimeSpeed
% PlotVar1DMaps.SmthTimeWindow, Xbinsize, PlotVar1DMaps.SmthSpatialWindow,delay, EXP.data.es.CircularMaze
% Phsbinsize, SmthPhsWindow
% 
% AnalysisOptions.addprop('animal');
% AnalysisOptions.animal = 'Mxxx';
end