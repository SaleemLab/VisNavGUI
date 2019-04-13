function VisNavGUI(savedParamsFile)
if nargin < 1
    savedParamsFile = [];
end

%Load or create Parameter structure used troughout the GUI
P = CreateParamsStructure(savedParamsFile);

%create and define fields of the data structure used troughout the GUI
EXP = CreateEXPStructure;

%create the GUI object
GUI = TMultigraph('VisNav');
GUI.addPage;
GUI.addPage;

GUI.RenamePage(P.PlotParamsMaps.Page,P.PlotParamsMaps.PageTitle);
GUI.Hdividepage(P.PlotParamsMaps.Page, 2, [1 9]);
SpatialMapsDialog = Tdialog(GUI.window{P.PlotParamsMaps.Page,P.PlotParamsMaps.DialogWindow});

GUI.RenamePage(P.PlotParamsDecoder.Page,P.PlotParamsDecoder.PageTitle);
GUI.Hdividepage(P.PlotParamsDecoder.Page, 2, [1 9]);
DecodingDialog = Tdialog(GUI.window{P.PlotParamsDecoder.Page,P.PlotParamsDecoder.DialogWindow});

GUI.RenamePage(P.PlotParamsBehavior.Page,P.PlotParamsBehavior.PageTitle);
GUI.Hdividepage(P.PlotParamsBehavior.Page, 2, [1 9]);
BehaviorDialog = Tdialog(GUI.window{P.PlotParamsBehavior.Page,P.PlotParamsBehavior.DialogWindow});


uimenu('Parent',GUI.FileMenu,'Label','load VR session',...
       'Callback', @(source,event)GUI_LoadVRMenu_Callback(source, event, GUI, EXP, P, SpatialMapsDialog, DecodingDialog, BehaviorDialog));
uimenu('Parent',GUI.FileMenu,'Label','save file',...
       'Callback',@(source,event)GUI_SaveMenu_Callback(source,event,EXP, P));
uimenu('Parent',GUI.FileMenu,'Label','load processed file',...
       'Callback',@(source,event)GUI_LoadProcessedMenu_Callback(source, event, GUI, EXP, P, SpatialMapsDialog, DecodingDialog, BehaviorDialog));

uimenu('Parent',GUI.AnalysisMenu,'Label','Analysis Options','Callback',@(source,event)GUI_AnalysisOptions_Callback(source,event,P));
uimenu('Parent',GUI.AnalysisMenu,'Label','Compute 1D maps','Callback',@(source,event)GUI_Compute1Dmaps_Callback(source,event,P));
uimenu('Parent',GUI.AnalysisMenu,'Label','Compute 2D maps','Callback',@(source,event)GUI_Compute2Dmaps_Callback(source,event,P));
uimenu('Parent',GUI.AnalysisMenu,'Label','Run Bayesian decoder','Callback',@(source,event)GUI_BayesDecoder_Callback(source,event,P));

uimenu('Parent',GUI.PlotMenu,'Label','Save Figure','Callback',@(source,event)GUI_SaveFigure_Callback(source,event));
uimenu('Parent',GUI.PlotMenu,'Label','Save Plot','Callback',@(source,event)GUI_SavePlot_Callback(source,event));
uimenu('Parent',GUI.PlotMenu,'Label','Add colorbar','Callback',@(source,event)GUI_Colorbar_Callback(source,event));
end

% % *************************Callback functions***************************** %

function GUI_LoadVRMenu_Callback(source, event, GUI, EXP, P, SpatialMapsDialog, DecodingDialog, BehaviorDialog)
%
%S = FAFF(arg);
% LoadandRun(EXP, animalname, iseries, iexplist, processedfiles);
if ~Fprocessed
end
LoadData_VisNav(EXP, P, S);
figstr = [num2str(LoadParams.LoadParams.animal) ' ' num2str(LoadParams.LoadParams.iseries)];
GUI.updateTitle(figstr);
InstallSpatialMapsDialog(GUI, EXP, P, SpatialMapsDialog);
InstallDecodingDialog(GUI, EXP, P, DecodingDialog);
InstallBehaviorDialog(GUI, EXP, P, BehaviorDialog);
end

function GUI_SaveMenu_Callback(source,event,EXP,P)
end

function GUI_LoadProcessedMenu_Callback(source, event, GUI, EXP, P, SpatialMapsDialog, DecodingDialog, BehaviorDialog)
dirname = '';
[filename, dirname] = uigetfile(dirname);

S = load([dirname filename]);
exp = TVRData;
exp.Copyobj(S.EXP);
CreateDataStructure(exp, EXP);

figstr = [num2str(EXP.animal) ' ' num2str(EXP.series)];
GUI.updateTitle(figstr);
InstallSpatialMapsDialog(GUI, EXP, P, SpatialMapsDialog);
% InstallDecodingDialog(GUI, EXP, P, DecodingDialog);
% InstallBehaviorDialog(GUI, EXP, P, BehaviorDialog);
end

function GUI_AnalysisOptions_Callback(source,event,P)
AnalysisOptionsDialog = Tdialog([],'Analysis Options',[0.3 0.1 0.4 0.8]);
end