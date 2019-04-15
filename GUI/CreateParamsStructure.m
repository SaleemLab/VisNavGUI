function P = CreateParamsStructure(P, savedParamsFile)
%create Parameter structure used troughout the GUI
if ~isempty(savedParamsFile)
    Psaved = load(savedParamsFile);
else
    Psaved = [];
end

P.addprop('DatabaseFilters');
P.DatabaseFilters = DefaultDatabaseFilters;
if ~isempty(Psaved)
    P.DatabaseFilters = updateFields(P.DatabaseFilters,Psaved.DatabaseFilters);
end

P.addprop('LoadParams');
P.LoadParams = DefaultLoadParams;
if ~isempty(Psaved)
    P.LoadParams = updateFields(P.LoadParams,Psaved.LoadParams);
end

P.addprop('AnalysisOptions');
P.AnalysisOptions = DefaultAnalysisOptions;
if ~isempty(Psaved)
    P.AnalysisOptions = updateFields(P.AnalysisOptions,Psaved.AnalysisOptions);
end

P.addprop('PlotParamsMaps');
P.PlotParamsMaps = DefaultPlotParamsMaps;
if ~isempty(Psaved)
    P.PlotParamsMaps = updateFields(P.PlotParamsMaps,Psaved.PlotParamsMaps);
end

P.addprop('PlotParamsDecoder');
P.PlotParamsDecoder = DefaultPlotParamsDecoder;
if ~isempty(Psaved)
    P.PlotParamsDecoder = updateFields(P.PlotParamsDecoder,Psaved.PlotParamsDecoder);
end

P.addprop('PlotParamsBehavior');
P.PlotParamsBehavior = DefaultPlotParamsBehavior;
if ~isempty(Psaved)
    P.PlotParamsBehavior = updateFields(P.PlotParamsBehavior,Psaved.PlotParamsBehavior);
end

P.addprop('DIRS');
P.DIRS = SetDirectories;
if ~isempty(Psaved)
    P.DIRS = updateFields(P.DIRS,Psaved.DIRS);
end

end