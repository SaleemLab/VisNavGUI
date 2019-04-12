function AnalysisOptions = DefaultAnalysisOptions
%Define structure with fields corresponding to loading parameters
AnalysisOptions.SpeedThreshold = 5;
LoadParams.iseries = '000';
LoadParams.explist = 101;
LoadParams.Shanknum = [0:7 0];
LoadParams.Shanksuffix = {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
LoadParams.Nthetaphsbins = 0;
LoadParams.samplerate = 60;
end