function LoadParams = DefaultLoadParams
%Define structure with fields corresponding to loading parameters
LoadParams.animal = 'Mxxx';
LoadParams.series = '000';
LoadParams.exp = 101;
LoadParams.Shanknum = [0:7 0];
LoadParams.Shanksuffix = {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
LoadParams.Nthetaphsbins = 0;
LoadParams.Samplerate = 60;
LoadParams.LoadSmthTime = 150;
end