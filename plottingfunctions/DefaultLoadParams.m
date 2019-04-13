function LoadParams = DefaultLoadParams
%Define structure with fields corresponding to loading parameters
LoadParams.animal = 'Mxxx';
LoadParams.series = '000';
LoadParams.exp = 101;
LoadParams.Channels{1} = 0:7;
LoadParams.Channels{2} = NaN;
LoadParams.Nthetaphsbins = 0;
LoadParams.Samplerate = 60;
LoadParams.LoadSmthTime = 150;
end