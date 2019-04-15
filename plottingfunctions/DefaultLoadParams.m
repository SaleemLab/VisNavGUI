function LoadParams = DefaultLoadParams
% Define structure with fields corresponding to loading parameters.
% Field names are:
% - animal: Animal name
% - series: series number
% - exp: list of experiment numbers
% - Channels: n x 1 cell array of channels to load where n = number of probes
% - Nthetaphsbins: Number of bins per cycle (used for theta binning for
%   now)
% - Samplerate: Sampling rate at which data should be resampled
% - LoadSmthTime: Time smoothing window used for some variables when loading 
%   the data (e.g. smthBallSpeed)
% - SynchType: Type of signal used to synch recordings from different sources
% - SynchSignalRef: Orgin of the Synchsignal used as a reference for all 
%   types of data

%Animal name
LoadParams.animal = 'Mxxx';

%Series number
LoadParams.series = '000';

%experiment list
LoadParams.exp = 101;

%Active channels as n x 1 cell array where n = number of probes
LoadParams.Channels{1} = 0:7;
LoadParams.Channels{2} = NaN;

%Number of bins per cycle (could be used for theta binning or something
%else
LoadParams.Nthetaphsbins = 0;

%Sample rate at which data should be resampled
LoadParams.Samplerate = 60;

%Time smoothing window used for some variables when loading the data (e.g.
%smthBallSpeed)
LoadParams.LoadSmthTime = 150;

%Type of signal used to synch recordings from different sources
LoadParams.SynchType = 'Photodiode';%'SynchPulse';%'Speed';%'Reward';

%Orgin of the Synchsignal used as a reference for all types of data
LoadParams.SynchSignalRef = 'VR';
end