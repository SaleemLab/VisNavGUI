function DIRS = SetDirectories
%set the path of directories where to find the data
%should call a common setDirs function
serverName = 'C:\';
DIRS.data           = fullfile(serverName,'Data','trodes');%directory with visual stimulation protocols
DIRS.spikes         = fullfile(serverName,'Data','Spikes');
DIRS.EyeCamera      = fullfile(serverName,'Data','EyeCamera');%directory with eye camera recordings
DIRS.EyeTrack       = fullfile(serverName,'Data','EyeTrack');%directory with eye tracking data
DIRS.xfiles         = fullfile(serverName,'Data','xfiles');%directory with visual stimulation protocols
DIRS.Cerebus        = fullfile(serverName,'Data','Cerebus');%directory of raw data recording for preprocessing LFP
DIRS.stimInfo       = fullfile(serverName,'Data','stimInfo');
DIRS.behavior       = fullfile(serverName,'Data','behavior');
DIRS.multichanspikes= fullfile(serverName,'Data','multichanspikes');%directory with sorted spikes
DIRS.ball           = fullfile(serverName,'Data','ball');%directory with behavior data
DIRS.expInfo        = fullfile(serverName,'Data','expInfo');

%then some additional fields on where to save the data can be added
%for instance:
DIRS.myfolder = 'D:\Data';
end