function waveform = nswaves(elec,varargin)
% NSWAVES   Gets wavesforms from .ns* file
%   W = NSWAVES(E) returns the vector W containing the raw trace
%   for electrode # E.
% 
% AZ 2010-03-16 Created

%% Initialization
% Load fid
global pepNEV
% Check if electrode # is valid
elecs = find(pepNEV.ids == elec);
if isempty(elecs)
   error('Electrode %u is out of range.\n',elec);
end
% Skip header
offset = 32 + 4*pepNEV.nchan;
% Load file into memory
% m = memmapfile(pepNEV.index.fid)


%% Read the data
all_data = fread(pepNEV.index.fid,inf,'int16=>double');
waveform = reshape(all_data,numel(elecs),numel(all_data)/numel(elecs));