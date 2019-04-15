function [t,w] = nevwaves(elec,varargin)
% NEVWAVES   Gets wavesforms from .nev file
%   [T,W] = NEVWAVES(E) returns timestamps T and waveforms W for each
%   spike. T is a vector containing a single timestamp per spike, W is a
%   matrix with columns containing the waveforms (DEFAULT length: 48) of 
%   each spike.
%
%   [T,W] = NEVWAVES(E, [FIRST LAST]) returns timepoints T and waveforms W
%   for only the FIRST to LAST spikes. If electrode has less than LAST
%   spikes, all spikes including the final spike will be returned.
%
%   [T,W] = NEVWAVES(E, NRAND) returns timepoints T and waveforms W
%   for a random sample of size NRAND spikes.
% 
%   [T,W] = NEVWAVES(E, [], NSAMPLESPERPACKET) returns timepoints T and waveforms W
%   for each spike, with a user-specified NSAMPLESPERPACKET.
%
% change log
% Dario Ringach created it
% 2008-08 LB corrected some matlab complaints and added some comments
% 2009-01 AZ allowed user to choose random sample of L # of spikes (didn't
%            work)
% 2010-01 AZ properly implemented random sampling
% 2010-05 AZ now allows waveforms of specified length to be loaded
% 2011-05 MC sped up in case we don't want the waveforms

%% Get waveforms for electrode
global pepNEV;

DoWaves = (nargout>1);

if elec+1 > length(pepNEV.index.idx)
    t=[];
    w=[];
    return;
end

idx = pepNEV.index.idx{elec+1};  %% find the indices

if nargin > 1
    L = varargin{1};
    
    if numel(L) == 2        % L = [FIRST LAST]
       idx = idx(max(1,L(1)):min(L(2),numel(idx)));
    elseif numel(L) == 1    % L = random sample, size L
       randIX = randperm(numel(idx));
       idx = sort(idx(randIX(1:min(L,numel(idx)))));
    end
    if nargin > 2
       nSamplesPerPacket = varargin{2};
    end
end
% Set default nSamplesPerPacket if empty/zero
if ~exist('nSamplesPerPacket','var') || isempty(nSamplesPerPacket) || ...
      nSamplesPerPacket < 1
   nSamplesPerPacket = 48;
end

if isempty(idx)
    t=[];
    w=[];
    return;
end

N = length(idx); % number of spikes

t = zeros(1,N);  % timestamps, one per spike
if DoWaves
    w = zeros(nSamplesPerPacket,N);
end

for j=1:N
    fseek(pepNEV.index.fid,pepNEV.index.offset+(idx(j)-1)*pepNEV.index.packetlen,'bof');
    t(j) = fread(pepNEV.index.fid,1,'uint32=>double',1);  %% timestamp
    if DoWaves
        fseek(pepNEV.index.fid,3,'cof');
        w(:,j) = fread(pepNEV.index.fid,nSamplesPerPacket,'int16=>double');  %% waveform
    end
end
