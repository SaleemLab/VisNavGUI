function t = nevtimes(electrodes,varargin)
% NEVTIMES   Gets spike times from .nev file
%   T = NEVTIMES(E) returns timestamps T for each spike on each of the electrodes in E.
%   T is a cellarray of vectors.
%
%   The function is based upon nevwaves, but since it allows to get the
%   times of the spikes on all the electodes at once, is much faster.
%   To do that, it reads big chanks [of ~1GB] of the nev file 
%   (although not the whole file) into RAM at once.
%
%   These probably work as in nevwaves:
%   T = NEVTIMES(E, [FIRST LAST])
%   T = NEVTIMES(E, NRAND)
%   T = NEVTIMES(E, [], NSAMPLESPERPACKET)
%
%   2011-09 MO

global pepNEV;
if max(electrodes)+1 > length(pepNEV.index.idx)
  t=[];
  return;
end;

t{max(electrodes)} = [];

CHUNK_SIZE = 2^28; % of uint32, i.e. ~ 1GB

fseek(pepNEV.index.fid,pepNEV.index.offset,'bof');
intsReadSoFar = 0;
count = CHUNK_SIZE;
while count == CHUNK_SIZE
  [nevContent, count] = fread(pepNEV.index.fid,CHUNK_SIZE,'uint32=>uint32');
  intsReadSoFar = intsReadSoFar + count;
  for elec = electrodes
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
    
    N = length(idx); % number of spikes
    
    if ~isempty(idx) && isempty(t{elec})
      t{elec} = zeros(1,N);  % preallocate the array, prior to using it...
    end;
    
    for j = 1:N
      pos = 1+(idx(j)-1)*pepNEV.index.packetlen/4; % each uint32 is 4 bytes...
      if pos <= intsReadSoFar-count
        continue; % pos refers to a portion that was already read
      elseif pos > intsReadSoFar
        continue; % we are trying to access a portion that was not read so far
      else
        t{elec}(j) = double(nevContent(pos-(intsReadSoFar-count)));
      end;
    end;
    
  end; % for elec = electrodes
end; % while count == CHUNK_SIZE
