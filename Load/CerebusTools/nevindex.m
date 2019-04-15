function r = nevindex(filename)

global pepNEV

fid = fopen(filename, 'r', 'ieee-le');
if fid == -1, 
    r = -1;
    return;
end;

% fseek(fid, 12, 'bof');
fseek(fid, 12, 'bof');
index.offset = fread(fid, 1, 'uint32'); % # bytes in headers
index.packetlen = fread(fid, 1, 'uint32'); % # bytes per data packet
    
TimeResOfTimeStamps = fread(fid,1,  'uint32=>double'); % AZ 2010-01-12
SamplingRateInKHZ   = fread(fid,1,  'uint32=>double')/1000;

% eof=0;
fseek(fid,index.offset+4,'bof');  %% because we are not reading timestamps
[packetid,index.N] = fread(fid,inf,'uint16=>uint8',index.packetlen-2);

index.fid = fid;

nelec = double(max(packetid));
for j = 0:nelec
    index.idx{j+1} = find(packetid==j); 
end

pepNEV.index = index;

[sync.timestamps,sync.digital] = nevdigin;

%% 2010-09-02 ND the below two corrections for the on/off times from the
%% photodiode work better if starting with alternating -1's and -2's.
idx = find(diff(sync.digital)==0)+1;
sync.timestamps(idx)=[];
sync.digital(idx) = [];

%% there is a funny bug that sometimes we get a -2 followed by -1 with less
%% than 1msec distance.   We remove them here...
% 
idx1 = find(sync.digital == -2);
% idx2 = find(diff(sync.timestamps)<=SamplingRateInKHZ); % 2010-07-19: was <=5
idx2 = find(diff(sync.timestamps)<=2*SamplingRateInKHZ); % 2010-09-02: allowed to be a 2 ms rather than 1 ms
idx = intersect(idx1,idx2);
sync.timestamps(idx)=[];
sync.digital(idx) = [];

%% 2010-07-19 AZ & MS sometimes we see two -1's 3 timepoints apart.  We
%% remove the later ones (and -2's 3 apart as well) here.
idx = find(diff(sync.digital)==0)+1;
sync.timestamps(idx)=[];
sync.digital(idx) = [];

%% Now check for sync from imaging computer

[isync,w] = nevwaves(131); %#ok<NASGU>
clear w;

%AZ20090605: less obscure display message
%ND20101020: but still not very clear - what is this supposed to be?
disp(['nevindex: isync size is ',num2str(size(isync))]);

pepNEV.sync = sync;
pepNEV.isync = isync;  %% it was empty before!

idx = findstr(filename,filesep);
idxfn = [filename(1:idx(end)-1) filesep 'index_' ...
   strtok(filename(idx(end)+1:end),'.') '.mat'];
save(idxfn,'pepNEV');

fclose(fid);
r=0;
return;
