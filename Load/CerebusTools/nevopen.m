function [r,SamplingRateInKHZ,nChans,nSamplesPerPacket] = nevopen(filename)

%% dlr

global pepNEV;

if(~isempty(pepNEV))
    try
        nevclose;
    end
    pepNEV = [];
end

idx = findstr(filename,filesep);
idxfn = [filename(1:idx(end)-1) filesep 'index_' strtok(filename(idx(end)+1:end),'.') '.mat'];
if exist(idxfn,'file')
    disp('NEV index found')    
    load(idxfn);
    global pepNEV
    fid = fopen(filename, 'r', 'ieee-le');
    if fid == -1,
        r = -1;
        return;
    end;
    pepNEV.index.fid = fid;
    % AZ 2010-02-05
%     nBytesInHeaders     = pepNEV.index.offset;
%     nBytesPerPacket     = pepNEV.index.packetlen;
    nSamplesPerPacket   = (pepNEV.index.packetlen - 8) / 2;
    nChans = numel(pepNEV.index.idx)-1;
    % AZ 2009-03-04
    fseek(fid,20,'bof'); % (standard header is 1st 336 bytes)
    TimeResOfTimeStamps = fread(fid,1,  'uint32=>double'); % AZ 2010-01-12
    SamplingRateInKHZ   = fread(fid,1,  'uint32=>double')/1000;
    % AZ 2010-01-12
    TimeIndexZeroGMT    = fread(fid,8,  'uint16=>double');
    CreatorApplication  = fread(fid,32, 'int8=>char'    );
    CreatorApplication  = CreatorApplication';
    CommentField        = fread(fid,256,'int8=>char'    );
    pepNEV.TimeIndex    = TimeIndexZeroGMT;
    if unique(double(CommentField))==0
       CommentField     = [];
    else
       CommentField     = CommentField';
    end
    
    r = 1;
    return;
else
    disp('Generating NEV index file');
end

r = nevindex(filename);  %% first index it...
if r == -1
    error('no such file: %s', filename);
end
[r,SamplingRateInKHZ,nChans,nSamplesPerPacket] = nevopen(filename);

return;
