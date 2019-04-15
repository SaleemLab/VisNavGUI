function [tt,ww] = UnitGetNEV(elec,index)
% UnitGetNEV reads MUA data from NEV file
%
% tt = UnitGetNEV(elec,index) reads only the spike times
%
% [tt,ww] = UnitGetNEV(elec,index) reads also the waveforms
%
% index comes from ...
%
% 2008-02 MC from nevwaves



if(elec+1>length(index.idx))
    tt=[];
    ww=[];
    return;
end

idx = index.idx{elec+1};  %% find the indices

if(isempty(idx))
    tt=[];
    ww=[];
    return;
end

nspikes = length(idx);

tt = zeros(1,nspikes);
ww = zeros(48,nspikes);

if nargout == 1
    % read only the spike times tt
    for ispike=1:nspikes
        fseek(index.fid,index.offset+(idx(ispike)-1)*index.packetlen,'bof');
        tt(ispike) = fread(index.fid,1,'uint32=>double',1);  %% timestamp
    end
else
    for ispike=1:nspikes
        fseek(index.fid,index.offset+(idx(ispike)-1)*index.packetlen,'bof');
        tt(ispike) = fread(index.fid,1,'uint32=>double',1);  %% timestamp
        fseek(index.fid,3,'cof');
        ww(:,ispike) = fread(index.fid,48,'int16=>double');
    end
end

tt = tt/30;


