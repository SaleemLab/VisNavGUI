function [t,w] = nevdigin

%% Get times and values for digin

global pepNEV;

idx = pepNEV.index.idx{1};
t = zeros(1,length(idx));
w = zeros(1,length(idx));

for j = 1:length(idx)
    fseek(pepNEV.index.fid,pepNEV.index.offset+(idx(j)-1)*pepNEV.index.packetlen,'bof');
    t(j) = fread(pepNEV.index.fid,1,'uint32=>double',1);  %% timestamp
    fseek(pepNEV.index.fid,3,'cof');
    w(:,j) = fread(pepNEV.index.fid,1,'int16=>double');
end
