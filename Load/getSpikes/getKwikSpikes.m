function chans = getKwikSpikes(DIRname, animal, iseries, iexp, igroup, addInfo)

basename = [animal '_s' num2str(iseries) '_1'];
basenamekwik = [animal '_s' num2str(iseries) '_1'];
if nargin>4
    if ~isempty(addInfo)
        basename = [animal '_s' num2str(iseries) '_' addInfo];
        basenamekwik = [animal '_s' num2str(iseries) '_' addInfo];
    end
end

if isfolder([DIRname filesep addInfo])
    load([DIRname filesep addInfo filesep basename]);
    if isdir([DIRname filesep addInfo filesep 'tet' num2str(igroup)])
        kwikFile = [DIRname filesep addInfo filesep 'tet' num2str(igroup) filesep basenamekwik '_tet' num2str(igroup) '.kwik'];
        clusterinfoFile = [DIRname filesep addInfo filesep 'tet' num2str(igroup) filesep basenamekwik '_tet' num2str(igroup) '_clusterinfo.mat'];
    else
        kwikFile = [DIRname filesep addInfo filesep basenamekwik '.kwik'];
        clusterinfoFile = [DIRname filesep addInfo filesep basenamekwik '_clusterinfo.mat'];
    end
    igroupfile = 0;
else
    load([DIRname filesep basename]);
    kwikFile = [DIRname filesep basenamekwik '.kwik'];
    igroupfile = igroup;
end

expIdx = find(SELECTED_EXPERIMENTS==iexp);
expEnds = cumsum(lims);
if expIdx>1
    startTime = expEnds(expIdx-1);
else
    startTime = 1;
end
endTime = expEnds(expIdx);

spkTimes = hdf5read(kwikFile, ['/channel_groups/' num2str(igroupfile) '/spikes/time_samples']);
spkClus = hdf5read(kwikFile, ['/channel_groups/' num2str(igroupfile) '/spikes/clusters/main']);
cellIDs = unique(spkClus);
ncells  = length(cellIDs);


spkClus(spkTimes>endTime) = [];
spkTimes(spkTimes>endTime) = [];

spkTimes = double(spkTimes);
spkTimes = spkTimes - startTime;
spkClus(spkTimes<0) = [];
spkTimes(spkTimes<0) = [];


temp =  h5info(kwikFile,  '/recordings/0/');
try
    sampleRate = double((temp.Attributes(7).Value));
catch
    sampleRate = double((temp.Attributes(3).Value));
end
spkTimes = double(spkTimes);%./sampleRate;

noise_list = [];

if exist(clusterinfoFile)
    load(clusterinfoFile);
else
    clus = [];
end
if numel(clus) ~= ncells
    warning('clusterinfo do not match with sorted kwik file: run createClusterinfo.m again');
    clus = [];
end

for icell = 1:ncells
    chans(icell).spiketimes = spkTimes(spkClus==cellIDs(icell));
    chans(icell).ichan = igroup;
    if ~isempty(clus)
        chans(icell).bestchan = clus(icell).bestchan;
        chans(icell).waveform = clus(icell).waveform(clus(icell).bestchan,:);
    else
        chans(icell).bestchan = NaN;
        chans(icell).waveform = NaN;
    end
    chans(icell).iexp = iexp;
    chans(icell).icell = cellIDs(icell);
    chans(icell).sampleRate = sampleRate;
    
    temp = h5info(kwikFile, ['/channel_groups/' num2str(igroupfile) '/clusters/main/' num2str(cellIDs(icell)) '/']);
    
    for idx = 1:length(temp.Attributes)
        if strcmp(temp.Attributes(idx).Name, 'cluster_group')
            ilabel = temp.Attributes(idx).Value;
            break
        end
        idx = idx + 1;
    end
    if idx <= length(temp.Attributes)
        ilabel = temp.Attributes(idx).Value;
        temp = h5info(kwikFile, ['/channel_groups/' num2str(igroupfile) '/cluster_groups/main']);
        labelType = eval(['temp.Groups(' num2str(ilabel+1) ').Attributes(end).Value']);
        if strcmp(labelType,'Noise')
            noise_list = [noise_list icell];
        end
        chans(icell).id = [basename '_c' num2str(igroup) '_' labelType '_cluster' num2str(cellIDs(icell))];
    else
        chans(icell).id = [basename '_c' num2str(igroup) '_none_cluster' num2str(cellIDs(icell))];
    end
end
chans(noise_list) = [];

