function Spk = LoadSpkData(Spkfilepath, sampleTimes, ZeroIdx, RateCorrection, Channels, clusterType, fieldnames)
%
load(Spkfilepath,'TempChans');
if ~iscell(TempChans)
    TempChans = {TempChans};
end

if nargin < 3
    ZeroIdx = 0;
end
if nargin < 4
    RateCorrection = 1;
end
if nargin < 5
    for iprobe = 1:numel(TempChans)
        Channels{iprobe} = 1:1000;
    end
end
if nargin < 6
    clusterType = '';
end
if nargin < 7
    fieldnames = {};
end
%would be nice to also have the cluster properties saved in that file

%Selecting the requested cells as indicated by Channels and clusterType and
%filling the fields of the Spk structure
icell = 0;
for iprobe = 1:numel(TempChans)
    if ~isnan(sum(Channels{iprobe}))
        for ichan = 1:length(TempChans{iprobe})
            if ismember(TempChans(ichan).ichan,Channels{iprobe})
                if isempty(clusterType) || ~isempty(strfind(TempChans(ichan).id,clusterType))
                    icell = icell + 1;
                    Spk.spikeIDs{icell} = TempChans(ichan).id;
                    Spk.chanIDs{icell} = TempChans(ichan).ichan;
                    Spk.ProbeIDs{icell} = iprobe;
                    Spk.spikeTimes{icell} = TempChans(ichan).spiketimes - ZeroIdx + 1;
                    Spk.spikeTimes{icell}(Spk.spikeTimes{icell}<0) = [];
                    for k = 1:numel(fieldnames)
                        Spk.(fieldnames{k}){icell} = TempChans(ichan).(fieldnames{k});
                    end
                end
            end
        end
    end
end

samplerate = mean(diff(sampleTimes));
Spk.spikeTrain = zeros(length(sampleTimes), numel(Spk.spikeTimes));
for icell = 1:numel(Spk.spikeTimes)
    Spk.spikeTrain(:,icell) = histcounts(Spk.spikeTimes{icell},sampleTimes+samplerate/2);
end
Spk.mua = sum(Spk.spikeTrain,2);
end