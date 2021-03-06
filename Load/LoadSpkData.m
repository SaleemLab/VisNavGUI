function Spk = LoadSpkData(Spkfilepath, sampleTimes, ZeroIdx, RateCorrection, Channels, clusterType, otherfields)
%Load the spike times in recording samples, convert them into predefined 
%recording units using ZeroIdx and RateCorrection 
%(t = RateCorrection*(spiketimes - ZeroIdx).
%inputs:
% - Spkfilepath: file path of the .mat file containing the spike sorting
%   data as a structure with fields as cell arrays corresponding at least
%   to spikeIDs, chanIDs, ProbeIDs, and spikeTimes.
% - sampleTimes: SampleTimes in synchronized recording units.
% - ZeroIdx: index offset to convert spike time units in synchronized 
%   recording units.
% - RateCorrection: relative sampling rate to convert spike time units in 
%   synchronized recording units.
% - Channels: chan id of clusters to load.
% - clusterType: cluster type to load (e.g. 'good')
% - otherfields: additional field names of the spkdata structure that needs
%   to be loaded into the output Spk structure
%
%Outputs:
% - Spk: contained spike data sampled in synchronized recording units, with
%   fields spikeIDs, chanIDs, ProbeIDs, spikeTimes, spikeTrain and any
%   additional fields requested in intput argument otherfields


S = load(Spkfilepath);
fields = fieldnames(S);
TempChans = S.(fields{1});

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
    otherfields = {};
end
%would be nice to also have the cluster properties saved in that file

%Selecting the requested cells as indicated by Channels and clusterType and
%filling the fields of the Spk structure
icell = 0;
for iprobe = 1:numel(TempChans)
    if ~isnan(sum(Channels{iprobe}))
        for ichan = 1:length(TempChans{iprobe}.spikeIDs)
            if ismember(TempChans{iprobe}.chanIDs{ichan},Channels{iprobe})
                if isempty(clusterType) || ~isempty(strfind(TempChans{iprobe}.spikeIDs{ichan},clusterType))
                    icell = icell + 1;
                    Spk.spikeIDs{icell} = TempChans{iprobe}.spikeIDs{ichan};
                    Spk.chanIDs{icell} = TempChans{iprobe}.chanIDs{ichan};
                    Spk.ProbeIDs{icell} = TempChans{iprobe}.ProbeIDs{ichan};
                    Spk.spikeTimes{icell} = (TempChans{iprobe}.spikeTimes{ichan} - ZeroIdx + 1)*RateCorrection;
                    Spk.spikeTimes{icell}(Spk.spikeTimes{icell}<0) = [];
                    for k = 1:numel(otherfields)
                        Spk.(otherfields{k}){icell} = TempChans{iprobe}.(otherfields{k}){ichan};
                    end
                end
            end
        end
    end
end

Spk.spikeTrain = zeros(length(sampleTimes), numel(Spk.spikeTimes));
sampleinterval = mean(diff(sampleTimes));
sampleTimes = [sampleTimes(1) - sampleinterval; sampleTimes(:)];
for icell = 1:numel(Spk.spikeTimes)
    Spk.spikeTrain(:,icell) = histcounts(Spk.spikeTimes{icell},sampleTimes+sampleinterval/2);
end
Spk.mua = sum(Spk.spikeTrain,2);
end