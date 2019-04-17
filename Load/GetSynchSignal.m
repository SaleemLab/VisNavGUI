function [SynchSignal, SynchUpSampling] = GetSynchSignal(filepath, SynchType)
m = matfile(filepath);
SynchSignal = m.(SynchType);
SynchUpSampling = m.UpSampling;
end
