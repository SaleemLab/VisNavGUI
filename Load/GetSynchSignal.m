function SynchSignal = GetSynchSignal(filepath, SynchType)
m = matfile(filepath,SynchType);
SynchSignal = m.(SynchType);
end