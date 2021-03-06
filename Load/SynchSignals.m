function [Signal_ZeroIdx, Signal_relativeRate] = SynchSignals(SynchRef, SynchSignal, SynchUpSamplingRef, SynchUpSampling)
%Returns the time of start of SynchRef in SynchSignal and the sampling rate
%of SynchSignal relative to SynchRef. These are expressed in units of the
%original sampling associated to SynchSignal, i.e. taking into account the
%upsampling of the Synch signals.

SynchRef = SynchRef(:) - mean(SynchRef);
SynchSignal = SynchSignal(:) - mean(SynchSignal);

%Finding pulses in the SynchSignals
%if the 1st value in the ref is lower than the mean, we look at up
%transitions, otherwise, we look at down transitions
if SynchRef(1) < 0
    SynchTimesRef = find(diff(sign(SynchRef)) > 0);
    SynchTimes = find(diff(sign(SynchSignal)) > 0);
else
    SynchTimesRef = find(diff(sign(SynchRef)) < 0);
    SynchTimes = find(diff(sign(SynchSignal)) < 0);
end

%Measuring the intervals between successive up transiton
SynchIntRef = diff(SynchTimesRef);
SynchInt = diff(SynchTimes);

%Selecting a pattern of n intervals in the middle of the reference signal,
%where n is a fraction of the total number of intervals
n = 1/20;
NintRef = numel(SynchIntRef);
Npattern = floor((NintRef)*n);
idxRef = floor((NintRef/2)-round(Npattern/2)):(floor(NintRef/2)+round(Npattern/2));
Npattern = numel(idxRef);
patternRef = SynchIntRef(idxRef);
patternRef = (patternRef - mean(patternRef));
patternRef = patternRef./sqrt(sum(patternRef.^2));

%creating a matrix with shifted versions of the sequence of intervals in
%the SynchSignal
SynchInt = repmat(SynchInt(:),[1 Npattern]);
for k = 0:(Npattern-1)
    SynchInt(:,k + 1) = circshift(SynchInt(:,k + 1), -k);
end
SynchInt = (SynchInt - repmat(mean(SynchInt,2), [1 Npattern]));
SynchInt = SynchInt./repmat(sqrt(sum(SynchInt.^2,2)), [1 Npattern]);

%Estimating the time of start of SynchRef in SynchSignal in units of number
%of pulses, by computing the normalized correlation between the pattern of 
%intervals selcted in SynchIntRef with the intervals in SynchInt.
SynchCorr = SynchInt*patternRef(:);
[~, correction] = max(SynchCorr);
SynchIntcorrection = correction - idxRef(1) + 1;

%Will now estimate the offset and sampling rate of SynchSignal 
%relative to SynchRef
%First, we select the pulse times (in recording units) that are common 
%between the two signals
SynchTimesRefStartIdx = max(1,-(SynchIntcorrection - 1) + 1);
SynchTimesStartIdx = max(1,SynchIntcorrection);
SynchTimesEndIdx = min(numel(SynchTimes),(SynchTimesStartIdx+numel(SynchTimesRef(SynchTimesRefStartIdx:end))-1));
t = SynchTimes(SynchTimesStartIdx:SynchTimesEndIdx);%-SynchTimes(SynchTimesStartIdx);
tref = SynchTimesRef(SynchTimesRefStartIdx:(SynchTimesRefStartIdx + numel(t) - 1)); %-SynchTimesRef(SynchTimesRefStartIdx);

%Second, we estimate the index offset and relative rate by fitting the 
%following equation: tref = p(1)*t + p(2)
foptions = fitoptions('Method','LinearLeastSquares');
g = fit(t,tref,'poly1',foptions);
p = coeffvalues(g);
if abs(numel(SynchSignal) - p(1)*numel(SynchSignal)) < 1
    p(1) = 1;
end

%Fimally, we convert p(1) and p(2) in units of the original recording 
%associated to SynchSignal, i.e. taking into account the upsampling of the Synch signals
Signal_relativeRate = p(1) * SynchUpSampling/SynchUpSamplingRef;
Signal_ZeroIdx = -floor(p(2) / SynchUpSampling);% + 1;
end