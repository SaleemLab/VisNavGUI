function es = LoadNavData(inputPath, sampleTimes, zeroCorrection, smthTimeWindow)

% Inherited from VRWheelLoad_SL
VRdata = load(inputPath);

%% To truncate the trials to only to max trial lengths
goodTimes = ~(VRdata.TRIAL.time==0); % non-zero times
goodTimes(:,1) = true;
[~,colGT,~] = find(goodTimes); % find columns with non-zero times
maxTrialLength = min([max(colGT) size(VRdata.TRIAL.posdata,2)]);

totTrials = min([min(find(sum(double(goodTimes'))>30,1,'last')) size(VRdata.TRIAL.traj,1)]); % the last trial that is longer than half second
if isfield(VRdata.TRIAL,'trialContr')
    totTrials = min([totTrials sum(~isnan(VRdata.TRIAL.trialContr))]);
end
if isempty(totTrials)
    totTrials = size(goodTimes,1);
end
goodTimes = goodTimes(1:totTrials,1:maxTrialLength);

% Truncating all relevant variables
VRdata.TRIAL.posdata  = VRdata.TRIAL.posdata(1:totTrials,1:maxTrialLength,:);
VRdata.TRIAL.traj     = VRdata.TRIAL.traj(1:totTrials,1:maxTrialLength);
VRdata.TRIAL.time     = VRdata.TRIAL.time(1:totTrials,1:maxTrialLength);
VRdata.TRIAL.balldata = VRdata.TRIAL.balldata(1:totTrials,1:maxTrialLength,:);
VRdata.TRIAL.trialIdx = VRdata.TRIAL.trialIdx(1:totTrials,1:maxTrialLength);
VRdata.TRIAL.trialBlanks = VRdata.TRIAL.trialBlanks(1:totTrials);
if isfield(VRdata.TRIAL,'lick')
    VRdata.TRIAL.lick  = VRdata.TRIAL.lick(1:totTrials,1:maxTrialLength);
end
if isfield(VRdata.TRIAL,'trialRL')
    VRdata.TRIAL.trialRL  = VRdata.TRIAL.trialRL(1:totTrials);
end

% Making the bad points NaNs
VRdata.TRIAL.time(~goodTimes)       = NaN;
VRdata.TRIAL.traj(~goodTimes)       = NaN;
VRdata.TRIAL.trialIdx(~goodTimes)   = NaN;
VRdata.TRIAL.time(:,1)              = NaN;
VRdata.TRIAL.time       = VRdata.TRIAL.time - min(VRdata.TRIAL.time(:)); %% ???????Should we be doing the zero correction here???????
VRdata.TRIAL.currTime   = VRdata.TRIAL.time - repmat(VRdata.TRIAL.time(:,2),1,size(VRdata.TRIAL.time,2));

VRdata.TRIAL.reward  = NaN*ones(size(VRdata.TRIAL.traj));
numRewards = length(VRdata.REWARD.count);

for irew = 1:numRewards
    VRdata.TRIAL.reward(VRdata.REWARD.TRIAL(irew),VRdata.REWARD.count(irew)) = VRdata.REWARD.TYPE(irew);
end

for n = 1:size(goodTimes,1)
    VRdata.TRIAL.posdata(n,min(find(~goodTimes(n,:))):end,:)  = NaN;
    VRdata.TRIAL.balldata(n,min(find(~goodTimes(n,:))):end,:)  = NaN;
end

%%
VRdata.traj = VRdata.TRIAL.traj;

trajspeed = zeros(size(VRdata.traj));
trajspeed(:,2:end) = (diff(VRdata.traj')');
trajspeed(trajspeed > (max(VRdata.TRIAL.traj(:))/2)) = NaN;
trajspeed(trajspeed < (-max(VRdata.TRIAL.traj(:))/2)) = NaN;
trajspeed(:,2:end) = trajspeed(:,2:end)./(diff(VRdata.TRIAL.time')');
% trajspeed = trajspeed

trajspeed(trajspeed<-50) = NaN;
trajspeed(trajspeed>200) = NaN;
% idx = find(trajspeed);
% trajspeed(idx(1)) = trajspeed(idx(1)) - 1000;

balldata = VRdata.TRIAL.balldata(:,:,3);

ballspeed = zeros(size(balldata));
ballspeed(:,2:end) = [balldata(:,2:end)./(diff(VRdata.TRIAL.time')')];

balldata(balldata<0) = 0;
distTrav = zeros(size(balldata));
totDistTrav = zeros(size(balldata));
trajPercent = zeros(size(VRdata.traj));
for itrial = 1:size(distTrav,1)
    distTrav(itrial,VRdata.traj(itrial,:)>0 & ~isnan(VRdata.traj(itrial,:))) = balldata(itrial,VRdata.traj(itrial,:)>0 & ~isnan(VRdata.traj(itrial,:)));
    totDistTrav(itrial,VRdata.traj(itrial,:)>=0 & ~isnan(VRdata.traj(itrial,:))) = balldata(itrial,VRdata.traj(itrial,:)>=0 & ~isnan(VRdata.traj(itrial,:)));
    
    trajPercent(itrial,:) = VRdata.traj(itrial,:)./(VRdata.TRIAL.trialRL(itrial));
end

distTrav = cumsum(distTrav,2);
distTrav(isnan(VRdata.traj)) = nan;

totDistTrav = cumsum(totDistTrav,2);
totDistTrav(isnan(VRdata.traj)) = nan;

if isfield(VRdata.EXP, 'wheelToVR')
    VRdata.trajspeed = trajspeed./VRdata.EXP.wheelToVR;
end

VRdata.ballspeed = ballspeed;
VRdata.distTrav  = distTrav;
VRdata.totDistTrav  = totDistTrav;
VRdata.trajPercent  = trajPercent;

numTrials = size(VRdata.traj,1);
% Truncating all the incomplete trial info
VRdata.TRIAL.trialContr     = VRdata.TRIAL.trialContr(1:numTrials);
VRdata.TRIAL.trialStart     = VRdata.TRIAL.trialStart(1:numTrials);
VRdata.TRIAL.trialGain      = VRdata.TRIAL.trialGain(1:numTrials);
VRdata.TRIAL.trialActive    = VRdata.TRIAL.trialActive(1:numTrials);
VRdata.TRIAL.trialRewPos    = VRdata.TRIAL.trialRewPos(1:numTrials);
VRdata.TRIAL.trialOutcome   = VRdata.TRIAL.trialOutcome(1:numTrials);
VRdata.TRIAL.trialRL        = VRdata.TRIAL.trialRL(1:numTrials);

VRdata.TRIAL.trialOutcome(VRdata.TRIAL.trialOutcome==0 & ...
    max(VRdata.TRIAL.currTime')>=VRdata.EXP.maxTrialDuration) = 5;

if iexp > 200
    [VRdata_o] = load([VRDIRname filesep VRdata.EXP.replayed_filename]);
else
    VRdata_o = [];
end

es = VREvenSample(VRdata, sampleTimes, sort(VRdata.TRIAL.time(:))-zeroCorrection); %% ???????Should we be doing the zero correction here???????
es.iexp = zeros(size(es.sampleTimes));
es.iexp(:) = iexp;