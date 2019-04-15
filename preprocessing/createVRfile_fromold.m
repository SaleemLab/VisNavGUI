function es = createVRfile_fromold(VRDIRname, animalname, series, iexp, samplerate, SmthTimeWindow)

VRfname = ls([VRDIRname filesep 'M*.mat']);

for f = 1:size(VRfname,1)
    VRdata = load([VRDIRname filesep VRfname(f,:)]);
    
    if nargin < 6
        SmthTimeWindow = 150;
    end
    
    %% To truncate the trials to only to max trial lengths
    VRdata.TRIAL.time(:,1) = [];
    goodTimes = ~(VRdata.TRIAL.time==0); % non-zero times
    goodTimes(:,1) = true;
    [~,colGT,~] = find(goodTimes); % find columns with non-zero times
    maxTrialLength = max(colGT);
    
    % JUL - 15.09.2015: changed the way the total number of trials is defined
    % by looking for the last trials that is longer than half a second instead
    % of :
    % totTrials = min([min(find(sum(double(goodTimes'))<=60)) size(VRdata.TRIAL.traj,1)]); % the first trial that is less than one second
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
    
    VRdata.TRIAL.eyeXpos = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    VRdata.TRIAL.eyeYpos = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    VRdata.TRIAL.pupilSize = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    VRdata.TRIAL.eyetime = zeros(size(VRdata.TRIAL.traj,1),size(VRdata.TRIAL.traj,2));
    
    VRdata.REWARD.count(VRdata.REWARD.TRIAL>totTrials) = [];
    VRdata.REWARD.TYPE(VRdata.REWARD.TRIAL>totTrials) = [];
    VRdata.REWARD.TRIAL(VRdata.REWARD.TRIAL>totTrials) = [];
    
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
    %     VRdata.TRIAL.time(:,1)              = NaN;
    VRdata.TRIAL.time       = VRdata.TRIAL.time - VRdata.TRIAL.time(1,1);%VRdata.TRIAL.time - VRdata.TRIAL.time(1,2);%repmat(VRdata.TRIAL.time(:,2),1,size(VRdata.TRIAL.time(2,:),2));
    VRdata.TRIAL.currTime   = VRdata.TRIAL.time - repmat(VRdata.TRIAL.time(:,1),1,size(VRdata.TRIAL.time,2));%VRdata.TRIAL.time - repmat(VRdata.TRIAL.time(:,2),1,size(VRdata.TRIAL.time,2));
    
    VRdata.TRIAL.eyeXpos(~goodTimes)    = NaN;
    VRdata.TRIAL.eyeYpos(~goodTimes)    = NaN;
    VRdata.TRIAL.pupilSize(~goodTimes)  = NaN;
    VRdata.TRIAL.eyetime(:,1)           = NaN;
    VRdata.TRIAL.eyetime    = VRdata.TRIAL.eyetime - VRdata.TRIAL.eyetime(1,1);%VRdata.TRIAL.eyetime - VRdata.TRIAL.eyetime(1,2);
    
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
    
    % try
    %     VRdata.trajspeed = fastsmooth(trajspeed,10,3,1);
    %     VRdata.ballspeed = fastsmooth(ballspeed,10,3,1);
    % catch
    if isfield(VRdata.EXP, 'wheelToVR')
        VRdata.trajspeed = trajspeed./VRdata.EXP.wheelToVR;
    else
        VRdata.trajspeed = trajspeed./2;
        % Dividing by 2 because I always use gain of 2. so that should be 1.
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
    
    if isfield(VRdata.EXP, 'CircularMaze')
        VRdata.CircularMaze = VRdata.EXP.CircularMaze;
    else
        VRdata.CircularMaze = false;
    end
    VRdata.SmthTimeWindow = SmthTimeWindow;
    
    es = makeVRevenSampled(VRdata, samplerate);
    es.animal = animalname;
    es.exp = iexp + f - 1;
    es.series = series;
    es.rewardTolerance = VRdata.EXP.rew_tol;
    
    save([VRDIRname filesep 'VRdata_' VRfname(f,:)], 'es', '-v7.3');
    
    %create Synchfile for VRdata
    Aperiodic = zeros(size(es.screenTimes));
    PhotoDiode = es.screenTimes;
    Licks = es.lick;
    RotEnc = es.ballspeed;
    save([VRDIRname filesep 'Synch_VR_' VRfname(f,:) '.mat'], 'Aperiodic', 'PhotoDiode', 'Licks', 'RotEnc', '-v7.3');
end
end

function es = makeVRevenSampled(VRdata, samplerate)


[es.screenTimes2, es.sortIDX] = sort(VRdata.TRIAL.time(:));
maxTime = min(find(isnan(es.screenTimes2)))-1;
if ~isempty(maxTime)
    es.screenTimes2 = es.screenTimes2(1:maxTime);
else
    maxTime = numel(es.screenTimes2);
end
if nargin < 4
    samplerate = 60;
end
% es.sampleRate = '60 Hz';
es.sampleRate = samplerate;%60;
es.screenTimes = es.screenTimes2;

evenSampleTime = (1/es.sampleRate):(1/es.sampleRate):max(es.screenTimes);
es.sampleTimes = evenSampleTime';

es.sortIDX = es.sortIDX(1:maxTime);

ballspeed = VRdata.ballspeed(es.sortIDX);
ballspeed(1) = 0;

trajspeed = VRdata.trajspeed(es.sortIDX);
trajspeed(1) = 0;

distTrav  = VRdata.distTrav(es.sortIDX);
distTrav(1)  = 0;

totDistTrav = VRdata.totDistTrav(es.sortIDX);
totDistTrav(1)  = 0;

trajPercent = VRdata.trajPercent(es.sortIDX);
trajPercent(1) = 0;

for i = 1:size(VRdata.TRIAL.trialIdx,1)
    VRdata.TRIAL.trialIdx(i,VRdata.TRIAL.trialIdx(i,:) == 0) = max(VRdata.TRIAL.trialIdx(i,:));
end

traj      = VRdata.traj(es.sortIDX);
currTime  = VRdata.TRIAL.currTime(es.sortIDX);
trialIdx  = VRdata.TRIAL.trialIdx(es.sortIDX);
lick      = VRdata.TRIAL.lick(es.sortIDX);
reward    = VRdata.TRIAL.reward(es.sortIDX);

eyeXpos   = VRdata.TRIAL.eyeXpos(es.sortIDX);
eyeYpos   = VRdata.TRIAL.eyeYpos(es.sortIDX);
pupilSize   = VRdata.TRIAL.pupilSize(es.sortIDX);

lickTimes = es.screenTimes(lick>0);
idx = zeros(1,numel(lickTimes));
Tsample = min(diff(es.sampleTimes));
for l = 1:numel(lickTimes)
    idx(l) = find((es.sampleTimes - lickTimes(l)) > -Tsample,1,'first');
end
es.lick      = zeros(size(es.sampleTimes));
es.lick(idx) = 1;

rewTimes = es.screenTimes(~isnan(reward));
[~,idx] = (min(abs((repmat(es.sampleTimes,1,length(rewTimes)))...
- (repmat(rewTimes',size(es.sampleTimes,1),1)))));
es.reward      = NaN*ones(size(es.sampleTimes));
try
    es.reward(idx) = VRdata.REWARD.TYPE;
catch
    disp('Rewards are wrong')
    try
        es.reward(idx) = VRdata.REWARD.TYPE(1:numel(rewTimes));
    catch
        rewTimes = es.screenTimes(~isnan(reward) & reward == 2);
        es.reward(idx) = VRdata.REWARD.TYPE(1:numel(rewTimes));
    end
end
roomlength = VRdata.EXP.l;

es.trialID       = round(interp1(es.screenTimes, trialIdx, es.sampleTimes, 'nearest'));
%JUL: changed interpolation method on traj and trajpercent to avoid
%articfact due to missing frame at the beginning of every trial
% es.traj          = interp1(es.screenTimes, traj,es.sampleTimes, 'nearest');
es.traj          = mod((interp1(es.screenTimes, unwrap(traj/roomlength*2*pi)*roomlength/(2*pi),es.sampleTimes, 'linear')),roomlength);%, 'previous')),roomlength);%
% es.traj          = mod((interp1(es.screenTimes, unwrap(traj/round(max(traj))*2*pi)*round(max(traj))/(2*pi),es.sampleTimes, 'linear')),round(max(traj)));
es.trajspeed     = interp1(es.screenTimes(~isnan(trajspeed)), trajspeed(~isnan(trajspeed)), es.sampleTimes, 'nearest');
es.ballspeed     = interp1(es.screenTimes(~isnan(ballspeed)), ballspeed(~isnan(ballspeed)), es.sampleTimes, 'nearest');
es.currTime      = interp1(es.screenTimes, currTime,  es.sampleTimes, 'linear');

es.distTrav      = interp1(es.screenTimes, distTrav,  es.sampleTimes, 'nearest');
es.totDistTrav      = interp1(es.screenTimes, totDistTrav,  es.sampleTimes, 'nearest');
% es.trajPercent      = interp1(es.screenTimes, trajPercent,  es.sampleTimes, 'nearest');
es.trajPercent          = mod((interp1(es.screenTimes, unwrap(trajPercent/roomlength*2*pi)*roomlength/(2*pi),es.sampleTimes, 'linear')),roomlength);%, 'previous')),roomlength);%

es.smthBallSpd   = NaN*ones(size(es.traj));
es.smthBallSpd(~isnan(es.ballspeed))   = smthInTime(es.ballspeed(~isnan(es.ballspeed)), es.sampleRate, VRdata.SmthTimeWindow, [], [], 'median');

es.smthBallAcc   = NaN*ones(size(es.traj));
es.smthBallAcc(~isnan(es.smthBallSpd))   = [diff(es.smthBallSpd(~isnan(es.smthBallSpd))) ; 0] ;
es.smthBallAcc(~isnan(es.smthBallAcc))   = smthInTime(es.smthBallAcc(~isnan(es.smthBallAcc)), es.sampleRate, VRdata.SmthTimeWindow);

es.smthTrajSpd   = NaN*ones(size(es.traj));
es.smthTrajSpd(~isnan(es.trajspeed))   = smthInTime(es.trajspeed(~isnan(es.trajspeed)), es.sampleRate, VRdata.SmthTimeWindow, [], [], 'median');

es.eyeXpos = interp1(es.screenTimes, eyeXpos, es.sampleTimes, 'linear');
es.eyeYpos = interp1(es.screenTimes, eyeYpos, es.sampleTimes, 'linear');
es.pupilSize = interp1(es.screenTimes, pupilSize, es.sampleTimes, 'linear');

es.screenTimes     = interp1(es.screenTimes, es.screenTimes, es.sampleTimes, 'previous');

es.contrast = NaN*ones(size(es.traj));
es.start    = NaN*ones(size(es.traj));
es.gain     = NaN*ones(size(es.traj));
es.blanks   = NaN*ones(size(es.traj));
es.active   = NaN*ones(size(es.traj));
es.rewardPos= NaN*ones(size(es.traj));
es.outcome  = NaN*ones(size(es.traj));
es.roomLength= NaN*ones(size(es.traj));

trialIDs = unique(es.trialID);
numTrials = length(trialIDs);

for itrial = 1:numTrials
    es.contrast(es.trialID==itrial)    = VRdata.TRIAL.trialContr(trialIDs(itrial));
    es.start(es.trialID==itrial)       = VRdata.TRIAL.trialStart(trialIDs(itrial));
    es.gain(es.trialID==itrial)        = VRdata.TRIAL.trialGain(trialIDs(itrial));
    es.blanks(es.trialID==itrial)      = VRdata.TRIAL.trialBlanks(trialIDs(itrial));
    es.active(es.trialID==itrial)      = VRdata.TRIAL.trialActive(trialIDs(itrial));
    es.rewardPos(es.trialID==itrial)   = VRdata.TRIAL.trialRewPos(trialIDs(itrial));
    es.outcome(es.trialID==itrial)     = VRdata.TRIAL.trialOutcome(trialIDs(itrial));
    es.roomLength(es.trialID==itrial)  = VRdata.TRIAL.trialRL(trialIDs(itrial));
end

if sum(isnan(es.gain)) ~= 0
    disp('WARNING: NaNs in gain vector');
    gainval = unique(es.gain(~isnan(es.gain)));
    es.gain(isnan(es.gain)) = 0.1*round((es.trajspeed(isnan(es.gain))./es.ballspeed(isnan(es.gain)))/0.1);
    es.gain(~ismember(es.gain,gainval)) = NaN;
end

if nargin==3
    disp('WARNING!!!! temporary fix. Check the photodiode signalling');
    evenSampleTime = evenSampleTime + screenTimes(1);
end

if VRdata.CircularMaze
    wrapfactor = 2;%1;% 
    disp(['wrapping up circular maze by a factor of ' num2str(wrapfactor)]); 
    roomlength = VRdata.EXP.l;
    es = VRwrap(es,roomlength,wrapfactor,true);
    es = VRredefineOutcome(es);
else
    es = redefineOutcome_AS(es);
end

es.CircularMaze = VRdata.CircularMaze;
end