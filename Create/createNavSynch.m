function createNavSynch(inputPath, outputPath)
% function createNavSynch(inputPath, outputPath)
% creates a file with the synch signals of VR only
% specifically it saves: 'screenTimes', 'lick','reward', 'balldata' (speed
% info), 'APSyncPulse'

% Aman, Apr 2019

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

VRdata.TRIAL.time     = VRdata.TRIAL.time(1:totTrials,1:maxTrialLength);
VRdata.TRIAL.balldata = VRdata.TRIAL.balldata(1:totTrials,1:maxTrialLength,:);

% Making the bad points NaNs
VRdata.TRIAL.time(~goodTimes)       = NaN;
VRdata.TRIAL.time(:,1)              = NaN;
VRdata.TRIAL.time       = VRdata.TRIAL.time - min(VRdata.TRIAL.time(:)); % % AS change in 04/2018  % VRdata.TRIAL.time(1,2); %repmat(VRdata.TRIAL.time(:,2),1,size(VRdata.TRIAL.time(2,:),2));

VRdata.TRIAL.reward  = NaN*ones(size(VRdata.TRIAL.traj));
numRewards = length(VRdata.REWARD.count);

for irew = 1:numRewards
    VRdata.TRIAL.reward(VRdata.REWARD.TRIAL(irew),VRdata.REWARD.count(irew)) = VRdata.REWARD.TYPE(irew);
end

for n = 1:size(goodTimes,1)
    VRdata.TRIAL.balldata(n,min(find(~goodTimes(n,:))):end,:)  = NaN;
end

[screenTimes2 sortIDX] = sort(VRdata.TRIAL.time(:));
maxTime = min(find(isnan(screenTimes2)))-1;
screenTimes2 = screenTimes2(1:maxTime);

%% define required variables
screenTimes = screenTimes2;
sortIDX = sortIDX(1:maxTime);

lick      = VRdata.TRIAL.lick(sortIDX);
reward    = VRdata.TRIAL.reward(sortIDX);

balldata = VRdata.TRIAL.balldata(:,:,3);
balldata(balldata<0) = 0;
balldata = balldata(sortIDX);

APSyncPulse = VRdata.TRIAL.balldata(:,:,4);
APSyncPulse = APSyncPulse(sortIDX);

lickTimes = screenTimes(lick>0);
[~,idx] = (min(abs((repmat(esSync.sampleTimes,1,length(lickTimes)))...
- (repmat(lickTimes',size(esSync.sampleTimes,1),1)))));
esSync.lick      = zeros(size(esSync.sampleTimes));
esSync.lick(idx) = 1;

rewTimes = esSync.screenTimes(~isnan(reward));
[~,idx] = (min(abs((repmat(esSync.sampleTimes,1,length(rewTimes)))...
- (repmat(rewTimes',size(esSync.sampleTimes,1),1)))));
esSync.reward      = NaN*ones(size(esSync.sampleTimes));

try
    VRdata.REWARD.TYPE;
catch
    display('Rewards are wrong')
end
evenSampleTime = (1/60):(1/60):max(screenTimes);
esSync.sampleTimes = evenSampleTime';
esSync.screenTimes = screenTimes;

esSync.balldata         = interp1(esSync.screenTimes, balldata,   esSync.sampleTimes);
esSync.APSyncPulse      = interp1(esSync.screenTimes, APSyncPulse,esSync.sampleTimes);

%%
save(outputPath, 'screenTimes', 'lick', 'reward', 'balldata', 'APSyncPulse', 'esSync');