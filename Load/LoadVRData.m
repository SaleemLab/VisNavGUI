function LoadVRData(EXP, P, animalname, iseries, iexplist, processedfiles)
P.LoadParams.animal = animalname;
P.LoadParams.series = iseries;
P.LoadParams.exp = iexplist;

obj.SpeedThresh = speed_th;
obj.SmthTimeWindow = SmthTimeWindow;
if isdir([P.DIRS.data2p filesep P.LoadParams.animal filesep num2str(P.LoadParams.series)])
    EXP.Nav = VRLoadMultipleExpts(P.LoadParams.animal, P.LoadParams.series, P.LoadParams.exp, '2PDATA', P.LoadParams.Shanknum, [], P.LoadParams.LoadSmthTime, P.LoadParams.Samplerate);
%     EXP.Spk = SpkLoadMultipleExpts(...);
%     obj.data.ephys = false;
%     obj.data.twophoton = true;
elseif isdir([DIRS.multichanspikes filesep P.LoadParams.animal filesep num2str(P.LoadParams.series)])
    EXP.Nav = VRLoadMultipleExpts(P.LoadParams.animal, P.LoadParams.series, P.LoadParams.exp, 'SPIKES', P.LoadParams.Shanknum, P.LoadParams.Shanksuffix, P.LoadParams.LoadSmthTime, P.LoadParams.Samplerate);
%     EXP.Spk = SpkLoadMultipleExpts(...);
%     obj.data.ephys = true;
%     obj.data.twophoton = false;
else
    EXP.Nav = VRLoadMultipleExpts(P.LoadParams.animal, P.LoadParams.series, P.LoadParams.exp, 'BEHAV_ONLY', P.LoadParams.Shanknum, P.LoadParams.Shanksuffix, P.LoadParams.LoadSmthTime, P.LoadParams.Samplerate);
%     obj.data.ephys = false;
%     obj.data.twophoton = false;
end


%if nthetaphsbins > 0, theta phase info must have been saved before by 
%preprocessing LFP with preprocessLFP.m
chnum = 34;%ch from which to measure theta phase. 
obj.data.es = VRbinbythetaphase(obj.data.es,nthetaphsbins,chnum);

% obj.data.es = getESDataSubset(obj.data.es, 'smthBallSpd', speed_th, []);
    
obj.CalculateSubsets();
end