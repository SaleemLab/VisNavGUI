function LoadData_VisNav(EXP, P, S)
P.LoadParams.animal = S.animalname;
P.LoadParams.series = S.series;
P.LoadParams.exp = S.explist;

%load synch files
for iexp = 1:numel(S.explist)
    if ~isempty(S.VR_path{iexp,2})
        Nav.SynchSignal = GetSynchTimes(S.VR_path{iexp,2}, signalType);
        if strcmp(P.LoadParams.SynchRef,'VR')
            SynchTimesRef = Nav.SynchTimes;
        end
    end
    if ~isempty(S.ephys_path{iexp,2})
        Spk.SynchSignal = GetSynchTimes(S.ephys_path{iexp,2}, signalType);
        Lfp.SynchSignal = GetSynchTimes(S.ephys_path{iexp,2}, signalType);
        if strcmp(P.LoadParams.SynchRef,'ephys')
            SynchTimesRef = Spk.SynchTimes;
        end
    end
    if ~isempty(S.vis_path{iexp,2})
        Vis.SynchSignal = GetSynchTimes(S.vis_path{iexp,2}, signalType);
        if strcmp(P.LoadParams.SynchRef,'BonV')
            SynchTimesRef = Vis.SynchTimes;
        end
    end
    if ~isempty(S.eye_path{iexp,2})
        Eye.SynchSignal = GetSynchTimes(S.eye_path{iexp,2}, signalType);
        if strcmp(P.LoadParams.SynchRef,'eye')
            SynchTimesRef = Eye.SynchTimes;
        end
    end
    
    %[Nav, Spk, Vis, Eye with zero correction, clockfactor] = SynchSignals(SynchTimesRef,Nav, Spk, Vis, Eye);
    
    if ~isempty(S.VR_path{iexp,1})
        %Nav = LoadNavData(S.VR_path{iexp,1},sampleTimes,Nav_zerocorrection);
    end
    if ~isempty(S.ephys_path{iexp,1})
        %Spk = getSpikes(S.ephys_path{iexp,1},sampleTimes,Spk_zerocorrection)
        %Lfp = getLFP(Lfp)
    end
    if ~isempty(S.vis_path{iexp,2})
        %Vis = getBonvision(Vis)
    end
    if ~isempty(S.eye_path{iexp,2})
        %Eye = getEyetracking(Eye)
    end
    EXP.Nav = combineTwoexpts(EXP.Nav, Nav);
    EXP.Spk = combineTwoexpts(EXP.Spk, Spk);
    EXP.Lfp = combineTwoexpts(EXP.Lfp, Lfp);
    EXP.Vis = combineTwoexpts(EXP.Vis, Vis);
    EXP.Eye = combineTwoexpts(EXP.Eye, Eye);
end

%Then possibly resample relative to another signal with a function similar
%to VRbinbythetaphase



switch type
    case 'BEHAV_ONLY'
        [~, ~, es] = VRWheelLoad(animal, iseries, expt_list(1), SmthTimeWindow);
        if length(expt_list>1)
            for iexp = 2:length(expt_list)
                [~, ~, esX] = VRWheelLoad(animal, iseries, expt_list(iexp));
                es = combineTwoVRexpts(es, esX);
            end
        end
    case 'SPIKES'
        if isempty(shank_list) || isempty(suffix)
            inp = inputdlg({'Enter the group to be loaded:','Enter the suffix of this group:'}...
                           ,'Enter Group number');
            igroup = str2num(inp{1}); 
            iaddinfo = cell(1,numel(igroup));
            for group_idx = 1:numel(igroup)
                iaddinfo{group_idx} = inp{2};
            end
        else
            igroup = shank_list;
            iaddinfo = suffix;
        end
        es = [];
        for iexp = 1:length(expt_list)
            if isempty(es)
                es = getVRspikes(animal,iseries,expt_list(iexp),100,0,1,0,igroup,iaddinfo,false,SmthTimeWindow,samplerate);
            else
                esX = getVRspikes(animal,iseries,expt_list(iexp),100,0,1,0,igroup,iaddinfo,false,SmthTimeWindow,samplerate);
                es = combineTwoVRexpts(es, esX);
            end
        end
    case '2PDATA'
        es = [];  
        Nplane = 4;%add to initial dialogfor loading the file
        igroup = 1:4;%shank_list;
        irecexp = 0;
        for iexp = 1:length(expt_list)
            fname = [animal '_' num2str(iseries) '_' num2str(expt_list(iexp))];
            dDIRname = [DIRS.data2p filesep animal filesep num2str(iseries)];
            if exist([dDIRname filesep fname '_screenTimes.mat'],'file')
                load([dDIRname filesep fname '_screenTimes']);
                if ~isempty(neuralFrameTimes)
                    irecexp = irecexp + 1;
                    if isempty(es)
                        es = getVR2pdata(animal,iseries,expt_list(iexp),irecexp,1,100,igroup,Nplane,SmthTimeWindow);
                    else
                        esX = getVR2pdata(animal,iseries,expt_list(iexp),irecexp,1,100,igroup,Nplane,SmthTimeWindow);
                        es = combineTwoVRexpts(es, esX);
                    end
                else
                    warning(['No recording during session #' num2str(expt_list(iexp))]);
                end
            else
                warning(['No screentimes during session #' num2str(expt_list(iexp))]);
            end
        end
end

EXP.Nav = VRLoadMultipleExpts(P.LoadParams.animal, P.LoadParams.series, P.LoadParams.exp, 'SPIKES', P.LoadParams.Shanknum, P.LoadParams.Shanksuffix, P.LoadParams.LoadSmthTime, P.LoadParams.Samplerate);

%if nthetaphsbins > 0, theta phase info must have been saved before by 
%preprocessing LFP with preprocessLFP.m
chnum = 34;%ch from which to measure theta phase. 
obj.data.es = VRbinbythetaphase(obj.data.es,nthetaphsbins,chnum);

% obj.data.es = getESDataSubset(obj.data.es, 'smthBallSpd', speed_th, []);
    
obj.CalculateSubsets();
end