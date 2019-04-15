function LoadData_VisNav(EXP, P, S)
% Load data from multiple experiments and synch them.
% Usage : LoadData_VisNav(EXP, P, S)
% Input arguments:
%  - EXP: TStructure object which will eventually contain the data in
%    different fields.
%       EXP.Nav: VR behavior
%       EXP.Spk: spiking data
%       EXP.Vis: BonVision stimulus parameters
%       EXP.Eye: eye tracking data
%  - P: parameter structure. See CreateParamsStructure.
%  - S: filepaths as returned by the FAFF function.
%       Each field of S is n x 2 cell array of path.
%       Rows correspond to path to different experiments. First 
%       column contains path to data files; second column contains path to
%       synchronization file.

%************Saving some metadata in P and EXP for convenience************%
P.LoadParams.animal = S.animalname;
P.LoadParams.series = S.series;
P.LoadParams.exp = S.explist;
EXP.animal = S.animalname;
EXP.series = S.series;
EXP.exp = S.explist;

%***********Synching and loading data from selected experiments*************%
for iexp = 1:numel(S.explist)
    
    %Loading Synchronization signals from the different data types
    if ~isempty(S.Nav_path{iexp,2})
        Nav_SynchSignal = GetSynchSignal(S.Nav_path{iexp,2}, P.LoadParams.SynchType);
        if strcmp(P.LoadParams.SynchSignalRef,'VR')
            SynchTimesRef = Nav_SynchTimes;
        end
    end
    if ~isempty(S.Spk_path{iexp,2})
        Spk_SynchSignal = GetSynchSignal(S.Spk_path{iexp,2}, P.LoadParams.SynchType);
        if strcmp(P.LoadParams.SynchSignalRef,'ephys')
            SynchTimesRef = Spk_SynchTimes;
        end
    end
    if ~isempty(S.Lfp_path{iexp,2})
        Lfp_SynchSignal = GetSynchSignal(S.Lfp_path{iexp,2}, P.LoadParams.SynchType);
        if strcmp(P.LoadParams.SynchSignalRef,'ephys')
            SynchTimesRef = Spk_SynchTimes;
        end
    end
    if ~isempty(S.Vis_path{iexp,2})
        Vis_SynchSignal = GetSynchSignal(S.Vis_path{iexp,2}, P.LoadParams.SynchType);
        if strcmp(P.LoadParams.SynchSignalRef,'BonV')
            SynchTimesRef = Vis_SynchTimes;
        end
    end
    if ~isempty(S.Eye_path{iexp,2})
        Eye_SynchSignal = GetSynchSignal(S.Eye_path{iexp,2}, P.LoadParams.SynchType);
        if strcmp(P.LoadParams.SynchSignalRef,'eye')
            SynchTimesRef = Eye_SynchTimes;
        end
    end
    
    %For each type of signal:
    %Synching signals from different data types to SynchTimesRef
    %Loading resampling and synching the data
    %Concatenating the data from successive experiments
    if ~isempty(S.Nav_path{iexp,1})
        [Nav_ZeroTime, Nav_relativeRate] = SynchSignals(SynchTimesRef, Nav_SynchSignal);
        Nav = LoadNavData(S.Nav_path{iexp,1}, P.LoadParams.LoadSmthTime, sampleTimes, Nav_ZeroTime, Nav_relativeRate);
        if ~isprop(EXP, 'Nav')
            EXP.addprop('Nav');
        end
        EXP.Nav = combineTwoexpts(EXP.Nav, Nav);
    end
    if ~isempty(S.Spk_path{iexp,1})
        [Spk_ZeroTime, Spk_relativeRate] = SynchSignals(SynchTimesRef, Spk_SynchSignal);
        Spk = LoadSpkData(S.Spk_path{iexp,1}, sampleTimes, Spk_ZeroTime, Spk_relativeRate, P.LoadParams.Channels, 'good', {'icell','id'});
        if ~isprop(EXP, 'Spk')
            EXP.addprop('Spk');
        end
        EXP.Spk = combineTwoexpts(EXP.Spk, Spk);
    end
    if ~isempty(S.Lfp_path{iexp,1})
        [Lfp_ZeroTime, Lfp_relativeRate] = SynchSignals(SynchTimesRef, Lfp_SynchSignal);
        Lfp = LoadLfpData(S.Lfp_path{iexp,1}, sampleTimes, Lfp_ZeroTime, Lfp_relativeRate, P.LoadParams.Channels);
        if ~isprop(EXP, 'Lfp')
           EXP.addprop('Lfp');
        end
        EXP.Lfp = combineTwoexpts(EXP.Lfp, Lfp);
    end
    if ~isempty(S.Vis_path{iexp,2})
        [Vis_ZeroTime, Vis_relativeRate] = SynchSignals(SynchTimesRef, Vis_SynchSignal);
        if ~isprop(EXP, 'Vis')
            EXP.addprop('Vis');
        end
        %Vis = getBonvision(Vis)
        %EXP.Vis = combineTwoexpts(EXP.Vis, Vis);
    end
    if ~isempty(S.Eye_path{iexp,2})
        [Eye_ZeroTime, Eye_relativeRate] = SynchSignals(SynchTimesRef, Eye_SynchSignal);
        if ~isprop(EXP, 'Eye')
            EXP.addprop('Eye');
        end
        %Eye = getEyetracking(Eye)
        %EXP.Eye = combineTwoexpts(EXP.Eye, Eye);
    end
end
end