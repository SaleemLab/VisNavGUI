function [EXP, PlotVarDecoder, DecodingDialog] = BatchBayesDecoder(EXP, DIRS, PlotVarDecoder, DecodingDialog)
strlistvarname = {'2p data','electrophys data'};
[varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
if ok && varnamesel == 1
    batch2p = true;
    expt = getExperimentList2p;
    strlistvarname = {'V1medial','V1lateral','PPC', 'AL', 'V1medialV1lateral', 'V1medialV1lateralPPCAL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area', 'SelectionMode', 'single', 'InitialValue', 1);
    area_str = strlistvarname{varnamesel};
elseif ok
    batch2p = false;
    expt = getExperimentList;
    area_str = 'CA1V1';
end
prompt = {'Speed Threshold';'nthetaphsbins';'nspdbins';'neyebins';'nphsbins';'Tsmthwin_field (ms)';'Spatial smth (%)';'Spatial win size (cm)';'Tsmthwin_dec (ms)';'# of X bins';'FGoodcluster';'FMUAcluster';'FoverwriteCellinfo';'latency correction';'# dec. phase bins';'# dec. Speed bins'};
dlg_title = 'Parameters';
num_lines = 1;
def = {'5';'0';'5';'1';'1';'15';'4';'1';'50';'100';'1';'0';'1';'0';'1';'5'};
nruns = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(nruns)
    SpeedThreshold = str2num(nruns{1});
    nthetaphsbins = str2num(nruns{2});
    nspeedbins = str2num(nruns{3});
    neyeXbins = str2num(nruns{4});
    nphsbins = str2num(nruns{5});
    Tsmthwin = str2num(nruns{6});
    Xsmthwin = str2num(nruns{7});
    Xbinsize = str2num(nruns{8});
    Tsmthwin_dec = str2num(nruns{9});
    nDecbins = str2num(nruns{10});
    FGoodcluster = logical(str2num(nruns{11}));
    FMUAcluster = logical(str2num(nruns{12}));
    FoverwriteCellinfo = logical(str2num(nruns{13}));
    latcorrectionlist = str2num(nruns{14});
    nthetaphsbinslist = str2num(nruns{15});
    nspdbinslist = str2num(nruns{16});
end
nDecbins = 100/Xbinsize;

nanimals = numel(expt);
Nperm_cellprop = 0;%100

FUnsortedcluster = 0;
maxRate = inf;
zth = -inf;
SSImin = -inf;

if FGoodcluster && ~FMUAcluster
    cellstr = 'goodonly';%'goodonly_unwrapped';%'goodonly_spktimes';%
elseif FGoodcluster && FMUAcluster
    cellstr = 'All';
end

filesuffix_EXP = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' 'Decwin' num2str(Tsmthwin_dec) '_' 'nDecbins' num2str(nDecbins) '_' num2str(nspeedbins) 'speedbins' '_' num2str(neyeXbins) 'eyebins' '_' num2str(nphsbins) 'thetabins' '_' cellstr];
filesuffix_cellprop = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' num2str(nthetaphsbins) 'thetabins' '_' cellstr];

for ianimal = 7:7%1:numel(expt)
    for iseries = 3:3%1:numel(expt(ianimal).series)
        if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
            EXP = TVRData;
            EXP.animal = expt(ianimal).animal;
            EXP.series = num2str(expt(ianimal).series{iseries});
            EXP.iseries = expt(ianimal).series{iseries};
            EXP.exptList = expt(ianimal).exp{iseries};
            if ~batch2p
                dDIRname = ['D:\DATA\batch'  filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.multichanspikes filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
            else
                dDIRname = ['D:\DATA\batch\2p'  filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];%[DIRS.data2p filesep EXP.animal filesep num2str(EXP.iseries) filesep 'processed'];
            end
            if ~isdir(dDIRname)
                mkdir(dDIRname)
            end
            
            savedfile = [dDIRname filesep 'EXP_' filesuffix_EXP '.mat'];
            savedfile_latency = [dDIRname filesep 'latency_' filesuffix_EXP '.mat'];
            savedfile_Thetaphsbins = [dDIRname filesep 'Thetaphsbins_' filesuffix_EXP '.mat'];
            savedfile_Speedbins = [dDIRname filesep 'Speedbins_' filesuffix_EXP '.mat'];
            disp([EXP.animal ' series ' num2str(EXP.iseries)]);
            
            if strcmp(EXP.animal,'M160114C_BALL') && EXP.iseries == 323
                shanknum = [0:15 0];
                suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
            else
                shanknum = [0:7 0];
                suffix =  {'CA1','CA1','CA1','CA1','CA1','CA1','CA1','CA1','V1'};
            end
            
            smthspeedwin = 150;
            samplerate = 60;
            EXP.LoadVRData(shanknum, suffix, SpeedThreshold, nthetaphsbins,smthspeedwin,samplerate);
            
            savedcellinfo = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_cellProperties.mat'];
            savedmaps1d = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d.mat'];
            savedmaps1d_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_Shf.mat'];
            savedmaps1d_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_2fold.mat'];
            savedmaps1dunwrapped = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dunwrapped.mat'];
            savedmaps1dunwrapped_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dunwrapped_Shf.mat'];
            savedmaps1dunwrapped_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dunwrapped_2fold.mat'];
            savedmaps1d_gain = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_gain.mat'];
            savedmaps1d_gain_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_gain_Shf.mat'];
            savedmaps2d_gain = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_gain.mat'];
            savedmaps2d_gain_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_gain_Shf.mat'];
            savedmaps1dspk = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dspk.mat'];
            savedmaps1dspk_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dspk_Shf.mat'];
            savedmaps1dspk_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dspk_2fold.mat'];
            savedmaps2d_spd = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_spd.mat'];
            savedmaps1d_phs = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_phs.mat'];
            savedmaps1d_phs_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_phs_Shf.mat'];
            savedmaps1d_phs_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1d_phs_2fold.mat'];
            savedmaps1dspk_phs = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dspk_phs.mat'];
            savedmaps1dspk_phs_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dspk_phs_Shf.mat'];
            savedmaps1dspk_phs_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps1dspk_phs_2fold.mat'];
            savedmaps2d_phs = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_phs.mat'];
            savedmaps2d_phs_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_phs_Shf.mat'];
            savedmaps2d_phs_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2d_phs_2fold.mat'];
            savedmaps2dspk_phs = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2dspk_phs.mat'];
            savedmaps2dspk_phs_Shf = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2dspk_phs_Shf.mat'];
            savedmaps2dspk_phs_2fold = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_maps2dspk_phs_2fold.mat'];
            savedVS = [dDIRname filesep EXP.animal '_' num2str(EXP.iseries) '_' filesuffix_cellprop '_VS.mat'];
            
            if exist(savedcellinfo) && ~FoverwriteCellinfo
                EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
                %we call definecellprop here to get the correction factor for the LFPphase
                if isfield(EXP.data.es,'LFPphase')
                    EXP.defineCellProp(Nperm_cellprop);
                    nProbe = numel(unique(EXP.CellInfo.Probe));
                    EXP.data.es.LFPphase2 = NaN(size(EXP.data.es.LFPphase,1),nProbe);
                    for iprobe = 1:nProbe
                        EXP.data.es.LFPphase2(:,iprobe) = mod(EXP.data.es.LFPphase - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);%- EXP.CellInfo.LFP2Spike_phscorrMUA{iprobe},360);
                    end
                    for icell = 1:size(EXP.data.es.spikePhase,2)
                        EXP.data.es.spikePhase(:,icell) = mod(EXP.data.es.spikePhase(:,icell) - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);
                        EXP.data.es.spikePhase(isnan(EXP.data.es.spikePhase(:,icell)),icell) = EXP.data.es.LFPphase2(isnan(EXP.data.es.spikePhase(:,icell)),1);
                    end
                end
                disp('loading cellInfo from saved file');
                S = load(savedcellinfo);
                EXP.CellInfo = S.CellInfo;
                S = load(savedVS);
                EXP.data.VStuning = S.VStuning;
                EXP.data.VSstim = S.VSstim;
                
                S = load(savedmaps1d);
                EXP.maps1d.trajPercent = S.maps1d;
                S = load(savedmaps1d_Shf);
                EXP.maps1d.trajPercent_Shf = S.maps1d;
                S = load(savedmaps1d_2fold);
                EXP.maps1d.trajPercent_2fold = S.maps1d;
                
                S = load(savedmaps1dunwrapped);
                EXP.maps1d.trajPercentunwrapped = S.maps1d;
                S = load(savedmaps1dunwrapped_Shf);
                EXP.maps1d.trajPercentunwrapped_Shf = S.maps1d;
                S = load(savedmaps1dunwrapped_2fold);
                EXP.maps1d.trajPercentunwrapped_2fold = S.maps1d;
                
                S = load(savedmaps1d_gain);
                EXP.maps1d.gain = S.maps1d;
                S = load(savedmaps1d_gain_Shf);
                EXP.maps1d.gain_Shf = S.maps1d;
                
                S = load(savedmaps2d_gain);
                EXP.maps2d.trajPercent_gain = S.maps2d;
                S = load(savedmaps2d_gain_Shf);
                EXP.maps2d.trajPercent_gain_Shf = S.maps2d;
                
%                 S = load(savedmaps1dspk);
%                 EXP.maps1d.spiketrajPercent = S.maps1d;
%                 S = load(savedmaps1dspk_Shf);
%                 EXP.maps1d.spiketrajPercent_Shf = S.maps1d;
%                 S = load(savedmaps1dspk_2fold);
%                 EXP.maps1d.spiketrajPercent_2fold = S.maps1d;
                
%                 S = load(savedmaps2d_spd);
%                 EXP.maps2d.trajPercent_smthBallSpd = S.maps2d;

                S = load(savedmaps1d_phs);
                EXP.maps1d.LFPphase2 = S.maps1d;
                S = load(savedmaps1d_phs_Shf);
                EXP.maps1d.LFPphase2_Shf = S.maps1d;
                S = load(savedmaps1d_phs_2fold);
                EXP.maps1d.LFPphase2_2fold = S.maps1d;
                
%                 S = load(savedmaps1dspk_phs);
%                 EXP.maps1d.spikePhase = S.maps1d;
%                 S = load(savedmaps1dspk_phs_Shf);
%                 EXP.maps1d.spikePhase_Shf = S.maps1d;
%                 S = load(savedmaps1dspk_phs_2fold);
%                 EXP.maps1d.spikePhase_2fold = S.maps1d;
                
                S = load(savedmaps2d_phs);
                EXP.maps2d.trajPercent_LFPphase2 = S.maps2d;
                S = load(savedmaps2d_phs_Shf);
                EXP.maps2d.trajPercent_LFPphase2_Shf = S.maps2d;
                S = load(savedmaps2d_phs_2fold);
                EXP.maps2d.trajPercent_LFPphase2_2fold = S.maps2d;
                
%                 S = load(savedmaps2dspk_phs);
%                 EXP.maps2d.spiketrajPercent_spikePhase = S.maps2d;
%                 S = load(savedmaps2dspk_phs_Shf);
%                 EXP.maps2d.spiketrajPercent_spikePhase_Shf = S.maps2d;
%                 S = load(savedmaps2dspk_phs_2fold);
%                 EXP.maps2d.spiketrajPercent_spikePhase_2fold = S.maps2d;
            else
                EXP.maps1d = [];
                EXP.maps2d = [];
                EXP.defineCellInfo(EXP.data.es.spikeIDs, EXP.data.es.chanIDs,EXP.data.es.ProbeIDs);
                %we call definecellprop here to get the correction factor for the LFPphase
                if isfield(EXP.data.es,'LFPphase')
                    EXP.defineCellProp(Nperm_cellprop);
                    nProbe = numel(unique(EXP.CellInfo.Probe));
                    EXP.data.es.LFPphase2 = NaN(size(EXP.data.es.LFPphase,1),nProbe);
                    for iprobe = 1:nProbe
                        EXP.data.es.LFPphase2(:,iprobe) = mod(EXP.data.es.LFPphase - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);%- EXP.CellInfo.LFP2Spike_phscorrMUA{iprobe},360);
                    end
                    for icell = 1:size(EXP.data.es.spikePhase,2)
                        EXP.data.es.spikePhase(:,icell) = mod(EXP.data.es.spikePhase(:,icell) - round(EXP.CellInfo.LFP2Spike_phscorrMUAnorm{1}),360);
                        EXP.data.es.spikePhase(isnan(EXP.data.es.spikePhase(:,icell)),icell) = EXP.data.es.LFPphase2(isnan(EXP.data.es.spikePhase(:,icell)),1);
                    end
                end
                EXP.CalculateStimTuning([], shanknum, suffix);
                VStuning = EXP.data.VStuning;
                VSstim = EXP.data.VSstim;
                save(savedVS,'VStuning','VSstim');
                delay = 0;
                tic
                EXP.Calculate1Dmaps('trajPercent', Tsmthwin, Xbinsize, Xsmthwin,delay, EXP.data.es.CircularMaze);
                toc
                tic
                EXP.Calculate1Dmaps('trajPercentunwrapped', Tsmthwin, Xbinsize, Xsmthwin,delay, EXP.data.es.CircularMaze);
                toc
                tic
                EXP.Calculate1Dmaps('gain', Tsmthwin, 0.1, 0, delay, false);
                toc
                tic
                EXP.Calculate2Dmaps('trajPercent', 'gain', Tsmthwin, Xbinsize, 0.1, Xsmthwin, 0,delay, EXP.data.es.CircularMaze, false);
                toc
%                 tic
%                 EXP.Calculate1Dmaps('spiketrajPercent', Tsmthwin, Xbinsize, Xsmthwin,delay, EXP.data.es.CircularMaze);
%                 toc

%                 Spdbinsize = 0.1;
%                 SmthSpdWindow = 0.1;
%                 EXP.Calculate2Dmaps('trajPercent', 'smthBallSpd', Tsmthwin, Xbinsize, Spdbinsize, Xsmthwin, SmthSpdWindow, delay, EXP.data.es.CircularMaze, false);
                if isfield(EXP.data.es,'LFPphase2')
                    Phsbinsize = 20;
                    SmthPhsWindow = 40;
                    tic
                    EXP.Calculate1Dmaps('LFPphase2', Tsmthwin, Phsbinsize, SmthPhsWindow, delay, true);
                    toc
%                     tic
%                     EXP.Calculate1Dmaps('spikePhase', Tsmthwin, Phsbinsize, SmthPhsWindow, delay, true);
%                     toc
                    tic
                    EXP.Calculate2Dmaps('trajPercent', 'LFPphase2', Tsmthwin, Xbinsize, Phsbinsize, Xsmthwin, SmthPhsWindow, delay, EXP.data.es.CircularMaze, true);
                    toc
%                     tic
%                     EXP.Calculate2Dmaps('spiketrajPercent', 'spikePhase', Tsmthwin, Xbinsize, Phsbinsize, Xsmthwin, SmthPhsWindow, delay, EXP.data.es.CircularMaze, true);
%                     toc
                end
                try
                maps1d = EXP.maps1d.trajPercent;
                save(savedmaps1d,'maps1d','-v7.3');
                maps1d = EXP.maps1d.trajPercent_Shf;
                save(savedmaps1d_Shf,'maps1d','-v7.3');
                maps1d = EXP.maps1d.trajPercent_2fold;
                save(savedmaps1d_2fold,'maps1d','-v7.3');
                
                maps1d = EXP.maps1d.trajPercentunwrapped;
                save(savedmaps1dunwrapped,'maps1d','-v7.3');
                maps1d = EXP.maps1d.trajPercentunwrapped_Shf;
                save(savedmaps1dunwrapped_Shf,'maps1d','-v7.3');
                maps1d = EXP.maps1d.trajPercentunwrapped_2fold;
                save(savedmaps1dunwrapped_2fold,'maps1d','-v7.3');
%                 maps1d = EXP.maps1d.spiketrajPercent;
%                 save(savedmaps1dspk,'maps1d','-v7.3');
%                 maps1d = EXP.maps1d.spiketrajPercent_Shf;
%                 save(savedmaps1dspk_Shf,'maps1d','-v7.3');
%                 maps1d = EXP.maps1d.spiketrajPercent_2fold;
%                 save(savedmaps1dspk_2fold,'maps1d','-v7.3');
                
%                 maps2d = EXP.maps2d.trajPercent_smthBallSpd;
%                 save(savedmaps2d_spd,'maps2d');

                maps1d = EXP.maps1d.gain;
                save(savedmaps1d_gain,'maps1d','-v7.3');
                maps1d = EXP.maps1d.gain_Shf;
                save(savedmaps1d_gain_Shf,'maps1d','-v7.3');
                
                maps2d = EXP.maps2d.trajPercent_gain;
                save(savedmaps2d_gain,'maps2d','-v7.3');
                maps2d = EXP.maps2d.trajPercent_gain_Shf;
                save(savedmaps2d_gain_Shf,'maps2d','-v7.3');

                maps1d = EXP.maps1d.LFPphase2;
                save(savedmaps1d_phs,'maps1d','-v7.3');
                maps1d = EXP.maps1d.LFPphase2_Shf;
                save(savedmaps1d_phs_Shf,'maps1d','-v7.3');
                maps1d = EXP.maps1d.LFPphase2_2fold;
                save(savedmaps1d_phs_2fold,'maps1d','-v7.3');
%                 maps1d = EXP.maps1d.spikePhase;
%                 save(savedmaps1dspk_phs,'maps1d','-v7.3');
%                 maps1d = EXP.maps1d.spikePhase_Shf;
%                 save(savedmaps1dspk_phs_Shf,'maps1d','-v7.3');
%                 maps1d = EXP.maps1d.spikePhase_2fold;
%                 save(savedmaps1dspk_phs_2fold,'maps1d','-v7.3');
                maps2d = EXP.maps2d.trajPercent_LFPphase2;
                save(savedmaps2d_phs,'maps2d','-v7.3');
                maps2d = EXP.maps2d.trajPercent_LFPphase2_Shf;
                save(savedmaps2d_phs_Shf,'maps2d','-v7.3');
                maps2d = EXP.maps2d.trajPercent_LFPphase2_2fold;
                save(savedmaps2d_phs_2fold,'maps2d','-v7.3');
%                 maps2d = EXP.maps2d.spiketrajPercent_spikePhase;
%                 save(savedmaps2dspk_phs,'maps2d','-v7.3');
%                 maps2d = EXP.maps2d.spiketrajPercent_spikePhase_Shf;
%                 save(savedmaps2dspk_phs_Shf,'maps2d','-v7.3');
%                 maps2d = EXP.maps2d.spiketrajPercent_spikePhase_2fold;
%                 save(savedmaps2dspk_phs_2fold,'maps2d','-v7.3');
                EXP.defineCellProp(Nperm_cellprop);
                CellInfo = EXP.CellInfo;
                save(savedcellinfo,'CellInfo','-v7.3');
                catch
                    keyboard
                end
            end
            
%             delay = 0;
%             Phsbinsize = 20;
%             SmthPhsWindow = 40;
%             EXP.Calculate2Dmaps('trajPercent', 'LFPphase2', Tsmthwin, Xbinsize, Phsbinsize, Xsmthwin, SmthPhsWindow, delay, EXP.data.es.CircularMaze, true);
%             maps2d = EXP.maps2d.trajPercent_LFPphase2;
%             save(savedmaps2d_phs,'maps2d');
%             
%             EXP.defineCellProp(Nperm_cellprop);
%             CellInfo = EXP.CellInfo;
%             save(savedcellinfo,'CellInfo');
            
            %             delay = 800;
            %             EXP.SimulPlaceFields(delay);
            %             EXP.Calculate1Dmaps('trajPercent',PlotVar1DMaps.SmthTimeWindow,delay);
            %             EXP.defineCellProp;
            
            traincont = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
            traingain = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
            trainroomlength = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
            trainoutcome = find(EXP.SubsetVal.outcome == 2);
            type = 'mean';
            Flookuptable = false;
            nruns = 1;
            goodidx_th = 30;
            Tsmth_win = Tsmthwin_dec;%20;%
            Xsmth_win = Xsmthwin/Xbinsize;
            numbins = 100/Xbinsize;
            thetachannel = 34;
            nthetabins = nphsbins;%1;%6;%
            FshuffleThetaPhs = false;
            FshuffleSpdbins = false;
            nspdbins = nspeedbins;%1;%
            neyebins = neyeXbins;%1;%
            Tsmth_field = Tsmthwin;
            speed_th = SpeedThreshold;
            kfold = 20;
            FoptiSmooth = 0;latcorrection = 0;
            alpha = 0;
            delta = 0;
            
            EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'FshuffleThetaPhs', FshuffleThetaPhs, 'FshuffleSpdbins', FshuffleSpdbins);
            
            %             numbinsX = nDecbins;
            %             EXP.RunBayesDecoder2D('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
            %                 'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBinsX', numbinsX, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'smth_spd', smth_spd, 'latcorrection', latcorrection,...
            %                 'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
            %                 'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin);
            
            %             maxtol = 0.1;
            %             EXP = BayesDecoderAverages(EXP,maxtol);%ComputeBayesAverage(EXP,Nperm_bayes);
            
            DecodingDialog.UpdateUIcontrol('Plots','String', PlotVarDecoder.PlotObjList, 'max', numel(PlotVarDecoder.PlotObjList), 'min', 0, 'Val', 1);
            PlotVarDecoder.ChosenObj = PlotVarDecoder.PlotObjList(1);
            DecodingDialog.UpdateUIcontrol('Contrast','String', strsplit(num2str(EXP.SubsetVal.contrast)), 'max', numel(EXP.SubsetVal.contrast), 'min', 0, 'Val', 1:numel(EXP.SubsetVal.contrast));
            contvalidx = find(EXP.SubsetVal.contrast > 0);
            PlotVarDecoder.ChosenContrast = contvalidx(:);
            DecodingDialog.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', 1:numel(EXP.SubsetVal.gain));
            gainvalidx = find(EXP.SubsetVal.gain > 0 & EXP.SubsetVal.gain < 1);
            PlotVarDecoder.ChosenGain = gainvalidx(:)';
            DecodingDialog.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', 1);
            PlotVarDecoder.ChosenRoomlength = (1:EXP.SubsetVal.roomlength);
            DecodingDialog.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', find(EXP.SubsetVal.outcome == 2));
            PlotVarDecoder.ChosenOutcome = find(EXP.SubsetVal.outcome == 2);
            
            save(savedfile,'EXP','-v7.3');
            
            
            if numel(latcorrectionlist)>1
                contval = 0.1:0.05:0.9;
                RLval = 1;
                outcomeVal = 2;
                Err = cell(2,3);
                
                nProbe = numel(EXP.Bayes.PosError0);
                nTimes = numel(EXP.Bayes.Xsmth0{1});
                Xrange = max(EXP.Bayes.Xsmth0{1});
                Prange = size(EXP.Bayes.PosError0{1},2);
                cont_list = find(ismember(EXP.SubsetVal.contrast, contval));
                RL_list = find(ismember(EXP.SubsetVal.roomlength, RLval));
                outcome_list = find(ismember(EXP.SubsetVal.outcome, outcomeVal));
                Xtidx = false(nTimes,Xrange);
                for xx = 1:Xrange
                    Xtidx(:,xx) = EXP.Bayes.Xsmth0{1} == xx;
                end
                Xsmthbin = 2;
                Ysmthbin = 2;
                
                for iprobe = 1:nProbe 
                    for g = [2 1 3]
                        Err{iprobe,g} = NaN(numel(latcorrectionlist),Prange,Xrange);
                    end
                end
                
                for ilat = 1:numel(latcorrectionlist)
                    latcorrection = latcorrectionlist(ilat);
                    EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                        'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                        'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                        'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'FshuffleThetaPhs', FshuffleThetaPhs, 'FshuffleSpdbins', FshuffleSpdbins);
                    
                    for g = [2 1 3]
                        tidx = false(size(EXP.Bayes.Xsmth0{1}));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    ttidx = EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                    tidx(ttidx) = true;
                                end
                            end
                        end
                        if sum(tidx) > 1
                            for iprobe = 1:numel(EXP.Bayes.PosError0)
                                for xx = 1:Xrange
                                    Err{iprobe,g}(ilat,:,xx) = nanmean(EXP.Bayes.PosError0{iprobe}(tidx & Xtidx(:,xx),:));
                                end
                                Err{iprobe,g}(ilat,:,:) = special_smooth2D(squeeze(Err{iprobe,g}(ilat,:,:)), [Ysmthbin/Prange Xsmthbin/Xrange],[EXP.data.es.CircularMaze EXP.data.es.CircularMaze]);
                            end
                        end
                    end
                end
                save(savedfile_latency,'Err','latcorrectionlist','-v7.3');
                latcorrection = 0;
            end
            
            Tsmth_win = Tsmthwin_dec;%20;%
            Xsmth_win = Xsmthwin/Xbinsize;
            numbins = 100/Xbinsize;
            thetachannel = 34;
            nthetabins = nphsbins;%1;%6;%
            FshuffleThetaPhs = false;
            FshuffleSpdbins = false;
            nspdbins = nspeedbins;%1;%
            neyebins = neyeXbins;%1;%
            Tsmth_field = Tsmthwin;
            speed_th = SpeedThreshold;
            kfold = 20;
            FoptiSmooth = 0;latcorrection = 0;
            alpha = 0;
            delta = 0;
            
            if numel(nthetaphsbinslist)>1 || nthetaphsbinslist>1
                contval = 0.1:0.05:0.9;
                RLval = 1;
                outcomeVal = 2;
%                 kfold = 1;
                nbiter_theta = 50;%10;
                
                nProbe = numel(EXP.Bayes.PosError0);
                nTimes = numel(EXP.Bayes.Xsmth0{1});
                Xrange = max(EXP.Bayes.Xsmth0{1});
                Prange = size(EXP.Bayes.PosError0{1},2);
                cont_list = find(ismember(EXP.SubsetVal.contrast, contval));
                RL_list = find(ismember(EXP.SubsetVal.roomlength, RLval));
                outcome_list = find(ismember(EXP.SubsetVal.outcome, outcomeVal));
                Xtidx = false(nTimes,Xrange);
                for xx = 1:Xrange
                    Xtidx(:,xx) = EXP.Bayes.Xsmth0{1} == xx;
                end
                
                for iprobe = 1:nProbe 
                    for g = [2 1 3]
                        maxdecErrX{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        meandecErrX{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        maxdecErrAll{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        meandecErrAll{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        for ibins = 1:numel(nthetaphsbinslist)
                            MaxDecodedPosition{iprobe,g,ibins} = [];
                            MeanDecodedPosition{iprobe,g,ibins} = [];
                            Xsmth0{iprobe,g,ibins} = [];
                            Phsbin{iprobe,g,ibins} = [];
                            Spdbin{iprobe,g,ibins} = [];
                        end
                        
                        maxdecErrX_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),nbiter_theta);
                        meandecErrX_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),nbiter_theta);
                        maxdecErrAll_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),nbiter_theta);
                        meandecErrAll_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),nbiter_theta);
                        for ibins = 1:numel(nthetaphsbinslist)
                            for iter = 1:nbiter_theta
                                MaxDecodedPosition_Shf{iprobe,g,ibins,iter} = [];
                                MeanDecodedPosition_Shf{iprobe,g,ibins,iter} = [];
                                Xsmth0_Shf{iprobe,g,ibins,iter} = [];
                                Phsbin_Shf{iprobe,g,ibins,iter} = [];
                                Spdbin_Shf{iprobe,g,ibins,iter} = [];
                            end
                        end
                    end
                end
                
                for ibins = 1:numel(nthetaphsbinslist)
                    nthetabins = nthetaphsbinslist(ibins);
                    FshuffleThetaPhs = false;
                    EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                        'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                        'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                        'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'FshuffleThetaPhs', FshuffleThetaPhs, 'FshuffleSpdbins', FshuffleSpdbins);
                    
                    for g = [2 1 3]
                        tidx = false(size(EXP.Bayes.Xsmth0{1}));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    ttidx = EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                    tidx(ttidx) = true;
                                end
                            end
                        end
                        
                        if sum(tidx) > 1
                            for iprobe = 1:numel(EXP.Bayes.PosError0)
                                maxdecErr1 = NaN(1,Xrange);
                                meandecErr1 = NaN(1,Xrange);
                                for i = 1:Xrange
                                    maxdecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                    meandecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                    maxdecErr1(i)= Prange/(2*pi)*circ_std(maxdecErr_temp);% nanmean(Prange/(2*pi)*abs(maxdecErr_temp));
                                    meandecErr1(i) = Prange/(2*pi)*circ_std(meandecErr_temp);% nanmean(Prange/(2*pi)*abs(meandecErr_temp));
                                end
                                maxdecErrX{iprobe,g}(ibins) = nanmean(maxdecErr1);
                                meandecErrX{iprobe,g}(ibins) = nanmean(meandecErr1);
                                
                                maxdecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                meandecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                maxdecErrAll{iprobe,g}(ibins) = Prange/(2*pi)*circ_std(maxdecErr2);% nanmean(Prange/(2*pi)*abs(maxdecErr2));
                                meandecErrAll{iprobe,g}(ibins) = Prange/(2*pi)*circ_std(meandecErr2);% nanmean(Prange/(2*pi)*abs(meandecErr2));
                                MaxDecodedPosition{iprobe,g,ibins} = EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx);
                                MeanDecodedPosition{iprobe,g,ibins} = EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx);
                                Phsbin{iprobe,g,ibins} = EXP.Bayes.Phsbin{1}(tidx);
                                Spdbin{iprobe,g,ibins} = EXP.Bayes.Spdbin(tidx);
                                Xsmth0{iprobe,g,ibins} = EXP.Bayes.Xsmth0{1}(tidx);                                
                            end
                        else
                            maxdecErrX{iprobe,g}(ibins) = NaN;
                            meandecErrX{iprobe,g}(ibins) = NaN;
                            maxdecErrAll{iprobe,g}(ibins) = NaN;
                            meandecErrAll{iprobe,g}(ibins) = NaN;
                        end
                    end
                    
                    if nthetabins == 1
                        Niter = 1;
                    else
                        Niter = nbiter_theta;
                    end
                    for iter = 1:Niter
                        FshuffleThetaPhs = true;
                        disp(' ')
                        disp(['Shuffled iteration #' num2str(iter)])
                        disp(' ')
                        EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                            'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                            'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                            'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'FshuffleThetaPhs', FshuffleThetaPhs, 'FshuffleSpdbins', FshuffleSpdbins);
                        
                        for g = [2 1 3]
                            tidx = false(size(EXP.Bayes.Xsmth0{1}));
                            for cont = 1:numel(cont_list)
                                for r = 1:numel(RL_list)
                                    for o = 1:numel(outcome_list)
                                        ttidx = EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                        tidx(ttidx) = true;
                                    end
                                end
                            end
                            
                            if sum(tidx) > 1
                                for iprobe = 1:numel(EXP.Bayes.PosError0)
                                    maxdecErr1 = NaN(1,Xrange);
                                    meandecErr1 = NaN(1,Xrange);
                                    for i = 1:Xrange
                                        maxdecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                        meandecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                        maxdecErr1(i)= Prange/(2*pi)*circ_std(maxdecErr_temp);% nanmean(Prange/(2*pi)*abs(maxdecErr_temp));
                                        meandecErr1(i) = Prange/(2*pi)*circ_std(meandecErr_temp);% nanmean(Prange/(2*pi)*abs(meandecErr_temp));
                                    end
                                    maxdecErrX_Shf{iprobe,g}(ibins,iter) = nanmean(maxdecErr1);
                                    meandecErrX_Shf{iprobe,g}(ibins,iter) = nanmean(meandecErr1);
                                    
                                    maxdecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                    meandecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                    maxdecErrAll_Shf{iprobe,g}(ibins,iter) = Prange/(2*pi)*circ_std(maxdecErr2);% nanmean(Prange/(2*pi)*abs(maxdecErr2));
                                    meandecErrAll_Shf{iprobe,g}(ibins,iter) = Prange/(2*pi)*circ_std(meandecErr2);% nanmean(Prange/(2*pi)*abs(meandecErr2));
                                    MaxDecodedPosition_Shf{iprobe,g,ibins,iter} = EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx);
                                    MeanDecodedPosition_Shf{iprobe,g,ibins,iter} = EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx);
                                    Phsbin_Shf{iprobe,g,ibins} = EXP.Bayes.Phsbin{1}(tidx);
                                    Spdbin_Shf{iprobe,g,ibins} = EXP.Bayes.Spdbin(tidx);
                                    Xsmth0_Shf{iprobe,g,ibins} = EXP.Bayes.Xsmth0{1}(tidx);
                                end
                            else
                                maxdecErrX_Shf{iprobe,g}(ibins,iter) = NaN;
                                meandecErrX_Shf{iprobe,g}(ibins,iter) = NaN;
                                maxdecErrAll_Shf{iprobe,g}(ibins,iter) = NaN;
                                meandecErrAll_Shf{iprobe,g}(ibins,iter) = NaN;
                            end
                        end
                    end
                end
                save(savedfile_Thetaphsbins,'maxdecErrX','meandecErrX','maxdecErrAll','meandecErrAll',...
                                            'maxdecErrX_Shf','meandecErrX_Shf','maxdecErrAll_Shf','meandecErrAll_Shf',...
                                            'MaxDecodedPosition','MeanDecodedPosition','Phsbin','Spdbin','Xsmth0',...
                                            'MaxDecodedPosition_Shf','MeanDecodedPosition_Shf','Phsbin_Shf','Spdbin_Shf','Xsmth0_Shf',...
                                            'nthetaphsbinslist','-v7.3');
                FshuffleThetaPhs = false;
%                 kfold = 20;
            end
            
            Tsmth_win = Tsmthwin_dec;%20;%
            Xsmth_win = Xsmthwin/Xbinsize;
            numbins = 100/Xbinsize;
            thetachannel = 34;
            nthetabins = nphsbins;%1;%6;%
            FshuffleThetaPhs = false;
            FshuffleSpdbins = false;
            nspdbins = nspeedbins;%1;%
            neyebins = neyeXbins;%1;%
            Tsmth_field = Tsmthwin;
            speed_th = SpeedThreshold;
            kfold = 20;
            FoptiSmooth = 0;latcorrection = 0;
            alpha = 0;
            delta = 0;
            
            if numel(nspdbinslist)>1
                contval = 0.1:0.05:0.9;
                RLval = 1;
                outcomeVal = 2;
%                 kfold = 1;
                
                nProbe = numel(EXP.Bayes.PosError0);
                nTimes = numel(EXP.Bayes.Xsmth0{1});
                Xrange = max(EXP.Bayes.Xsmth0{1});
                Prange = size(EXP.Bayes.PosError0{1},2);
                cont_list = find(ismember(EXP.SubsetVal.contrast, contval));
                RL_list = find(ismember(EXP.SubsetVal.roomlength, RLval));
                outcome_list = find(ismember(EXP.SubsetVal.outcome, outcomeVal));
                Xtidx = false(nTimes,Xrange);
                for xx = 1:Xrange
                    Xtidx(:,xx) = EXP.Bayes.Xsmth0{1} == xx;
                end
                
                for iprobe = 1:nProbe 
                    for g = [2 1 3]
                        maxdecErrX{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        meandecErrX{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        maxdecErrAll{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        meandecErrAll{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        for ibins = 1:numel(nthetaphsbinslist)
                            MaxDecodedPosition{iprobe,g,ibins} = [];
                            MeanDecodedPosition{iprobe,g,ibins} = [];
                            Xsmth0{iprobe,g,ibins} = [];
                        end
                        
                        maxdecErrX_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        meandecErrX_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        maxdecErrAll_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        meandecErrAll_Shf{iprobe,g} = NaN(numel(nthetaphsbinslist),1);
                        for ibins = 1:numel(nthetaphsbinslist)
                            MaxDecodedPosition_Shf{iprobe,g,ibins} = [];
                            MeanDecodedPosition_Shf{iprobe,g,ibins} = [];
                            Xsmth0_Shf{iprobe,g,ibins} = [];
                        end
                    end
                end
                
                for ibins = 1:numel(nspdbinslist)
                    nspdbins = nspdbinslist(ibins);
                    FshuffleSpdbins = false;
                    EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                        'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                        'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                        'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'FshuffleThetaPhs', FshuffleThetaPhs, 'FshuffleSpdbins', FshuffleSpdbins);
                    
                    for g = [2 1 3]
                        tidx = false(size(EXP.Bayes.Xsmth0{1}));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    ttidx = EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                    tidx(ttidx) = true;
                                end
                            end
                        end
                        
                        if sum(tidx) > 1
                            for iprobe = 1:numel(EXP.Bayes.PosError0)
                                maxdecErr1 = NaN(1,Xrange);
                                meandecErr1 = NaN(1,Xrange);
                                for i = 1:Xrange
                                    maxdecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                    meandecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                    maxdecErr1(i)= Prange/(2*pi)*circ_std(maxdecErr_temp);% nanmean(Prange/(2*pi)*abs(maxdecErr_temp));
                                    meandecErr1(i) = Prange/(2*pi)*circ_std(meandecErr_temp);% nanmean(Prange/(2*pi)*abs(meandecErr_temp));
                                end
                                maxdecErrX{iprobe,g}(ibins) = nanmean(maxdecErr1);
                                meandecErrX{iprobe,g}(ibins) = nanmean(meandecErr1);
                                
                                maxdecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                meandecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                maxdecErrAll{iprobe,g}(ibins) = Prange/(2*pi)*circ_std(maxdecErr2);% nanmean(Prange/(2*pi)*abs(maxdecErr2));
                                meandecErrAll{iprobe,g}(ibins) = Prange/(2*pi)*circ_std(meandecErr2);% nanmean(Prange/(2*pi)*abs(meandecErr2));
                                MaxDecodedPosition{iprobe,g,ibins} = EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx);
                                MeanDecodedPosition{iprobe,g,ibins} = EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx);
                                Xsmth0{iprobe,g,ibins} = EXP.Bayes.Xsmth0{1}(tidx);                                
                            end
                        else
                            maxdecErrX{iprobe,g}(ibins) = NaN;
                            meandecErrX{iprobe,g}(ibins) = NaN;
                            maxdecErrAll{iprobe,g}(ibins) = NaN;
                            meandecErrAll{iprobe,g}(ibins) = NaN;
                        end
                    end
                    
                    FshuffleSpdbins = true;
                    EXP.RunBayesDecoder('trajPercent','spikeTrain','train_contrast', traincont, 'train_gain', traingain, 'train_roomlength', trainroomlength, 'train_outcome', trainoutcome, 'type', type,'Flookuptable',Flookuptable,...
                        'Tsmth_win', Tsmth_win, 'Xsmth_win', Xsmth_win, 'numBins', numbins, 'nthetaphsbins', nthetabins, 'thetaChannel', thetachannel, 'nspdbins', nspdbins, 'neyebins', neyebins, 'Tsmth_field', Tsmth_field, 'latcorrection', latcorrection,...
                        'alpha', alpha, 'delta', delta, 'nruns', nruns, 'error_th', goodidx_th,'kfold', kfold, 'FoptiSmooth', FoptiSmooth,...
                        'speed_th', speed_th, 'FGoodcluster', FGoodcluster, 'FUnsortedcluster', FUnsortedcluster, 'FMUAcluster', FMUAcluster, 'maxRate', maxRate, 'zth', zth, 'SSImin', SSImin, 'FshuffleThetaPhs', FshuffleThetaPhs, 'FshuffleSpdbins', FshuffleSpdbins);
                    
                    for g = [2 1 3]
                        tidx = false(size(EXP.Bayes.Xsmth0{1}));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    ttidx = EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                    tidx(ttidx) = true;
                                end
                            end
                        end
                        
                        if sum(tidx) > 1
                            for iprobe = 1:numel(EXP.Bayes.PosError0)
                                maxdecErr1 = NaN(1,Xrange);
                                meandecErr1 = NaN(1,Xrange);
                                for i = 1:Xrange
                                    maxdecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                    meandecErr_temp = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx & Xtidx(:,i)),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx & Xtidx(:,i)));
                                    maxdecErr1(i)= Prange/(2*pi)*circ_std(maxdecErr_temp);% nanmean(Prange/(2*pi)*abs(maxdecErr_temp));
                                    meandecErr1(i) = Prange/(2*pi)*circ_std(meandecErr_temp);% nanmean(Prange/(2*pi)*abs(meandecErr_temp));
                                end
                                maxdecErrX_Shf{iprobe,g}(ibins) = nanmean(maxdecErr1);
                                meandecErrX_Shf{iprobe,g}(ibins) = nanmean(meandecErr1);
                                
                                maxdecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                meandecErr2 = circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx),2*pi/Prange*EXP.Bayes.Xsmth0{1}(tidx));
                                maxdecErrAll_Shf{iprobe,g}(ibins) = Prange/(2*pi)*circ_std(maxdecErr2);% nanmean(Prange/(2*pi)*abs(maxdecErr2));
                                meandecErrAll_Shf{iprobe,g}(ibins) = Prange/(2*pi)*circ_std(meandecErr2);% nanmean(Prange/(2*pi)*abs(meandecErr2));
                                MaxDecodedPosition_Shf{iprobe,g,ibins} = EXP.Bayes.MaxDecodedPosition0{iprobe}(tidx);
                                MeanDecodedPosition_Shf{iprobe,g,ibins} = EXP.Bayes.MeanDecodedPosition0{iprobe}(tidx);
                                Xsmth0_Shf{iprobe,g,ibins} = EXP.Bayes.Xsmth0{1}(tidx);
                            end
                        else
                            maxdecErrX_Shf{iprobe,g}(ibins) = NaN;
                            meandecErrX_Shf{iprobe,g}(ibins) = NaN;
                            maxdecErrAll_Shf{iprobe,g}(ibins) = NaN;
                            meandecErrAll_Shf{iprobe,g}(ibins) = NaN;
                        end
                    end
                end
                save(savedfile_Speedbins,'maxdecErrX','meandecErrX','maxdecErrAll','meandecErrAll',...
                                            'maxdecErrX_Shf','meandecErrX_Shf','maxdecErrAll_Shf','meandecErrAll_Shf',...
                                            'MaxDecodedPosition','MeanDecodedPosition','Xsmth0',...
                                            'MaxDecodedPosition_Shf','MeanDecodedPosition_Shf','Xsmth0_Shf',...
                                            'nspdbinslist','-v7.3');
                FshuffleSpdbins = false;
%                 kfold = 20;
            end
        end
    end
end
end

function mat_oo = special_smooth2D(mat_i,win,Fcircular)
% lambdaSmooth = 2;
% nonvalid = isnan(mat_i);
% mat_i(nonvalid) = 1;%1
% G = smooth1D(repmat(mat_i,3,3),lambdaSmooth);
% H = smooth1D(G',lambdaSmooth)';
% mat_o = H(size(mat_i,1)+1:2*size(mat_i,1),size(mat_i,2)+1:2*size(mat_i,2));
% mat_o(nonvalid) = NaN;

mat_oo = mat_i;
for k = 1:size(mat_i,3)
    mat_o = mat_i(:,:,k);
    mat_o(isnan(mat_i(:,:,k))) = 0;
    if Fcircular(1)
        mat_o = repmat(mat_o,[3 1]);
    end
    if Fcircular(2)
        mat_o = repmat(mat_o,[1 3]);
    end
    if sum(isnan(win)) > 0
        if isnan(win(1)) && ~isnan(win(2))
            for i = 1:size(mat_o,1)
                mat_o(i,:) = special_smooth_1d(mat_o(i,:), win(2), [], size(mat_i,2));
            end
        end
        if isnan(win(2)) && ~isnan(win(1))
            for j = 1:size(mat_o,2)
                mat_o(:,j) = special_smooth_1d(mat_o(:,j), win(1), [], size(mat_i,1));
            end
        end
    else
        mat_o = special_smooth_2d(mat_o, win, [], [], [size(mat_i,1) size(mat_i,2)]);
    end
    if Fcircular(1)
        mat_o = mat_o(floor(size(mat_o,1)/3)+1:2*floor(size(mat_o,1)/3),:);
    end
    if Fcircular(2)
        mat_o = mat_o(:,floor(size(mat_o,2)/3)+1:2*floor(size(mat_o,2)/3));
    end
    
    mat_o(isnan(mat_i(:,:,k))) = NaN;
    mat_oo(:,:,k) = mat_o;
end
end