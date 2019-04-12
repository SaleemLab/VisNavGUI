function cellprop = BatchCellInfo(batch2p)
SetDirs;
cellprop = [];

if nargin < 1
    strlistvarname = {'2p data','electrophys data'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
    if ok && varnamesel == 1
        batch2p = true;
    elseif ok
        batch2p = false;
    end
end

if batch2p
    expt = getExperimentList2p;
    datadir = 'D:\DATA\batch\2p';%DIRS.data2p;
    
    cellprop.Tsmthwin = 250;%250;%250;%250;%150;%300;%40;%120;%50
    cellprop.Xsmthwin = 2;%1;%
    cellprop.SpeedThreshold = 1;
    cellprop.nthetaphsbins = 0;%1;%
    cellprop.cellstr = 'goodonly';
    
    suffix = ['Twin' num2str(cellprop.Tsmthwin) '_Xwin' num2str(cellprop.Xsmthwin) '_spdth' num2str(cellprop.SpeedThreshold) '_' num2str(cellprop.nthetaphsbins) 'thetabins_' cellprop.cellstr];
    disp(suffix)
    suffix_cellprop = [suffix '_cellProperties'];
    suffix_map1d = [suffix '_maps1d'];
else
    expt = getExperimentList;
    datadir = 'D:\DATA\batch';%datadir = DIRS.multichanspikes;
    
    cellprop.Tsmthwin = 16;%15;%16;%250;%150;%300;%
    cellprop.Xsmthwin = 4;%2;%1;%
    cellprop.SpeedThreshold = 5;
    cellprop.nthetaphsbins = 0;%1;%
    cellprop.cellstr = 'goodonly';%'goodonly_unwrapped';
    suffix = ['Twin' num2str(cellprop.Tsmthwin) '_Xwin' num2str(cellprop.Xsmthwin) '_spdth' num2str(cellprop.SpeedThreshold) '_' num2str(cellprop.nthetaphsbins) 'thetabins_' cellprop.cellstr];
    disp(suffix)
    suffix_cellprop = [suffix '_cellProperties'];
    suffix_map1d = [suffix '_maps1d'];
    suffix_map2dtheta = [suffix '_maps2d_phs'];
    suffix_lfpspkcoherence = [cellprop.cellstr '_LFPSpikeCoherence'];
end

savedfilename_cellprop = ['D:\DATA\batch\All\SingleCells\cellprop' '_Twin' num2str(cellprop.Tsmthwin) '_Xwin' num2str(cellprop.Xsmthwin) '_spdth' num2str(cellprop.SpeedThreshold) '_' num2str(cellprop.nthetaphsbins) 'thetabins_' cellprop.cellstr '.mat'];


nanimal = numel(expt);

nbXpos = 10; nbYpos = 9; nVStimepts = 61;

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        if exist([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_cellprop '.mat'],'file')
            S = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_cellprop '.mat']);
            Smaps1d = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_map1d '.mat']);
            Smaps2dtheta = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_map2dtheta '.mat']);
            Slfpspkcoherence.resCA1V1 = [];
%             Slfpspkcoherence = load([dDIRname filesep animalname '_' num2str(expt(ianimal).series{iseries}) '_' suffix_lfpspkcoherence '.mat']);
            if ~isfield(cellprop,'animal')
                cellprop.animal = ianimal*ones(1,S.CellInfo.NumCells);
                cellprop.iseries = iseries*ones(1,S.CellInfo.NumCells);
                cellprop.CellList = S.CellInfo.CellList;
                cellprop.Shank = S.CellInfo.Shank;
                cellprop.MUAcluster = S.CellInfo.MUAcluster;
                cellprop.Goodcluster = S.CellInfo.Goodcluster;
                cellprop.Finterneuron = S.CellInfo.Finterneuron;
                cellprop.Unsortedcluster = S.CellInfo.Unsortedcluster;
                if isfield(S.CellInfo,'rateautocorr')
                    cellprop.Spatialmodulation = S.CellInfo.Spatialmodulation;
                    cellprop.burstindex = S.CellInfo.burstindex;
                    cellprop.minFR = S.CellInfo.minFR;
                    cellprop.rateautocorr = S.CellInfo.rateautocorr;
%                     cellprop.min2maxSpkwf = S.CellInfo.min2maxSpkwf;
                    min2maxSpkwf = NaN(size(S.CellInfo.min2maxSpkwf));
                    asymmetrywf = NaN(size(S.CellInfo.min2maxSpkwf));
                    min0 = NaN(size(S.CellInfo.min2maxSpkwf));
                    max1 = NaN(size(S.CellInfo.min2maxSpkwf));
                    max2 = NaN(size(S.CellInfo.min2maxSpkwf));
                    for icell = 1:size(S.CellInfo.waveform,1)
                        [pks,minidx] = findpeaks(-S.CellInfo.waveform(icell,:));%find(S.CellInfo.waveform(icell,:) == min(S.CellInfo.waveform(icell,:)));
                        pks = pks(minidx >= 10 & minidx <= size(S.CellInfo.waveform,2)-10);
                        minidx = minidx(minidx >= 10 & minidx <= size(S.CellInfo.waveform,2)-10);
                        minidx = minidx(pks == max(pks));
                        if minidx >= 10 && minidx <= size(S.CellInfo.waveform,2)-10
                            minidx = minidx(1);
                            mini =  min(S.CellInfo.waveform(icell,minidx));
                            min0(icell) = mini(1);
                            maxi = max(S.CellInfo.waveform(icell,1:(minidx-1)));
                            max1(icell) = maxi(1);
                            [maxi,maxidx] = max(S.CellInfo.waveform(icell,(minidx+1):end));
                            max2(icell) = maxi(1);
                            asymmetrywf(icell) = (max2(icell) - max1(icell))/(max2(icell) + max1(icell));
                            min2maxSpkwf(icell) = maxidx;
                        end
                    end
                    cellprop.min2maxSpkwf = min2maxSpkwf;
                    cellprop.asymmetrywf = asymmetrywf;
                    cellprop.minwf = min0;
                    cellprop.max1wf = max1;
                    cellprop.max2wf = max2;
                    
                    cellprop.spkautocorr = S.CellInfo.spkautocorr;
                end
                cellprop.Probe = S.CellInfo.Probe;
                cellprop.bestchan = S.CellInfo.bestchan;
                if strcmp(expt(ianimal).area{iseries},'V1medial')
                    cellprop.Cellpos2p = 1*ones(1,S.CellInfo.NumCells);
                elseif strcmp(expt(ianimal).area{iseries},'V1lateral')
                    cellprop.Cellpos2p = 2*ones(1,S.CellInfo.NumCells);
                elseif strcmp(expt(ianimal).area{iseries},'PPC')
                    cellprop.Cellpos2p = 3*ones(1,S.CellInfo.NumCells);
                elseif strcmp(expt(ianimal).area{iseries},'AL')
                    cellprop.Cellpos2p = 4*ones(1,S.CellInfo.NumCells);
                else
                    cellprop.Cellpos2p = zeros(1,S.CellInfo.NumCells);
                end
                if ~isempty(S.CellInfo.RFXpos)
                    cellprop.RFXpos = S.CellInfo.RFXpos;
                    cellprop.Xpos = S.CellInfo.Xpos;
                    cellprop.XposZmax = S.CellInfo.XposZmax;
                    cellprop.globalXposrep = S.CellInfo.globalXposrep;
                    [~,xposmax] = max(mean(S.CellInfo.RFXpos(S.CellInfo.XposZmax>2,:),1));
                    cellprop.XposPop = xposmax*ones(1,S.CellInfo.NumCells);

                    cellprop.RFYpos = S.CellInfo.RFYpos;
                    cellprop.Ypos = S.CellInfo.Ypos;
                    cellprop.YposZmax = S.CellInfo.YposZmax;
                    cellprop.globalYposrep = S.CellInfo.globalYposrep;
                    [~,yposmax] = max(mean(S.CellInfo.RFYpos(S.CellInfo.YposZmax>2,:),1));
                    cellprop.YposPop = yposmax*ones(1,S.CellInfo.NumCells);
                    
                    cellprop.VStime = (-20:40)*0.05;
                else
                    cellprop.RFXpos = NaN(S.CellInfo.NumCells,nbXpos);
                    cellprop.Xpos = NaN(1,S.CellInfo.NumCells);
                    cellprop.XposZmax = NaN(1,S.CellInfo.NumCells);
                    cellprop.globalXposrep = NaN(S.CellInfo.NumCells,nVStimepts);
                    cellprop.XposPop = NaN(1,S.CellInfo.NumCells);

                    cellprop.RFYpos = NaN(S.CellInfo.NumCells,nbYpos);
                    cellprop.Ypos = NaN(1,S.CellInfo.NumCells);
                    cellprop.YposZmax = NaN(1,S.CellInfo.NumCells);
                    cellprop.globalYposrep = NaN(S.CellInfo.NumCells,nVStimepts);
                    cellprop.YposPop = NaN(1,S.CellInfo.NumCells);
                    
                    cellprop.VStime = (-20:40)*0.05;
                end
                
                for iprobe = 1:2
                    if ~isempty(S.CellInfo.LFP2Spike_phscorrMUA{iprobe})
                        cellprop.phsfieldMUA{iprobe} = S.CellInfo.phsfieldMUA{iprobe}';
                        cellprop.phsfieldMUASE{iprobe} = S.CellInfo.phsfieldMUASE{iprobe}';
                        cellprop.phsfieldZMUA{iprobe} = S.CellInfo.phsfieldZMUA{iprobe};
                        cellprop.phsfieldPosMUA{iprobe} = S.CellInfo.phsfieldPosMUA{iprobe};
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = S.CellInfo.LFP2Spike_phscorrMUA{iprobe};
                        
                        cellprop.phsfieldMUAnorm{iprobe} = S.CellInfo.phsfieldMUAnorm{iprobe};
                        cellprop.phsfieldPosMUAnorm{iprobe} = S.CellInfo.phsfieldPosMUAnorm{iprobe};
                        phsCOM = getCircularAverage(S.CellInfo.phsfieldMUA{iprobe}(:),0,1)/numel(S.CellInfo.phsfieldMUA{iprobe})*360 - 180;
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = 360/(2*pi)*circ_dist(S.CellInfo.LFP2Spike_phscorrMUAnorm{1}/360*2*pi,phsCOM/360*2*pi)*ones(1,S.CellInfo.NumCells);
                                                
                        cellprop.phsfieldMean{iprobe} = S.CellInfo.phsfieldMean{iprobe};
                        cellprop.phsfieldPosMean{iprobe} = S.CellInfo.phsfieldPosMean{iprobe};
                        cellprop.LFP2Spike_phscorrMean{iprobe} = S.CellInfo.LFP2Spike_phscorrMean{iprobe};
                    else
                        cellprop.phsfieldMUA{iprobe} = NaN(size(S.CellInfo.phsfieldMUA{1}'));
                        cellprop.phsfieldMUASE{iprobe} = NaN(size(S.CellInfo.phsfieldMUASE{1}'));
                        cellprop.phsfieldZMUA{iprobe} = NaN(size(S.CellInfo.phsfieldZMUA{1}));
                        cellprop.phsfieldPosMUA{iprobe} = NaN(size(S.CellInfo.phsfieldPosMUA{1}));
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = NaN(size(S.CellInfo.LFP2Spike_phscorrMUA{1}));
                        
                        cellprop.phsfieldMUAnorm{iprobe} = NaN(size(S.CellInfo.phsfieldMUAnorm{1}));
                        cellprop.phsfieldPosMUAnorm{iprobe} = NaN(size(S.CellInfo.phsfieldPosMUAnorm{1}));
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = NaN(1,S.CellInfo.NumCells);
                        
                        cellprop.phsfieldMean{iprobe} = NaN(size(S.CellInfo.phsfieldMean{1}));
                        cellprop.phsfieldPosMean{iprobe} = NaN(size(S.CellInfo.phsfieldPosMean{1}));
                        cellprop.LFP2Spike_phscorrMean{iprobe} = NaN(size(S.CellInfo.LFP2Spike_phscorrMean{1}));
                    end
                end
                
                allcontidx = size(S.CellInfo.field,1);
                for g = [2 1 3]
                    cellprop.field{g} = S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldSE{g} = S.CellInfo.fieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldShf{g} = S.CellInfo.fieldShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldShfSE{g} = S.CellInfo.fieldShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldZ{g} = S.CellInfo.fieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldPos{g} = S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldPosSE{g} = S.CellInfo.fieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldCOM{g} = S.CellInfo.fieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldAmp{g} = S.CellInfo.fieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldAmpSE{g} = S.CellInfo.fieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldShfAmp{g} = S.CellInfo.fieldShfAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldShfAmpSQ{g,1} = quantile(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.fieldShfAmpSQ{g,2} = quantile(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.fieldShfAmpSQ{g,3} = quantile(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.fieldXcorr{g} = S.CellInfo.fieldXcorr{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.rate{g} = S.CellInfo.rate{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldMax{g} = S.CellInfo.fieldMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldMin{g} = S.CellInfo.fieldMin{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.SSI{g} = S.CellInfo.SSI{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.SpatialInfo{g} = S.CellInfo.SpatialInfo{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.SpatialInfoPerSpike{g} = S.CellInfo.SpatialInfoPerSpike{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    AmpZ = abs(S.CellInfo.fieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2})./S.CellInfo.fieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    Amp = abs(S.CellInfo.fieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    ShfAmp = abs(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(ShfAmp,1);
                    ShfAmpZ = sum(repmat(Amp,[nShf 1])< ShfAmp,1)/nShf;
                    cellprop.fieldAmpZ{g} = AmpZ;
                    cellprop.fieldShfAmpZ{g} = ShfAmpZ;
                    
                    cellprop.rategain{g} = S.CellInfo.rategain{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.rategainSE{g} = S.CellInfo.rategainSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    rategainZ = abs(S.CellInfo.rategain{allcontidx,g,1,S.CellInfo.outcomeVal == 2})./S.CellInfo.rategainSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    rategaindiff = abs(S.CellInfo.rategain{allcontidx,g,1,S.CellInfo.outcomeVal == 2} - S.CellInfo.rategain{allcontidx,2,1,S.CellInfo.outcomeVal == 2});
                    Shfrategaindiff = abs(S.CellInfo.rateShfgainiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}-S.CellInfo.rateShfgainiter{allcontidx,2,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(Shfrategaindiff,1);
                    ShfrategaindiffZ = sum(repmat(rategaindiff,[nShf 1])< Shfrategaindiff,1)/nShf;
                    cellprop.rategainZ{g} = rategainZ;
                    cellprop.rategainShfdiffZ{g} = ShfrategaindiffZ;
                    
                    
                    cellprop.fieldXgainShift{g} = S.CellInfo.fieldShift{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldXgainShiftSE{g} = S.CellInfo.fieldShiftSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    ShiftZ = S.CellInfo.fieldShift{allcontidx,g,1,S.CellInfo.outcomeVal == 2}./S.CellInfo.fieldShiftSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    Shift = S.CellInfo.fieldShift{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    ShfShift = S.CellInfo.fieldShfShiftiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}';
                    nShf = size(ShfShift,1);
                    if g == 1
                        ShfShiftZ = sum(repmat(Shift,[nShf 1])> ShfShift,1)/nShf;
                    else
                        ShfShiftZ = sum(repmat(Shift,[nShf 1])< ShfShift,1)/nShf;
                    end
                    cellprop.fieldXgainShiftZ{g} = ShiftZ;
                    cellprop.fieldShfXgainShiftZ{g} = ShfShiftZ;
                    
                    fieldset1 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    fieldset2 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    cellprop.field_set1{g} = fieldset1;
                    cellprop.field_set2{g} = fieldset2;
                    reliabilityCorr = sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).*...
                                               (fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])),2)./...
                                               sqrt(sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).^2,2).*...
                                               sum((fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])).^2,2));
                    cellprop.reliabilityCorr{g} = reliabilityCorr;
                    
                    cellprop.fieldunwrap{g} = S.CellInfo.unwrapfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldunwrapSE{g} = S.CellInfo.unwrapfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldunwrap_set1{g} = squeeze(S.CellInfo.unwrapfield_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    cellprop.fieldunwrap_set2{g} = squeeze(S.CellInfo.unwrapfield_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    
                    cellprop.fieldunwraphalfdiff{g} = S.CellInfo.unwrapfieldhaldiff{allcontidx,g,1,S.CellInfo.outcomeVal == 2}./max(S.CellInfo.unwrapfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2},[],2)';
                    cellprop.fieldunwraphalfdiffSE{g} = S.CellInfo.unwrapfieldhaldiffSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}./max(S.CellInfo.unwrapfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2},[],2)';
                    halfdiffZ = abs(S.CellInfo.unwrapfieldhaldiff{allcontidx,g,1,S.CellInfo.outcomeVal == 2})./S.CellInfo.unwrapfieldhaldiffSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    halfdiff = abs(S.CellInfo.unwrapfieldhaldiff{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    Shfhalfdiff = abs(S.CellInfo.unwrapfieldShfhaldiffiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(Shfhalfdiff,1);
                    ShfhalfdiffZ = sum(repmat(halfdiff,[nShf 1])< Shfhalfdiff,1)/nShf;
                    cellprop.fieldunwraphalfdiffZ{g} = halfdiffZ;
                    ShfhalfdiffZ(isnan(halfdiffZ)) = NaN;
                    cellprop.fieldShfunwraphalfdiffZ{g} = ShfhalfdiffZ;
                    
                    
%                     pqthreshold = 1;
%                     fieldsize = NaN(1,size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1));
%                     for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                         [~,fieldX] = findfield(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),pqthreshold);
%                         fieldsize(icell) = numel(fieldX);
%                     end
%                     cellprop.fieldsize{g} = fieldsize;
                    cellprop.fieldsize{g} = S.CellInfo.fieldsize{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.fieldsizeSE{g} = S.CellInfo.fieldsizeSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    fieldsize_ahead = NaN(size(S.CellInfo.fieldsize{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    fieldsize_behind = NaN(size(S.CellInfo.fieldsize{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    numBinsX = size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                    for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
                        [~,fieldX] = findfield(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),1);
                        if ~isempty(fieldX)
                            fieldsize_ahead(icell) = numBinsX/(2*pi)*min(circ_dist(2*pi/numBinsX*fieldX,2*pi/numBinsX*repmat(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell),[1 numel(fieldX)])));
                            fieldsize_behind(icell) = numBinsX/(2*pi)*max(circ_dist(2*pi/numBinsX*fieldX,2*pi/numBinsX*repmat(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell),[1 numel(fieldX)])));
                        end
                    end
                    cellprop.fieldsize_ahead{g} = fieldsize_ahead;
                    cellprop.fieldsize_behind{g} = fieldsize_behind;
                    
                    cellprop.field2dXPhstheta{g} = S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXPhsthetaSE{g} = S.CellInfo.field2dXPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXthetapos{g} = S.CellInfo.field2dXthetapos{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXthetaposSE{g} = S.CellInfo.field2dXthetaposSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
%                     cellprop.field2dXthetaposNorm{g} = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
%                     cellprop.field2dXthetaposSENorm{g} = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dPhstheta{g} = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dPhsthetaSE{g} = S.CellInfo.field2dPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dXcorrtheta{g} = S.CellInfo.field2dXcorrtheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaSE{g} = S.CellInfo.field2dXcorrthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamax{g} = S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxSE{g} = S.CellInfo.field2dXcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxAmp{g} = S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxAmpSE{g} = S.CellInfo.field2dXcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxsinAmp{g} = S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxsinAmpSE{g} = S.CellInfo.field2dXcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxOffset{g} = S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetamaxOffsetSE{g} = S.CellInfo.field2dXcorrthetamaxOffsetSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dXcorrthetaShf{g} = S.CellInfo.field2dXcorrthetaShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaShfSE{g} = S.CellInfo.field2dXcorrthetaShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaShfmax{g} = S.CellInfo.field2dXcorrthetaShfmax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaShfmaxSE{g} = S.CellInfo.field2dXcorrthetaShfmaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaShfmaxAmp{g} = S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaShfmaxAmpSQ{g,1} = quantile(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dXcorrthetaShfmaxAmpSQ{g,2} = quantile(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dXcorrthetaShfmaxAmpSQ{g,3} = quantile(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.field2dXcorrthetaShfmaxsinAmp{g} = S.CellInfo.field2dXcorrthetaShfmaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,1} = quantile(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,2} = quantile(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,3} = quantile(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.field2dXcorrthetaShfmaxOffset{g} = S.CellInfo.field2dXcorrthetaShfmaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    vecAmp = S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpSE = S.CellInfo.field2dXcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.ZAmpthetapos{g} = abs(vecAmp)./vecAmpSE;
                    vecAmp = S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpSE = S.CellInfo.field2dXcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.ZsinAmpthetapos{g} = abs(vecAmp)./vecAmpSE;
                    vecAmp = abs(S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecAmpShf = abs(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecAmpShf,1);
                    cellprop.ZAmpthetaShfpos{g} = sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf;
                    vecAmp = S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpShf = abs(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecAmpShf,1);
                    cellprop.ZsinAmpthetaShfpos{g} = sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf;
                    
                    cellprop.field2dXRefcorrtheta{g} = S.CellInfo.field2dXRefcorrtheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaSE{g} = S.CellInfo.field2dXRefcorrthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamax{g} = S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxSE{g} = S.CellInfo.field2dXRefcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxAmp{g} = S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxAmpSE{g} = S.CellInfo.field2dXRefcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxsinAmp{g} = S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxsinAmpSE{g} = S.CellInfo.field2dXRefcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxOffset{g} = S.CellInfo.field2dXRefcorrthetamaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetamaxOffsetSE{g} = S.CellInfo.field2dXRefcorrthetamaxOffsetSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dXRefcorrthetaShf{g} = S.CellInfo.field2dXRefcorrthetaShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaShfSE{g} = S.CellInfo.field2dXRefcorrthetaShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaShfmax{g} = S.CellInfo.field2dXRefcorrthetaShfmax{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaShfmaxSE{g} = S.CellInfo.field2dXRefcorrthetaShfmaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaShfmaxAmp{g} = S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,1} = quantile(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,2} = quantile(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,3} = quantile(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.field2dXRefcorrthetaShfmaxsinAmp{g} = S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,1} = quantile(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,2} = quantile(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,3} = quantile(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.field2dXRefcorrthetaShfmaxOffset{g} = S.CellInfo.field2dXRefcorrthetaShfmaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    vecAmp = S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpSE = S.CellInfo.field2dXRefcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.ZAmpRefthetapos{g} = abs(vecAmp)./vecAmpSE;
                    vecAmp = S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpSE = S.CellInfo.field2dXRefcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.ZsinAmpRefthetapos{g} = abs(vecAmp)./vecAmpSE;
                    vecAmp = abs(S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecAmpShf = abs(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecAmpShf,1);
                    cellprop.ZAmpRefthetaShfpos{g} = sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf;
                    vecAmp = abs(S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecAmpShf = abs(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecAmpShf,1);
                    cellprop.ZsinAmpRefthetaShfpos{g} = sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf;
                    
                    cellprop.field2dslopeXY{g} = S.CellInfo.field2dslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dphi0XY{g} = S.CellInfo.field2dphi0XY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2drhoXY{g} = S.CellInfo.field2drhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dslopeXYSE{g} = S.CellInfo.field2dslopeXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dphi0XYSE{g} = S.CellInfo.field2dphi0XYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2drhoXYSE{g} = S.CellInfo.field2drhoXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    cellprop.field2dShfslopeXY{g} = S.CellInfo.field2dShfslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dShfphi0XY{g} = S.CellInfo.field2dShfphi0XY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dShfrhoXY{g} = S.CellInfo.field2dShfrhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dShfslopeXYSQ{g,1} = quantile(S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dShfslopeXYSQ{g,2} = quantile(S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dShfslopeXYSQ{g,3} = quantile(S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.field2dShfphi0XYSQ{g,1} = quantile(S.CellInfo.field2dShfphi0XYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dShfphi0XYSQ{g,2} = quantile(S.CellInfo.field2dShfphi0XYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dShfphi0XYSQ{g,3} = quantile(S.CellInfo.field2dShfphi0XYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.field2dShfrhoXYSQ{g,1} = quantile(S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.field2dShfrhoXYSQ{g,2} = quantile(S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.field2dShfrhoXYSQ{g,3} = quantile(S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    
                    vec = S.CellInfo.field2dslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecSE = S.CellInfo.field2dslopeXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2dslopeXYZ{g} = abs(vec)./vecSE;
                    vec = S.CellInfo.field2drhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecSE = S.CellInfo.field2drhoXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.field2drhoXYZ{g} = abs(vec)./vecSE;
                    vec = (S.CellInfo.field2dslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecShf = (S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecShf,1);
                    cellprop.field2dShfslopeXYZ{g} = min(sum(repmat(vec,[nShf 1]) < vecShf,1),sum(repmat(vec,[nShf 1]) > vecShf,1))/nShf;
                    vec = (S.CellInfo.field2drhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecShf = (S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecShf,1);
                    cellprop.field2dShfrhoXYZ{g} = min(sum(repmat(vec,[nShf 1]) < vecShf,1),sum(repmat(vec,[nShf 1]) > vecShf,1))/nShf;
                    
                    cellprop.field2dXcorrthetamax_set1{g} = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    cellprop.field2dXcorrthetamax_set2{g} = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    thetafieldset1 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                    thetafieldset2 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                    cellprop.field2dXcorrthetamax_set1{g} = thetafieldset1;
                    cellprop.field2dXcorrthetamax_set2{g} = thetafieldset2;
                    thetareliabilityCorr = sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).*...
                        (thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])),2)./...
                        sqrt(sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).^2,2).*...
                        sum((thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])).^2,2));
                    cellprop.thetareliabilityCorr{g} = thetareliabilityCorr;
                    
                    cellprop.phsfieldPos{g} = S.CellInfo.phsfieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldPosSE{g} = S.CellInfo.phsfieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldCOM{g} = S.CellInfo.phsfieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldCOMSE{g} = S.CellInfo.phsfieldCOMSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfield{g} = S.CellInfo.phsfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldSE{g} = S.CellInfo.phsfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldAmp{g} = S.CellInfo.phsfieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldAmpSE{g} = S.CellInfo.phsfieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldsinAmp{g} = S.CellInfo.phsfieldsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldsinAmpSE{g} = S.CellInfo.phsfieldsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldsinOffset{g} = S.CellInfo.phsfieldsinOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldsinOffsetSE{g} = S.CellInfo.phsfieldsinOffsetSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldShf{g} = S.CellInfo.phsfieldShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldShfSE{g} = S.CellInfo.phsfieldShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldShfAmp{g} = S.CellInfo.phsfieldShfAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldShfAmpSQ{g,1} = quantile(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.phsfieldShfAmpSQ{g,2} = quantile(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.phsfieldShfAmpSQ{g,3} = quantile(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    cellprop.phsfieldShfsinAmp{g} = S.CellInfo.phsfieldShfsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldShfsinAmpSQ{g,1} = quantile(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95);
                    cellprop.phsfieldShfsinAmpSQ{g,2} = quantile(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975);
                    cellprop.phsfieldShfsinAmpSQ{g,3} = quantile(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99);
                    
%                     cellprop.phsfieldZ{g} = S.CellInfo.phsfieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsmodulation{g} = S.CellInfo.phsmodulation{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    vecAmp = S.CellInfo.phsfieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpSE = S.CellInfo.phsfieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldAmpZ{g} = abs(vecAmp)./vecAmpSE;
                    vecAmp = S.CellInfo.phsfieldsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    vecAmpSE = S.CellInfo.phsfieldsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.phsfieldsinAmpZ{g} = abs(vecAmp)./vecAmpSE;
                    vecAmp = abs(S.CellInfo.phsfieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecAmpShf = abs(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecAmpShf,1);
                    cellprop.phsfieldShfAmpZ{g} = sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf;
                    vecAmp = abs(S.CellInfo.phsfieldsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                    vecAmpShf = abs(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                    nShf = size(vecAmpShf,1);
                    cellprop.phsfieldShfsinAmpZ{g} = sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf;
                    
                    ZXposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    ZXposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    XposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    XposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    ZXPhsthetaMax = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    ZXPhsthetaMin = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    PhsthetaMod = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXmin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXmax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXahead = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    phaseXbehind = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                    nphsbins = size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                    for icell = 1:size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
                        Xpostheta = S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);%S.CellInfo.field2dXRefcorrthetamaxNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        XposthetaSE = S.CellInfo.field2dXRefcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);%S.CellInfo.field2dXRefcorrthetamaxSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        [~, imax] = max(Xpostheta);
                        ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
                        [~, imin] = min(Xpostheta);
                        ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
                        
                        Xpostheta = S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);%
                        [XposthetaMax(icell), ~] = max(Xpostheta);
                        [XposthetaMin(icell), ~] = min(Xpostheta);
                        
                        Phstheta = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        PhsthetaSE = S.CellInfo.field2dPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                        [~, imax] = max(Phstheta);
                        ZXPhsthetaMax(icell) = abs((Phstheta(imax) - mean(Phstheta))./PhsthetaSE(imax));
                        [~, imin] = min(Phstheta);
                        ZXPhsthetaMin(icell) = abs((Phstheta(imin) - mean(Phstheta))./PhsthetaSE(imin));
                        
                        Phstheta = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:)./mean(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));
                        [PhsthetaMod(icell), ~] = max(Phstheta);
                        
                        Xpostheta = squeeze(S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));%
                        Xpostheta = Xpostheta - mean(Xpostheta);
                        imax = getCircularAverage(Xpostheta',0,0.01,0.05);%[~, imax] = max(Xpostheta);
                        phaseXmax(icell) = imax;
                        imin = getCircularAverage(-Xpostheta',0,0.01,0.05);%[~, imin] = min(Xpostheta);
                        phaseXmin(icell) = imin;
                        
                        Xpostheta_ahead = Xpostheta - mean(Xpostheta);
                        Xpostheta_ahead(Xpostheta_ahead>0) = 0;
                        phaseXahead(icell) = getCircularAverage(abs(Xpostheta_ahead)',0,1);
                        Xpostheta_behind = Xpostheta - mean(Xpostheta);
                        Xpostheta_behind(Xpostheta_behind<0) = 0;
                        phaseXbehind(icell) = getCircularAverage(abs(Xpostheta_behind)',0,1);
                    end
                    cellprop.ZXmaxthetapos{g} = ZXposthetaMax;
                    cellprop.ZXminthetapos{g} = ZXposthetaMin;
                    cellprop.Xmaxthetapos{g} = XposthetaMax;
                    cellprop.Xminthetapos{g} = XposthetaMin;
                    cellprop.ZPhsmaxtheta{g} = ZXPhsthetaMax;
                    cellprop.ZPhsmintheta{g} = ZXPhsthetaMin;
                    cellprop.PhsthetaMod{g} = PhsthetaMod;
                    cellprop.PhsXmaxthetapos{g} = phaseXmax;
                    cellprop.PhsXminthetapos{g} = phaseXmin;
                    cellprop.PhsXaheadthetapos{g} = phaseXahead;
                    cellprop.PhsXbehindthetapos{g} = phaseXbehind;
                    
%                     ZXposthetaMax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     ZXposthetaMin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     phaseXmin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     phaseXmax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                     for icell = 1:size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                         Xpostheta = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                         XposthetaSE = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                         [~, imax] = max(Xpostheta);
%                         phaseXmax(icell) = imax;
%                         ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
%                         [~, imin] = min(Xpostheta);
%                         ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
%                         phaseXmin(icell) = imin;
%                     end
%                     cellprop.ZXmaxthetaposNorm{g} = ZXposthetaMax;
%                     cellprop.ZXminthetaposNorm{g} = ZXposthetaMin;
%                     cellprop.PhsXmaxthetaposNorm{g} = phaseXmax;
%                     cellprop.PhsXminthetaposNorm{g} = phaseXmin;
                    
                    cellprop.PhaseRho{g} = S.CellInfo.PhaseRho{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.PhasePval{g} = S.CellInfo.PhasePval{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.PhaseRayleighPval{g} = S.CellInfo.PhaseRayleighPval{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    cellprop.PhaseRayleighZ{g} = S.CellInfo.PhaseRayleighZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                    
                    
%                     maxpos = [];
%                     for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                         map = repmat(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),[1 3]);
%                         xorig = 0.5:size(map,2)-0.5;
%                         xinterp = 0.1:0.1:size(map,2);
%                         mapinterp = interp1(xorig,map,xinterp,'spline');
%                         mapinterp = mapinterp((numel(mapinterp)/3+1):(2*numel(mapinterp)/3));
%                         [~,imax] = max(mapinterp);
%                         maxpos(icell) = xinterp(imax);
%                     end
%                     cellprop.fieldPos{g} = maxpos(:)';

                    maxcorr{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    maxcorrstd{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    earlyphs = 1:9;
                    latephs = 10:18;
                    maxcorr_early{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    maxcorr_late{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                    if ~isempty(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model)
                        for icell = 1:numel(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning)
                            map = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                            mapbase = Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                            
                            map = map - nanmean(map);
                            map = map./sqrt(sum(map.^2));
                            mapbase = mapbase - nanmean(mapbase);
                            mapbase = mapbase./sqrt(sum(mapbase.^2));
                            fieldXcorr_all = zeros(1,numel(mapbase));
                            xshiftlim = floor(numel(mapbase)/2);
                            outcorrXrange = [1:floor(numel(mapbase)/4) floor(numel(mapbase)*3/4)+1:numel(mapbase)];
                            ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                            for xshift = -xshiftlim:xshiftlim-1
                                ishift = ishift + 1;
                                fieldXcorr_all(ishift) = map*circshift(mapbase,xshift)';
                            end
                            fieldXcorr_all(outcorrXrange) = 0;
                            
                            maxcorr_all = getCircularAverage(fieldXcorr_all(:),0,0.1,0.05);
                            maxcorr{g}(icell) = maxcorr_all;
                            stdmaxcorr = 0;
                            kfold = size(Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel,1);
                            %SEM computed using Jacknife method
                            for i = 1:kfold
                                map_iter = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%map;%
                                mapbase_iter = mapbase;%Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%
                                
                                map_iter = map_iter - nanmean(map_iter);
                                map_iter = map_iter./sqrt(sum(map_iter.^2));
                                mapbase_iter = mapbase_iter - nanmean(mapbase_iter);
                                mapbase_iter = mapbase_iter./sqrt(sum(mapbase_iter.^2));
                                fieldXcorr_iter = zeros(1,numel(mapbase_iter));
                                xshiftlim = floor(numel(mapbase_iter)/2);
                                outcorrXrange = [1:floor(numel(mapbase)/4) floor(numel(mapbase)*3/4)+1:numel(mapbase)];
                                ishift = 0;%floor(numel(mapbase_iter)/2)-xshiftlim-1;
                                for xshift = -xshiftlim:xshiftlim-1
                                    ishift = ishift + 1;
                                    fieldXcorr_iter(ishift) = map_iter*circshift(mapbase_iter,xshift)';
                                end
                                fieldXcorr_iter(outcorrXrange) = 0;
%                                 if sum(isnan(fieldXcorr_iter)) == 0
%                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp,'spline');
%                                 else
%                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp);
%                                 end
                                maxcorr_iter = getCircularAverage(fieldXcorr_iter(:),0,0.1,0.05);
                                numbins = numel(fieldXcorr_iter);
%                                 [~ , maxcorr_iter] = max(fieldXcorr_iter);
%                                 maxcorr_iter = xinterp(maxcorr_iter) - xinterp(floor(size(fieldXcorr_iter,2)/2) + 1);
                                stdmaxcorr = stdmaxcorr + (kfold - 1)/kfold*(numbins/(2*pi)*circ_dist(2*pi/numbins*maxcorr_iter,2*pi/numbins*maxcorr_all)).^2;
                            end
                            stdmaxcorr = (stdmaxcorr).^0.5;
                            maxcorrstd{g}(icell) = stdmaxcorr;
                            
                            %doing the same for early and late theta phase
                            %now
                            numBinsY = size(S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                            map_centered = squeeze(S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:,:));
                            phsoffset = S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell);
                            if isnan(phsoffset)
                                phsoffset = 0;
                            end
                            map_centered = circshift(map_centered,numBinsY/2 - floor(phsoffset/360*numBinsY));
                            mapearly = mean(map_centered(earlyphs,:),1);
                            mapbase_centered = squeeze(S.CellInfo.field2dXPhstheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(icell,:,:));
                            phsoffset = S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(icell);
                            if isnan(phsoffset)
                                phsoffset = 0;
                            end
                            mapbase_centered = circshift(mapbase_centered,numBinsY/2 - floor(phsoffset/360*numBinsY));
                            mapearlybase = mean(mapbase_centered(earlyphs,:),1);
                            
                            mapearly = mapearly - nanmean(mapearly);
                            mapearly = mapearly./sqrt(sum(mapearly.^2));
                            mapearlybase = mapearlybase - nanmean(mapearlybase);
                            mapearlybase = mapearlybase./sqrt(sum(mapearlybase.^2));
                            fieldXcorr_allearly = zeros(1,numel(mapearlybase));
                            xshiftlim = floor(numel(mapearlybase)/2);
                            outcorrXrange = [1:floor(numel(mapearlybase)/4) floor(numel(mapearlybase)*3/4)+1:numel(mapearlybase)];
                            ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                            for xshift = -xshiftlim:xshiftlim-1
                                ishift = ishift + 1;
                                fieldXcorr_allearly(ishift) = mapearly*circshift(mapearlybase,xshift)';
                            end
                            fieldXcorr_allearly(outcorrXrange) = 0;
                            maxcorr_allearly = getCircularAverage(fieldXcorr_allearly(:),0,0.1,0.05);
                            maxcorr_early{g}(icell) = maxcorr_allearly;
                            
                            maplate = mean(map_centered(latephs,:),1);
                            maplatebase = mean(mapbase_centered(latephs,:),1);
                            
                            maplate = maplate - nanmean(maplate);
                            maplate = maplate./sqrt(sum(maplate.^2));
                            maplatebase = maplatebase - nanmean(maplatebase);
                            maplatebase = maplatebase./sqrt(sum(maplatebase.^2));
                            fieldXcorr_alllate = zeros(1,numel(maplatebase));
                            xshiftlim = floor(numel(maplatebase)/2);
                            outcorrXrange = [1:floor(numel(maplatebase)/4) floor(numel(maplatebase)*3/4)+1:numel(maplatebase)];
                            ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                            for xshift = -xshiftlim:xshiftlim-1
                                ishift = ishift + 1;
                                fieldXcorr_alllate(ishift) = maplate*circshift(maplatebase,xshift)';
                            end
                            fieldXcorr_alllate(outcorrXrange) = 0;
                            maxcorr_alllate = getCircularAverage(fieldXcorr_alllate(:),0,0.1,0.05);
                            maxcorr_late{g}(icell) = maxcorr_alllate;
                            
                        end
                    end
                    cellprop.fieldXcorrMax{g} = maxcorr{g};
                    cellprop.fieldXcorrMaxSE{g} = maxcorrstd{g};
                    cellprop.fieldXcorrMax_early{g} = maxcorr_early{g};
                    cellprop.fieldXcorrMax_late{g} = maxcorr_late{g};
                    
                    if ~isempty(Slfpspkcoherence.resCA1V1)
                        if isempty(Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll)
                            lfpspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll';
                            lfpPhsspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll';
                            lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                        else
                            lfpspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll))';
                            lfpPhsspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_PhsCohSpecChAll))';
                            lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                        end
                        cellprop.lfpSpkCoherence{g} = lfpspkcoherence;
                        cellprop.lfpSpkThetaPhsCoherence{g} = lfpspkthetaPhscoherence;
                    end
                end
            else
                cellprop.animal = [cellprop.animal ianimal*ones(1,S.CellInfo.NumCells)];
                cellprop.iseries = [cellprop.iseries iseries*ones(1,S.CellInfo.NumCells)];
                cellprop.CellList = [cellprop.CellList S.CellInfo.CellList];
                cellprop.Shank = [cellprop.Shank S.CellInfo.Shank];
                cellprop.MUAcluster = [cellprop.MUAcluster S.CellInfo.MUAcluster];
                cellprop.Goodcluster = [cellprop.Goodcluster S.CellInfo.Goodcluster];
                cellprop.Finterneuron = [cellprop.Finterneuron S.CellInfo.Finterneuron];
                cellprop.Unsortedcluster = [cellprop.Unsortedcluster S.CellInfo.Unsortedcluster];
                if isfield(S.CellInfo,'rateautocorr')
                    cellprop.Spatialmodulation = [cellprop.Spatialmodulation S.CellInfo.Spatialmodulation];
                    cellprop.burstindex = [cellprop.burstindex S.CellInfo.burstindex];
                    cellprop.minFR = [cellprop.minFR S.CellInfo.minFR];
                    cellprop.rateautocorr = [cellprop.rateautocorr S.CellInfo.rateautocorr];
%                     cellprop.min2maxSpkwf = [cellprop.min2maxSpkwf S.CellInfo.min2maxSpkwf];
                    min2maxSpkwf = NaN(size(S.CellInfo.min2maxSpkwf));
                    asymmetrywf = NaN(size(S.CellInfo.min2maxSpkwf));
                    min0 = NaN(size(S.CellInfo.min2maxSpkwf));
                    max1 = NaN(size(S.CellInfo.min2maxSpkwf));
                    max2 = NaN(size(S.CellInfo.min2maxSpkwf));
                    for icell = 1:size(S.CellInfo.waveform,1)
                        [pks,minidx] = findpeaks(-S.CellInfo.waveform(icell,:));%find(S.CellInfo.waveform(icell,:) == min(S.CellInfo.waveform(icell,:)));
                        pks = pks(minidx >= 10 & minidx <= size(S.CellInfo.waveform,2)-10);
                        minidx = minidx(minidx >= 10 & minidx <= size(S.CellInfo.waveform,2)-10);
                        minidx = minidx(pks == max(pks));
                        if minidx >= 10 && minidx <= size(S.CellInfo.waveform,2)-10
                            minidx = minidx(1);
                            mini =  min(S.CellInfo.waveform(icell,minidx));
                            min0(icell) = mini(1);
                            maxi = max(S.CellInfo.waveform(icell,1:(minidx-1)));
                            max1(icell) = maxi(1);
                            [maxi,maxidx] = max(S.CellInfo.waveform(icell,(minidx+1):end));
                            max2(icell) = maxi(1);
                            asymmetrywf(icell) = (max2(icell) - max1(icell))/(max2(icell) + max1(icell));
                            min2maxSpkwf(icell) = maxidx;
                        end
                    end
                    cellprop.min2maxSpkwf = [cellprop.min2maxSpkwf min2maxSpkwf];
                    cellprop.asymmetrywf = [cellprop.asymmetrywf  asymmetrywf];
                    cellprop.minwf = [cellprop.minwf  min0];
                    cellprop.max1wf = [cellprop.max1wf  max1];
                    cellprop.max2wf = [cellprop.max2wf  max2];
                    cellprop.spkautocorr = cat(1,cellprop.spkautocorr,S.CellInfo.spkautocorr);
                end
                cellprop.Probe = [cellprop.Probe S.CellInfo.Probe];
                cellprop.bestchan = [cellprop.bestchan S.CellInfo.bestchan];
                if strcmp(expt(ianimal).area{iseries},'V1medial')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 1*ones(1,S.CellInfo.NumCells)];
                elseif strcmp(expt(ianimal).area{iseries},'V1lateral')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 2*ones(1,S.CellInfo.NumCells)];
                elseif strcmp(expt(ianimal).area{iseries},'PPC')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 3*ones(1,S.CellInfo.NumCells)];
                elseif strcmp(expt(ianimal).area{iseries},'AL')
                    cellprop.Cellpos2p = [cellprop.Cellpos2p 4*ones(1,S.CellInfo.NumCells)];
                else
                    cellprop.Cellpos2p = [cellprop.Cellpos2p zeros(1,S.CellInfo.NumCells)];
                end
                if ~isempty(S.CellInfo.RFXpos)
                    cellprop.RFXpos = cat(1,cellprop.RFXpos,S.CellInfo.RFXpos);
                    cellprop.Xpos = [cellprop.Xpos S.CellInfo.Xpos];
                    cellprop.XposZmax = [cellprop.XposZmax S.CellInfo.XposZmax];
                    cellprop.globalXposrep = cat(1,cellprop.globalXposrep,S.CellInfo.globalXposrep);
                    [~,xposmax] = max(mean(S.CellInfo.RFXpos(S.CellInfo.XposZmax>2 & S.CellInfo.Probe == 2,:),1));
                    cellprop.XposPop = [cellprop.XposPop xposmax*ones(1,S.CellInfo.NumCells)];

                    cellprop.RFYpos = cat(1,cellprop.RFYpos,S.CellInfo.RFYpos);
                    cellprop.Ypos = [cellprop.Ypos S.CellInfo.Ypos];
                    cellprop.YposZmax = [cellprop.YposZmax S.CellInfo.YposZmax];
                    cellprop.globalYposrep = cat(1,cellprop.globalYposrep,S.CellInfo.globalYposrep);
                    [~,yposmax] = max(mean(S.CellInfo.RFYpos(S.CellInfo.XposZmax>2 & S.CellInfo.Probe == 2,:),1));
                    cellprop.YposPop = [cellprop.YposPop yposmax*ones(1,S.CellInfo.NumCells)];
                else
                    cellprop.RFXpos = cat(1,cellprop.RFXpos,NaN(S.CellInfo.NumCells,nbXpos));
                    cellprop.Xpos = [cellprop.Xpos NaN(1,S.CellInfo.NumCells)];
                    cellprop.XposZmax = [cellprop.XposZmax NaN(1,S.CellInfo.NumCells)];
                    cellprop.globalXposrep = cat(1,cellprop.globalXposrep,NaN(S.CellInfo.NumCells,nVStimepts));
                    cellprop.XposPop = [cellprop.XposPop NaN(1,S.CellInfo.NumCells)];

                    cellprop.RFYpos = cat(1,cellprop.RFYpos,NaN(S.CellInfo.NumCells,nbYpos));
                    cellprop.Ypos = [cellprop.Ypos NaN(1,S.CellInfo.NumCells)];
                    cellprop.YposZmax = [cellprop.YposZmax NaN(1,S.CellInfo.NumCells)];
                    cellprop.globalYposrep = cat(1,cellprop.globalYposrep,NaN(S.CellInfo.NumCells,nVStimepts));
                    cellprop.YposPop = [cellprop.YposPop NaN(1,S.CellInfo.NumCells)];
                end
                
                for iprobe = 1:2
                    if ~isempty(S.CellInfo.LFP2Spike_phscorrMUA{iprobe})
                        cellprop.phsfieldMUA{iprobe} = cat(1,cellprop.phsfieldMUA{iprobe},S.CellInfo.phsfieldMUA{iprobe}');
                        cellprop.phsfieldMUASE{iprobe} = cat(1,cellprop.phsfieldMUASE{iprobe},S.CellInfo.phsfieldMUASE{iprobe}');
                        cellprop.phsfieldZMUA{iprobe} = cat(1,cellprop.phsfieldZMUA{iprobe},S.CellInfo.phsfieldZMUA{iprobe});
                        cellprop.phsfieldPosMUA{iprobe} = cat(1,cellprop.phsfieldPosMUA{iprobe},S.CellInfo.phsfieldPosMUA{iprobe});
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUA{iprobe},S.CellInfo.LFP2Spike_phscorrMUA{iprobe});
                        
                        cellprop.phsfieldMUAnorm{iprobe} = cat(1,cellprop.phsfieldMUAnorm{iprobe},S.CellInfo.phsfieldMUAnorm{iprobe});
                        cellprop.phsfieldPosMUAnorm{iprobe} = cat(1,cellprop.phsfieldPosMUAnorm{iprobe},S.CellInfo.phsfieldPosMUAnorm{iprobe});
                        phsCOM = getCircularAverage(S.CellInfo.phsfieldMUA{iprobe}(:),0,1)/numel(S.CellInfo.phsfieldMUA{iprobe})*360 - 180;
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = cat(2,cellprop.LFP2Spike_phscorrMUAnorm{iprobe},360/(2*pi)*circ_dist(S.CellInfo.LFP2Spike_phscorrMUAnorm{1}/360*2*pi,phsCOM/360*2*pi)*ones(1,S.CellInfo.NumCells));
                        
%                         cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUAnorm{iprobe},S.CellInfo.LFP2Spike_phscorrMUAnorm{iprobe});
                        
                        cellprop.phsfieldMean{iprobe} = cat(1,cellprop.phsfieldMean{iprobe},S.CellInfo.phsfieldMean{iprobe});
                        cellprop.phsfieldPosMean{iprobe} = cat(1,cellprop.phsfieldPosMean{iprobe},S.CellInfo.phsfieldPosMean{iprobe});
                        cellprop.LFP2Spike_phscorrMean{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMean{iprobe},S.CellInfo.LFP2Spike_phscorrMean{iprobe});
                    else
                        cellprop.phsfieldMUA{iprobe} = cat(2,cellprop.phsfieldMUA{iprobe},NaN(size(S.CellInfo.phsfieldMUA{1}')));
                        cellprop.phsfieldMUASE{iprobe} = cat(2,cellprop.phsfieldMUASE{iprobe},NaN(size(S.CellInfo.phsfieldMUASE{1}')));
                        cellprop.phsfieldZMUA{iprobe} = cat(1,cellprop.phsfieldZMUA{iprobe},NaN(size(S.CellInfo.phsfieldZMUA{1})));
                        cellprop.phsfieldPosMUA{iprobe} = cat(1,cellprop.phsfieldPosMUA{iprobe},NaN(size(S.CellInfo.phsfieldPosMUA{1})));
                        cellprop.LFP2Spike_phscorrMUA{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMUA{iprobe},NaN(size(S.CellInfo.LFP2Spike_phscorrMUA{1})));
                        
                        cellprop.phsfieldMUAnorm{iprobe} = cat(2,cellprop.phsfieldMUAnorm{iprobe},NaN(size(S.CellInfo.phsfieldMUAnorm{1})));
                        cellprop.phsfieldPosMUAnorm{iprobe} = cat(1,cellprop.phsfieldPosMUAnorm{iprobe},NaN(size(S.CellInfo.phsfieldPosMUAnorm{1})));
                        cellprop.LFP2Spike_phscorrMUAnorm{iprobe} = cat(2,cellprop.LFP2Spike_phscorrMUAnorm{iprobe},NaN(1,S.CellInfo.NumCells));
                        
                        cellprop.phsfieldMean{iprobe} = cat(2,cellprop.phsfieldMean{iprobe},NaN(size(S.CellInfo.phsfieldMean{1})));
                        cellprop.phsfieldPosMean{iprobe} = cat(1,cellprop.phsfieldPosMean{iprobe},NaN(size(S.CellInfo.phsfieldPosMean{1})));
                        cellprop.LFP2Spike_phscorrMean{iprobe} = cat(1,cellprop.LFP2Spike_phscorrMean{iprobe},NaN(size(S.CellInfo.LFP2Spike_phscorrMean{1})));
                    end
                end
                
                allcontidx = size(S.CellInfo.field,1);
                for g = 1:3
                    if g <= size(S.CellInfo.field,2)
                        cellprop.field{g} = cat(1,cellprop.field{g}, S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldSE{g} = cat(1,cellprop.fieldSE{g}, S.CellInfo.fieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldShf{g} = cat(1,cellprop.fieldShf{g}, S.CellInfo.fieldShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldShfSE{g} = cat(1,cellprop.fieldShfSE{g}, S.CellInfo.fieldShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldZ{g} = cat(2,cellprop.fieldZ{g}, S.CellInfo.fieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldXcorr{g} = cat(1,cellprop.fieldXcorr{g}, S.CellInfo.fieldXcorr{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldPos{g} = cat(2,cellprop.fieldPos{g}, S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldPosSE{g} = cat(2,cellprop.fieldPosSE{g}, S.CellInfo.fieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldCOM{g} = cat(2,cellprop.fieldCOM{g}, S.CellInfo.fieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldAmp{g} = cat(2,cellprop.fieldAmp{g},S.CellInfo.fieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldAmpSE{g} = cat(2,cellprop.fieldAmpSE{g},S.CellInfo.fieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldShfAmp{g} = cat(2,cellprop.fieldShfAmp{g},S.CellInfo.fieldShfAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldShfAmpSQ{g,1} = cat(2,cellprop.fieldShfAmpSQ{g,1},quantile(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.fieldShfAmpSQ{g,2} = cat(2,cellprop.fieldShfAmpSQ{g,2},quantile(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.fieldShfAmpSQ{g,3} = cat(2,cellprop.fieldShfAmpSQ{g,3},quantile(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.rate{g} = cat(2,cellprop.rate{g}, S.CellInfo.rate{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldMax{g} = cat(2,cellprop.fieldMax{g}, S.CellInfo.fieldMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldMin{g} = cat(2,cellprop.fieldMin{g}, S.CellInfo.fieldMin{allcontidx,g,1,S.CellInfo.outcomeVal == 2});

                        cellprop.SSI{g} = cat(2,cellprop.SSI{g}, S.CellInfo.SSI{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.SpatialInfo{g} = cat(2,cellprop.SpatialInfo{g}, S.CellInfo.SpatialInfo{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.SpatialInfoPerSpike{g} = cat(2,cellprop.SpatialInfoPerSpike{g}, S.CellInfo.SpatialInfoPerSpike{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        AmpZ = abs(S.CellInfo.fieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2})./S.CellInfo.fieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        Amp = abs(S.CellInfo.fieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        ShfAmp = abs(S.CellInfo.fieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(ShfAmp,1);
                        ShfAmpZ = sum(repmat(Amp,[nShf 1])< ShfAmp,1)/nShf;
                        cellprop.fieldAmpZ{g} = cat(2,cellprop.fieldAmpZ{g},AmpZ);
                        cellprop.fieldShfAmpZ{g} = cat(2,cellprop.fieldShfAmpZ{g},ShfAmpZ);
                        
                        cellprop.rategain{g} = cat(2,cellprop.rategain{g},S.CellInfo.rategain{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.rategainSE{g} = cat(2,cellprop.rategainSE{g},S.CellInfo.rategainSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        rategainZ = abs(S.CellInfo.rategain{allcontidx,g,1,S.CellInfo.outcomeVal == 2})./S.CellInfo.rategainSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        rategaindiff = abs(S.CellInfo.rategain{allcontidx,g,1,S.CellInfo.outcomeVal == 2} - S.CellInfo.rategain{allcontidx,2,1,S.CellInfo.outcomeVal == 2});
                        Shfrategaindiff = abs(S.CellInfo.rateShfgainiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}-S.CellInfo.rateShfgainiter{allcontidx,2,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(Shfrategaindiff,1);
                        ShfrategaindiffZ = sum(repmat(rategaindiff,[nShf 1])< Shfrategaindiff,1)/nShf;
                        cellprop.rategainZ{g} = cat(2,cellprop.rategainZ{g},rategainZ);
                        cellprop.rategainShfdiffZ{g} = cat(2,cellprop.rategainShfdiffZ{g},ShfrategaindiffZ);
                        
                        cellprop.fieldXgainShift{g} = cat(2,cellprop.fieldXgainShift{g},S.CellInfo.fieldShift{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldXgainShiftSE{g} = cat(2,cellprop.fieldXgainShiftSE{g},S.CellInfo.fieldShiftSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        ShiftZ = S.CellInfo.fieldShift{allcontidx,g,1,S.CellInfo.outcomeVal == 2}./S.CellInfo.fieldShiftSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        Shift = S.CellInfo.fieldShift{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        ShfShift = S.CellInfo.fieldShfShiftiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}';
                        nShf = size(ShfShift,1);
                        if g == 1
                            ShfShiftZ = sum(repmat(Shift,[nShf 1])> ShfShift,1)/nShf;
                        else
                            ShfShiftZ = sum(repmat(Shift,[nShf 1])< ShfShift,1)/nShf;
                        end
                        cellprop.fieldXgainShiftZ{g} = cat(2,cellprop.fieldXgainShiftZ{g},ShiftZ);
                        cellprop.fieldShfXgainShiftZ{g} = cat(2,cellprop.fieldShfXgainShiftZ{g},ShfShiftZ);
                        
                        fieldset1 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                        fieldset2 = squeeze(S.CellInfo.field_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                        cellprop.field_set1{g} = cat(1,cellprop.field_set1{g}, fieldset1);
                        cellprop.field_set2{g} = cat(1,cellprop.field_set2{g}, fieldset2);
                        reliabilityCorr = sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).*...
                                               (fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])),2)./...
                                               sqrt(sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).^2,2).*...
                                               sum((fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])).^2,2));
                        cellprop.reliabilityCorr{g} = cat(1,cellprop.reliabilityCorr{g},reliabilityCorr);
                        
                        cellprop.fieldunwrap{g} = cat(1,cellprop.fieldunwrap{g},S.CellInfo.unwrapfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldunwrapSE{g} = cat(1,cellprop.fieldunwrapSE{g},S.CellInfo.unwrapfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldunwrap_set1{g} = cat(1,cellprop.fieldunwrap_set1{g},squeeze(S.CellInfo.unwrapfield_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:)));
                        cellprop.fieldunwrap_set2{g} = cat(1,cellprop.fieldunwrap_set2{g},squeeze(S.CellInfo.unwrapfield_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:)));
                        
                        cellprop.fieldunwraphalfdiff{g} = cat(2,cellprop.fieldunwraphalfdiff{g},S.CellInfo.unwrapfieldhaldiff{allcontidx,g,1,S.CellInfo.outcomeVal == 2}./max(S.CellInfo.unwrapfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2},[],2)');
                        cellprop.fieldunwraphalfdiffSE{g} = cat(2,cellprop.fieldunwraphalfdiffSE{g},S.CellInfo.unwrapfieldhaldiffSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}./max(S.CellInfo.unwrapfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2},[],2)');
                        halfdiffZ = abs(S.CellInfo.unwrapfieldhaldiff{allcontidx,g,1,S.CellInfo.outcomeVal == 2})./S.CellInfo.unwrapfieldhaldiffSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        halfdiff = abs(S.CellInfo.unwrapfieldhaldiff{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        Shfhalfdiff = abs(S.CellInfo.unwrapfieldShfhaldiffiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(Shfhalfdiff,1);
                        ShfhalfdiffZ = sum(repmat(halfdiff,[nShf 1])< Shfhalfdiff,1)/nShf;
                        ShfhalfdiffZ(isnan(halfdiffZ)) = NaN;
                        cellprop.fieldunwraphalfdiffZ{g} = cat(2,cellprop.fieldunwraphalfdiffZ{g},halfdiffZ);
                        cellprop.fieldShfunwraphalfdiffZ{g} = cat(2,cellprop.fieldShfunwraphalfdiffZ{g},ShfhalfdiffZ);
                        
%                         pqthreshold = 1;
%                         fieldsize = NaN(1,size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1));
%                         for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                             [~,fieldX] = findfield(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),pqthreshold);
%                             fieldsize(icell) = numel(fieldX);
%                         end
%                         cellprop.fieldsize{g} = cat(2,cellprop.fieldsize{g},fieldsize);
                        cellprop.fieldsize{g} = cat(2,cellprop.fieldsize{g},S.CellInfo.fieldsize{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.fieldsizeSE{g} = cat(2,cellprop.fieldsize{g},S.CellInfo.fieldsizeSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        fieldsize_ahead = NaN(size(S.CellInfo.fieldsize{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        fieldsize_behind = NaN(size(S.CellInfo.fieldsize{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        numBinsX = size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                        for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
                            [~,fieldX] = findfield(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),1);
                            if ~isempty(fieldX)
                                fieldsize_ahead(icell) = numBinsX/(2*pi)*min(circ_dist(2*pi/numBinsX*fieldX,2*pi/numBinsX*repmat(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell),[1 numel(fieldX)])));
                                fieldsize_behind(icell) = numBinsX/(2*pi)*max(circ_dist(2*pi/numBinsX*fieldX,2*pi/numBinsX*repmat(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell),[1 numel(fieldX)])));
                            end
                        end
                        cellprop.fieldsize_ahead{g} = cat(2,cellprop.fieldsize_ahead{g},fieldsize_ahead);
                        cellprop.fieldsize_behind{g} = cat(2,cellprop.fieldsize_behind{g},fieldsize_behind);
                        
%                         cellprop.thetaPhaseMax{g} = cat(2,cellprop.thetaPhaseMax{g},S.CellInfo.PhaseMax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXPhstheta{g} = cat(1,cellprop.field2dXPhstheta{g},S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXPhsthetaSE{g} = cat(1,cellprop.field2dXPhsthetaSE{g},S.CellInfo.field2dXPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXthetapos{g} = cat(1,cellprop.field2dXthetapos{g},S.CellInfo.field2dXthetapos{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXthetaposSE{g} = cat(1,cellprop.field2dXthetaposSE{g},S.CellInfo.field2dXthetaposSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
%                         cellprop.field2dXthetaposNorm{g} = cat(1,cellprop.field2dXthetaposNorm{g},S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
%                         cellprop.field2dXthetaposSENorm{g} = cat(1,cellprop.field2dXthetaposSENorm{g},S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        cellprop.field2dXcorrtheta{g} = cat(1,cellprop.field2dXcorrtheta{g},S.CellInfo.field2dXcorrtheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaSE{g} = cat(1,cellprop.field2dXcorrthetaSE{g},S.CellInfo.field2dXcorrthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamax{g} = cat(1,cellprop.field2dXcorrthetamax{g},S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxSE{g} = cat(1,cellprop.field2dXcorrthetamaxSE{g},S.CellInfo.field2dXcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxAmp{g} = cat(2,cellprop.field2dXcorrthetamaxAmp{g},S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxAmpSE{g} = cat(2,cellprop.field2dXcorrthetamaxAmpSE{g},S.CellInfo.field2dXcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxsinAmp{g} = cat(2,cellprop.field2dXcorrthetamaxsinAmp{g},S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxsinAmpSE{g} = cat(2,cellprop.field2dXcorrthetamaxsinAmpSE{g},S.CellInfo.field2dXcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxOffset{g} = cat(2,cellprop.field2dXcorrthetamaxOffset{g},S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetamaxOffsetSE{g} = cat(2,cellprop.field2dXcorrthetamaxOffsetSE{g},S.CellInfo.field2dXcorrthetamaxOffsetSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        cellprop.field2dXcorrthetaShf{g} = cat(1,cellprop.field2dXcorrthetaShf{g},S.CellInfo.field2dXcorrthetaShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaShfSE{g} = cat(1,cellprop.field2dXcorrthetaShfSE{g},S.CellInfo.field2dXcorrthetaShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaShfmax{g} = cat(1,cellprop.field2dXcorrthetaShfmax{g},S.CellInfo.field2dXcorrthetaShfmax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaShfmaxSE{g} = cat(1,cellprop.field2dXcorrthetaShfmaxSE{g},S.CellInfo.field2dXcorrthetaShfmaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaShfmaxAmp{g} = cat(2,cellprop.field2dXcorrthetaShfmaxAmp{g},S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaShfmaxAmpSQ{g,1} = cat(2,cellprop.field2dXcorrthetaShfmaxAmpSQ{g,1},quantile(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dXcorrthetaShfmaxAmpSQ{g,2} = cat(2,cellprop.field2dXcorrthetaShfmaxAmpSQ{g,2},quantile(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dXcorrthetaShfmaxAmpSQ{g,3} = cat(2,cellprop.field2dXcorrthetaShfmaxAmpSQ{g,3},quantile(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.field2dXcorrthetaShfmaxsinAmp{g} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmp{g},S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,1} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,1},quantile(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,2} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,2},quantile(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,3} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,3},quantile(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.field2dXcorrthetaShfmaxOffset{g} = cat(2,cellprop.field2dXcorrthetaShfmaxOffset{g},S.CellInfo.field2dXcorrthetaShfmaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        vecAmp = S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpSE = S.CellInfo.field2dXcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.ZAmpthetapos{g} = cat(2,cellprop.ZAmpthetapos{g},abs(vecAmp)./vecAmpSE);
                        vecAmp = S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpSE = S.CellInfo.field2dXcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.ZsinAmpthetapos{g} = cat(2,cellprop.ZsinAmpthetapos{g},abs(vecAmp)./vecAmpSE);
                        vecAmp = abs(S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecAmpShf = abs(S.CellInfo.field2dXcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecAmpShf,1);
                        cellprop.ZAmpthetaShfpos{g} = cat(2,cellprop.ZAmpthetaShfpos{g},sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf);
                        vecAmp = S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpShf = abs(S.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecAmpShf,1);
                        cellprop.ZsinAmpthetaShfpos{g} = cat(2,cellprop.ZsinAmpthetaShfpos{g},sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf);
                        
                        cellprop.field2dXRefcorrtheta{g} = cat(1,cellprop.field2dXRefcorrtheta{g},S.CellInfo.field2dXRefcorrtheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaSE{g} = cat(1,cellprop.field2dXRefcorrthetaSE{g},S.CellInfo.field2dXRefcorrthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamax{g} = cat(1,cellprop.field2dXRefcorrthetamax{g},S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxSE{g} = cat(1,cellprop.field2dXRefcorrthetamaxSE{g},S.CellInfo.field2dXRefcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxAmp{g} = cat(2,cellprop.field2dXRefcorrthetamaxAmp{g},S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxAmpSE{g} = cat(2,cellprop.field2dXRefcorrthetamaxAmpSE{g},S.CellInfo.field2dXRefcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxsinAmp{g} = cat(2,cellprop.field2dXRefcorrthetamaxsinAmp{g},S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxsinAmpSE{g} = cat(2,cellprop.field2dXRefcorrthetamaxsinAmpSE{g},S.CellInfo.field2dXRefcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxOffset{g} = cat(2,cellprop.field2dXRefcorrthetamaxOffset{g},S.CellInfo.field2dXRefcorrthetamaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetamaxOffsetSE{g} = cat(2,cellprop.field2dXRefcorrthetamaxOffsetSE{g},S.CellInfo.field2dXRefcorrthetamaxOffsetSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        cellprop.field2dXRefcorrthetaShf{g} = cat(1,cellprop.field2dXRefcorrthetaShf{g},S.CellInfo.field2dXRefcorrthetaShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaShfSE{g} = cat(1,cellprop.field2dXRefcorrthetaShfSE{g},S.CellInfo.field2dXRefcorrthetaShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaShfmax{g} = cat(1,cellprop.field2dXRefcorrthetaShfmax{g},S.CellInfo.field2dXRefcorrthetaShfmax{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaShfmaxSE{g} = cat(1,cellprop.field2dXRefcorrthetaShfmaxSE{g},S.CellInfo.field2dXRefcorrthetaShfmaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaShfmaxAmp{g} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmp{g},S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,1} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,1},quantile(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,2} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,2},quantile(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,3} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,3},quantile(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmp{g} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmp{g},S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,1} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,1},quantile(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,2} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,2},quantile(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,3} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,3},quantile(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.field2dXRefcorrthetaShfmaxOffset{g} = cat(2,cellprop.field2dXRefcorrthetaShfmaxOffset{g},S.CellInfo.field2dXRefcorrthetaShfmaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        vecAmp = S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpSE = S.CellInfo.field2dXRefcorrthetamaxAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.ZAmpRefthetapos{g} = cat(2,cellprop.ZAmpRefthetapos{g},abs(vecAmp)./vecAmpSE);
                        vecAmp = S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpSE = S.CellInfo.field2dXRefcorrthetamaxsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.ZsinAmpRefthetapos{g} = cat(2,cellprop.ZsinAmpRefthetapos{g},abs(vecAmp)./vecAmpSE);
                        vecAmp = abs(S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecAmpShf = abs(S.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecAmpShf,1);
                        cellprop.ZAmpRefthetaShfpos{g} = cat(2,cellprop.ZAmpRefthetaShfpos{g},sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf);
                        vecAmp = abs(S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecAmpShf = abs(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecAmpShf,1);
                        cellprop.ZsinAmpRefthetaShfpos{g} = cat(2,cellprop.ZsinAmpRefthetaShfpos{g},sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf);
                        
                        cellprop.field2dslopeXY{g} = cat(2,cellprop.field2dslopeXY{g},S.CellInfo.field2dslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dphi0XY{g} = cat(2,cellprop.field2dphi0XY{g},S.CellInfo.field2dphi0XY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2drhoXY{g} = cat(2,cellprop.field2drhoXY{g},S.CellInfo.field2drhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dslopeXYSE{g} = cat(2,cellprop.field2dslopeXYSE{g},S.CellInfo.field2dslopeXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dphi0XYSE{g} = cat(2,cellprop.field2dphi0XYSE{g},S.CellInfo.field2dphi0XYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2drhoXYSE{g} = cat(2,cellprop.field2drhoXYSE{g},S.CellInfo.field2drhoXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        cellprop.field2dShfslopeXY{g} = cat(2,cellprop.field2dShfslopeXY{g},S.CellInfo.field2dShfslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dShfphi0XY{g} = cat(2,cellprop.field2dShfphi0XY{g},S.CellInfo.field2dShfphi0XY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dShfrhoXY{g} = cat(2,cellprop.field2dShfrhoXY{g},S.CellInfo.field2dShfrhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dShfslopeXYSQ{g,1} = cat(2,cellprop.field2dShfslopeXYSQ{g,1},quantile(S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dShfslopeXYSQ{g,2} = cat(2,cellprop.field2dShfslopeXYSQ{g,2},quantile(S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dShfslopeXYSQ{g,3} = cat(2,cellprop.field2dShfslopeXYSQ{g,3},quantile(S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.field2dShfphi0XYSQ{g,1} = cat(2,cellprop.field2dShfphi0XYSQ{g,1},quantile(S.CellInfo.field2dShfphi0XYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dShfphi0XYSQ{g,2} = cat(2,cellprop.field2dShfphi0XYSQ{g,2},quantile(S.CellInfo.field2dShfphi0XYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dShfphi0XYSQ{g,3} = cat(2,cellprop.field2dShfphi0XYSQ{g,3},quantile(S.CellInfo.field2dShfphi0XYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.field2dShfrhoXYSQ{g,1} = cat(2,cellprop.field2dShfrhoXYSQ{g,1},quantile(S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.field2dShfrhoXYSQ{g,2} = cat(2,cellprop.field2dShfrhoXYSQ{g,2},quantile(S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.field2dShfrhoXYSQ{g,3} = cat(2,cellprop.field2dShfrhoXYSQ{g,3},quantile(S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        
                        vec = S.CellInfo.field2dslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecSE = S.CellInfo.field2dslopeXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.field2dslopeXYZ{g} = cat(2,cellprop.field2dslopeXYZ{g},abs(vec)./vecSE);
                        vec = S.CellInfo.field2drhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecSE = S.CellInfo.field2drhoXYSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.field2drhoXYZ{g} = cat(2,cellprop.field2drhoXYZ{g},abs(vec)./vecSE);
                        vec = (S.CellInfo.field2dslopeXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecShf = (S.CellInfo.field2dShfslopeXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecShf,1);
                        cellprop.field2dShfslopeXYZ{g} = cat(2,cellprop.field2dShfslopeXYZ{g},min(sum(repmat(vec,[nShf 1]) < vecShf,1),sum(repmat(vec,[nShf 1]) > vecShf,1))/nShf);
                        vec = (S.CellInfo.field2drhoXY{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecShf = (S.CellInfo.field2dShfrhoXYiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecShf,1);
                        cellprop.field2dShfrhoXYZ{g} = cat(2,cellprop.field2dShfrhoXYZ{g},min(sum(repmat(vec,[nShf 1]) < vecShf,1),sum(repmat(vec,[nShf 1]) > vecShf,1))/nShf);
                        
                        thetafieldset1 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,1,:));
                        thetafieldset2 = squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(:,2,:));
                        cellprop.field2dXcorrthetamax_set1{g} = cat(1,cellprop.field2dXcorrthetamax_set1{g},thetafieldset1);
                        cellprop.field2dXcorrthetamax_set2{g} = cat(1,cellprop.field2dXcorrthetamax_set2{g},thetafieldset2);
                        thetareliabilityCorr = sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).*...
                                               (thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])),2)./...
                                               sqrt(sum((thetafieldset1-repmat(nanmean(thetafieldset1,2),[1 size(thetafieldset1,2)])).^2,2).*...
                                               sum((thetafieldset2-repmat(nanmean(thetafieldset2,2),[1 size(thetafieldset2,2)])).^2,2));
                        cellprop.thetareliabilityCorr{g} = cat(1,cellprop.thetareliabilityCorr{g},thetareliabilityCorr);                        
                        
                        cellprop.phsfieldPos{g} = cat(2,cellprop.phsfieldPos{g},S.CellInfo.phsfieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldPosSE{g} = cat(2,cellprop.phsfieldPosSE{g},S.CellInfo.phsfieldPosSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldCOM{g} = cat(2,cellprop.phsfieldCOM{g},S.CellInfo.phsfieldCOM{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldCOMSE{g} = cat(2,cellprop.phsfieldCOMSE{g},S.CellInfo.phsfieldCOMSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfield{g} = cat(1,cellprop.phsfield{g},S.CellInfo.phsfield{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldSE{g} = cat(1,cellprop.phsfieldSE{g},S.CellInfo.phsfieldSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldAmp{g} = cat(2,cellprop.phsfieldAmp{g},S.CellInfo.phsfieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldAmpSE{g} = cat(2,cellprop.phsfieldAmpSE{g},S.CellInfo.phsfieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldsinAmp{g} = cat(2,cellprop.phsfieldsinAmp{g},S.CellInfo.phsfieldsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldsinAmpSE{g} = cat(2,cellprop.phsfieldsinAmpSE{g},S.CellInfo.phsfieldsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldsinOffset{g} = cat(2,cellprop.phsfieldsinOffset{g},S.CellInfo.phsfieldsinOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldsinOffsetSE{g} = cat(2,cellprop.phsfieldsinOffsetSE{g},S.CellInfo.phsfieldsinOffsetSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldShf{g} = cat(1,cellprop.phsfieldShf{g},S.CellInfo.phsfieldShf{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldShfSE{g} = cat(1,cellprop.phsfieldShfSE{g},S.CellInfo.phsfieldShfSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldShfAmp{g} = cat(2,cellprop.phsfieldShfAmp{g},S.CellInfo.phsfieldShfAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldShfAmpSQ{g,1} = cat(2,cellprop.phsfieldShfAmpSQ{g,1},quantile(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.phsfieldShfAmpSQ{g,2} = cat(2,cellprop.phsfieldShfAmpSQ{g,2},quantile(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.phsfieldShfAmpSQ{g,3} = cat(2,cellprop.phsfieldShfAmpSQ{g,3},quantile(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
                        cellprop.phsfieldShfsinAmp{g} = cat(2,cellprop.phsfieldShfsinAmp{g},S.CellInfo.phsfieldShfsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsfieldShfsinAmpSQ{g,1} = cat(2,cellprop.phsfieldShfsinAmpSQ{g,1},quantile(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.95));
                        cellprop.phsfieldShfsinAmpSQ{g,2} = cat(2,cellprop.phsfieldShfsinAmpSQ{g,2},quantile(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.975));
                        cellprop.phsfieldShfsinAmpSQ{g,3} = cat(2,cellprop.phsfieldShfsinAmpSQ{g,3},quantile(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2}',0.99));
%                         cellprop.phsfieldZ{g} = cat(2,cellprop.phsfieldZ{g},S.CellInfo.phsfieldZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.phsmodulation{g} = cat(2,cellprop.phsmodulation{g},S.CellInfo.phsmodulation{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        vecAmp = S.CellInfo.phsfieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpSE = S.CellInfo.phsfieldAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.phsfieldAmpZ{g} = cat(2,cellprop.phsfieldAmpZ{g},abs(vecAmp)./vecAmpSE);
                        vecAmp = S.CellInfo.phsfieldsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        vecAmpSE = S.CellInfo.phsfieldsinAmpSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2};
                        cellprop.phsfieldsinAmpZ{g} = cat(2,cellprop.phsfieldsinAmpZ{g},abs(vecAmp)./vecAmpSE);
                        vecAmp = abs(S.CellInfo.phsfieldAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecAmpShf = abs(S.CellInfo.phsfieldShfAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecAmpShf,1);
                        cellprop.phsfieldShfAmpZ{g} = cat(2,cellprop.phsfieldShfAmpZ{g},sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf);
                        vecAmp = abs(S.CellInfo.phsfieldsinAmp{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        vecAmpShf = abs(S.CellInfo.phsfieldShfsinAmpiter{allcontidx,g,1,S.CellInfo.outcomeVal == 2})';
                        nShf = size(vecAmpShf,1);
                        cellprop.phsfieldShfsinAmpZ{g} = cat(2,cellprop.phsfieldShfsinAmpZ{g},sum(repmat(vecAmp,[nShf 1]) < vecAmpShf,1)/nShf);
                        
                        ZXposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        ZXposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        XposthetaMax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        XposthetaMin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        ZXPhsthetaMax = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        ZXPhsthetaMin = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        PhsthetaMod = zeros(size(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXmin = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXmax = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXahead = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        phaseXbehind = zeros(size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
                        nphsbins = size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                        for icell = 1:size(S.CellInfo.field2dXcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
                            Xpostheta = S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);%S.CellInfo.field2dXRefcorrthetamaxNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            XposthetaSE = S.CellInfo.field2dXRefcorrthetamaxSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);%S.CellInfo.field2dXRefcorrthetamaxSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            [~, imax] = max(Xpostheta);
                            ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
                            [~, imin] = min(Xpostheta);
                            ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
                            
                            Xpostheta = S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);%
                            [XposthetaMax(icell), ~] = max(Xpostheta);
                            [XposthetaMin(icell), ~] = min(Xpostheta);
                            
                            Phstheta = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            PhsthetaSE = S.CellInfo.field2dPhsthetaSE{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
                            [~, imax] = max(Phstheta);
                            ZXPhsthetaMax(icell) = abs((Phstheta(imax) - mean(Phstheta))./PhsthetaSE(imax));
                            [~, imin] = min(Phstheta);
                            ZXPhsthetaMin(icell) = abs((Phstheta(imin) - mean(Phstheta))./PhsthetaSE(imin));
                            
                            Phstheta = S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:)./mean(S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));
                            [PhsthetaMod(icell), ~] = max(Phstheta);
                            
                            Xpostheta = squeeze(S.CellInfo.field2dXRefcorrthetamax{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:));%
                            Xpostheta = Xpostheta - mean(Xpostheta);
                            imax = getCircularAverage(Xpostheta',0,0.01,0.05);%[~, imax] = max(Xpostheta);
                            phaseXmax(icell) = imax;
                            imin = getCircularAverage(-Xpostheta',0,0.01,0.05);%[~, imin] = min(Xpostheta);
                            phaseXmin(icell) = imin;
                            
                            Xpostheta_ahead = Xpostheta - mean(Xpostheta);
                            Xpostheta_ahead(Xpostheta_ahead>0) = 0;
                            phaseXahead(icell) = getCircularAverage(abs(Xpostheta_ahead)',0,1);
                            Xpostheta_behind = Xpostheta - mean(Xpostheta);
                            Xpostheta_behind(Xpostheta_behind<0) = 0;
                            phaseXbehind(icell) = getCircularAverage(abs(Xpostheta_behind)',0,1);
                        end
                        cellprop.ZXmaxthetapos{g} = cat(1,cellprop.ZXmaxthetapos{g},ZXposthetaMax);
                        cellprop.ZXminthetapos{g} = cat(1,cellprop.ZXminthetapos{g},ZXposthetaMin);
                        cellprop.Xmaxthetapos{g} = cat(1,cellprop.Xmaxthetapos{g},XposthetaMax);
                        cellprop.Xminthetapos{g} = cat(1,cellprop.Xminthetapos{g},XposthetaMin);
                        cellprop.ZPhsmaxtheta{g} = cat(1,cellprop.ZPhsmaxtheta{g},ZXPhsthetaMax);
                        cellprop.ZPhsmintheta{g} = cat(1,cellprop.ZPhsmintheta{g},ZXPhsthetaMin);
                        cellprop.PhsthetaMod{g} = cat(1,cellprop.PhsthetaMod{g},PhsthetaMod);
                        cellprop.PhsXmaxthetapos{g} = cat(1,cellprop.PhsXmaxthetapos{g},phaseXmax);
                        cellprop.PhsXminthetapos{g} = cat(1,cellprop.PhsXminthetapos{g},phaseXmin);
                        cellprop.PhsXaheadthetapos{g} = cat(1,cellprop.PhsXaheadthetapos{g},phaseXahead);
                        cellprop.PhsXbehindthetapos{g} = cat(1,cellprop.PhsXbehindthetapos{g},phaseXbehind);
                        
%                         ZXposthetaMax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         ZXposthetaMin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         phaseXmin = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         phaseXmax = zeros(size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1),1);
%                         for icell = 1:size(S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                             Xpostheta = S.CellInfo.field2dXthetaposNorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                             XposthetaSE = S.CellInfo.field2dXthetaposSENorm{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:);
%                             [~, imax] = max(Xpostheta);
%                             phaseXmax(icell) = imax;
%                             ZXposthetaMax(icell) = abs((Xpostheta(imax) - mean(Xpostheta))./XposthetaSE(imax));
%                             [~, imin] = min(Xpostheta);
%                             ZXposthetaMin(icell) = abs((Xpostheta(imin) - mean(Xpostheta))./XposthetaSE(imin));
%                             phaseXmin(icell) = imin;
%                         end
%                         cellprop.ZXmaxthetaposNorm{g} = cat(1,cellprop.ZXmaxthetaposNorm{g},ZXposthetaMax);
%                         cellprop.ZXminthetaposNorm{g} = cat(1,cellprop.ZXminthetaposNorm{g},ZXposthetaMin);
%                         cellprop.PhsXmaxthetaposNorm{g} = cat(1,cellprop.PhsXmaxthetaposNorm{g},phaseXmax);
%                         cellprop.PhsXminthetaposNorm{g} =cat(1,cellprop.PhsXminthetaposNorm{g},phaseXmin);
                        
                        cellprop.PhaseRho{g} = cat(2,cellprop.PhaseRho{g},S.CellInfo.PhaseRho{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.PhasePval{g} = cat(2,cellprop.PhasePval{g},S.CellInfo.PhasePval{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.field2dPhstheta{g} = cat(1,cellprop.field2dPhstheta{g},S.CellInfo.field2dPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.PhaseRayleighPval{g} = cat(2,cellprop.PhaseRayleighPval{g},S.CellInfo.PhaseRayleighPval{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        cellprop.PhaseRayleighZ{g} = cat(2,cellprop.PhaseRayleighZ{g},S.CellInfo.PhaseRayleighZ{allcontidx,g,1,S.CellInfo.outcomeVal == 2});
                        
                        
%                         maxpos = [];
%                         for icell = 1:size(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2},1)
%                             map = repmat(S.CellInfo.field{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:),[1 3]);
%                             xorig = 0.5:size(map,2)-0.5;
%                             xinterp = 0.1:0.1:size(map,2);
%                             mapinterp = interp1(xorig,map,xinterp,'spline');
%                             mapinterp = mapinterp((numel(mapinterp)/3+1):(2*numel(mapinterp)/3));
%                             [~,imax] = max(mapinterp);
%                             maxpos(icell) = xinterp(imax);
%                         end
%                         cellprop.fieldPos{g} = cat(2,cellprop.fieldPos{g},maxpos(:)');
                        
                        maxcorr{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        maxcorrstd{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        earlyphs = 1:9;
                        latephs = 10:18;
                        maxcorr_early{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        maxcorr_late{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,g,1,S.CellInfo.outcomeVal == 2}));
                        if ~isempty(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model)
                            for icell = 1:numel(Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning)
                                map = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;
                                mapbase = Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).meanrespModel;

                                map = map - nanmean(map);
                                map = map./sqrt(sum(map.^2));
                                mapbase = mapbase - nanmean(mapbase);
                                mapbase = mapbase./sqrt(sum(mapbase.^2));
                                fieldXcorr_all = zeros(1,numel(mapbase));
                                xshiftlim = floor(numel(mapbase)/2);
                                outcorrXrange = [1:floor(numel(mapbase)/4) floor(numel(mapbase)*3/4)+1:numel(mapbase)];
                                ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                                for xshift = -xshiftlim:xshiftlim-1
                                    ishift = ishift + 1;
                                    fieldXcorr_all(ishift) = map*circshift(mapbase,xshift)';
                                end
                                fieldXcorr_all(outcorrXrange) = 0;

                                maxcorr_all = getCircularAverage(fieldXcorr_all(:),0,0.1,0.05);
                                maxcorr{g}(icell) = maxcorr_all;
                                stdmaxcorr = 0;
                                kfold = size(Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel,1);
                                %SEM computed using Jacknife method
                                for i = 1:kfold
                                    map_iter = Smaps1d.maps1d{allcontidx,g,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%map;%
                                    mapbase_iter = mapbase;%Smaps1d.maps1d{allcontidx,2,1,S.CellInfo.outcomeVal == 2}.model.tuning(icell).respModel(i,:);%
                                    
                                    map_iter = map_iter - nanmean(map_iter);
                                    map_iter = map_iter./sqrt(sum(map_iter.^2));
                                    mapbase_iter = mapbase_iter - nanmean(mapbase_iter);
                                    mapbase_iter = mapbase_iter./sqrt(sum(mapbase_iter.^2));
                                    fieldXcorr_iter = zeros(1,numel(mapbase_iter));
                                    xshiftlim = floor(numel(mapbase_iter)/2);
                                    outcorrXrange = [1:floor(numel(mapbase)/4) floor(numel(mapbase)*3/4)+1:numel(mapbase)];
                                    ishift = 0;%floor(numel(mapbase_iter)/2)-xshiftlim-1;
                                    for xshift = -xshiftlim:xshiftlim-1
                                        ishift = ishift + 1;
                                        fieldXcorr_iter(ishift) = map_iter*circshift(mapbase_iter,xshift)';
                                    end
                                    fieldXcorr_iter(outcorrXrange) = 0;
                                    %                                 if sum(isnan(fieldXcorr_iter)) == 0
                                    %                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp,'spline');
                                    %                                 else
                                    %                                     fieldXcorr_iter = interp1(xorig,fieldXcorr_iter,xinterp);
                                    %                                 end
                                    maxcorr_iter = getCircularAverage(fieldXcorr_iter(:),0,0.1,0.05);
                                    numbins = numel(fieldXcorr_iter);
                                    %                                 [~ , maxcorr_iter] = max(fieldXcorr_iter);
                                    %                                 maxcorr_iter = xinterp(maxcorr_iter) - xinterp(floor(size(fieldXcorr_iter,2)/2) + 1);
                                    stdmaxcorr = stdmaxcorr + (kfold - 1)/kfold*(numbins/(2*pi)*circ_dist(2*pi/numbins*maxcorr_iter,2*pi/numbins*maxcorr_all)).^2;
                                end
                                stdmaxcorr = (stdmaxcorr).^0.5;
                                maxcorrstd{g}(icell) = stdmaxcorr;

                                %doing the same for early and late theta phase
                                %now
                                numBinsY = size(S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2},2);
                                map_centered = squeeze(S.CellInfo.field2dXPhstheta{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell,:,:));
                                phsoffset = S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,g,1,S.CellInfo.outcomeVal == 2}(icell);
                                if isnan(phsoffset)
                                    phsoffset = 0;
                                end
                                map_centered = circshift(map_centered,numBinsY/2 - floor(phsoffset/360*numBinsY));
                                mapearly = mean(map_centered(earlyphs,:),1);
                                mapbase_centered = squeeze(S.CellInfo.field2dXPhstheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(icell,:,:));
                                phsoffset = S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(icell);
                                if isnan(phsoffset)
                                    phsoffset = 0;
                                end
                                mapbase_centered = circshift(mapbase_centered,numBinsY/2 - floor(phsoffset/360*numBinsY));
                                mapearlybase = mean(mapbase_centered(earlyphs,:),1);

                                mapearly = mapearly - nanmean(mapearly);
                                mapearly = mapearly./sqrt(sum(mapearly.^2));
                                mapearlybase = mapearlybase - nanmean(mapearlybase);
                                mapearlybase = mapearlybase./sqrt(sum(mapearlybase.^2));
                                fieldXcorr_allearly = zeros(1,numel(mapearlybase));
                                xshiftlim = floor(numel(mapearlybase)/2);
                                outcorrXrange = [1:floor(numel(mapearlybase)/4) floor(numel(mapearlybase)*3/4)+1:numel(mapearlybase)];
                                ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                                for xshift = -xshiftlim:xshiftlim-1
                                    ishift = ishift + 1;
                                    fieldXcorr_allearly(ishift) = mapearly*circshift(mapearlybase,xshift)';
                                end
                                fieldXcorr_allearly(outcorrXrange) = 0;
                                maxcorr_allearly = getCircularAverage(fieldXcorr_allearly(:),0,0.1,0.05);
                                maxcorr_early{g}(icell) = maxcorr_allearly;

                                maplate = mean(map_centered(latephs,:),1);
                                maplatebase = mean(mapbase_centered(latephs,:),1);

                                maplate = maplate - nanmean(maplate);
                                maplate = maplate./sqrt(sum(maplate.^2));
                                maplatebase = maplatebase - nanmean(maplatebase);
                                maplatebase = maplatebase./sqrt(sum(maplatebase.^2));
                                fieldXcorr_alllate = zeros(1,numel(maplatebase));
                                xshiftlim = floor(numel(maplatebase)/2);
                                outcorrXrange = [1:floor(numel(maplatebase)/4) floor(numel(maplatebase)*3/4)+1:numel(maplatebase)];
                                ishift = 0;%floor(numel(mapbase)/2)-xshiftlim-1;
                                for xshift = -xshiftlim:xshiftlim-1
                                    ishift = ishift + 1;
                                    fieldXcorr_alllate(ishift) = maplate*circshift(maplatebase,xshift)';
                                end
                                fieldXcorr_alllate(outcorrXrange) = 0;
                                maxcorr_alllate = getCircularAverage(fieldXcorr_alllate(:),0,0.1,0.05);
                                maxcorr_late{g}(icell) = maxcorr_alllate;

                            end
                        end
                        cellprop.fieldXcorrMax{g} = cat(2, cellprop.fieldXcorrMax{g}, maxcorr{g});
                        cellprop.fieldXcorrMaxSE{g} = cat(2, cellprop.fieldXcorrMaxSE{g}, maxcorrstd{g});
                        cellprop.fieldXcorrMax_early{g} = cat(2, cellprop.fieldXcorrMax_early{g}, maxcorr_early{g});
                        cellprop.fieldXcorrMax_late{g} = cat(2, cellprop.fieldXcorrMax_late{g}, maxcorr_late{g});
                        
                        if ~isempty(Slfpspkcoherence.resCA1V1)
                            if isempty(Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll)
                                lfpspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll';
                                lfpPhsspkcoherence = Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll';
                                lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                            else
                                lfpspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_CohSpecChAll))';
                                lfpPhsspkcoherence = (cat(2,Slfpspkcoherence.resCA1V1(1).spk_PhsCohSpecChAll,Slfpspkcoherence.resCA1V1(2).spk_PhsCohSpecChAll))';
                                lfpspkthetaPhscoherence = circ_mean(lfpPhsspkcoherence(:,Slfpspkcoherence.resCA1V1(1).f>=6 & Slfpspkcoherence.resCA1V1(1).f<=9),[],2);
                            end
                            cellprop.lfpSpkCoherence{g} = cat(1,cellprop.lfpSpkCoherence{g},lfpspkcoherence);
                            cellprop.lfpSpkThetaPhsCoherence{g} = cat(1,cellprop.lfpSpkThetaPhsCoherence{g},lfpspkthetaPhscoherence);
                        end
                    else
                        cellprop.field{g} = cat(1,cellprop.field{g}, NaN(size(S.CellInfo.field{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldSE{g} = cat(1,cellprop.fieldSE{g}, NaN(size(S.CellInfo.fieldSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShf{g} = cat(1,cellprop.fieldShf{g}, NaN(size(S.CellInfo.fieldShf{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfSE{g} = cat(1,cellprop.fieldShfSE{g}, NaN(size(S.CellInfo.fieldShfSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldZ{g} = cat(2,cellprop.fieldZ{g}, NaN(size(S.CellInfo.fieldZ{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldXcorr{g} = cat(1,cellprop.fieldXcorr{g}, NaN(size(S.CellInfo.fieldXcorr{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldPos{g} = cat(2,cellprop.fieldPos{g}, NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldPosSE{g} = cat(2,cellprop.fieldPosSE{g}, NaN(size(S.CellInfo.fieldPosSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldCOM{g} = cat(2,cellprop.fieldCOM{g}, NaN(size(S.CellInfo.fieldCOM{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldAmp{g} = cat(2,cellprop.fieldAmp{g},NaN(size(S.CellInfo.fieldAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldAmpSE{g} = cat(2,cellprop.fieldAmpSE{g},NaN(size(S.CellInfo.fieldAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfAmp{g} = cat(2,cellprop.fieldShfAmp{g},NaN(size(S.CellInfo.fieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfAmpSQ{g,1} = cat(2,cellprop.fieldShfAmpSQ{g,1},NaN(size(S.CellInfo.fieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfAmpSQ{g,2} = cat(2,cellprop.fieldShfAmpSQ{g,2},NaN(size(S.CellInfo.fieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfAmpSQ{g,3} = cat(2,cellprop.fieldShfAmpSQ{g,3},NaN(size(S.CellInfo.fieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.rate{g} = cat(2,cellprop.rate{g}, NaN(size(S.CellInfo.rate{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldMax{g} = cat(2,cellprop.fieldMax{g}, NaN(size(S.CellInfo.fieldMax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldMin{g} = cat(2,cellprop.fieldMin{g}, NaN(size(S.CellInfo.fieldMin{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));

                        cellprop.SSI{g} = cat(2,cellprop.SSI{g}, NaN(size(S.CellInfo.SSI{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.SpatialInfo{g} = cat(2,cellprop.SpatialInfo{g}, NaN(size(S.CellInfo.SpatialInfo{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.SpatialInfoPerSpike{g} = cat(2,cellprop.SpatialInfoPerSpike{g}, NaN(size(S.CellInfo.SpatialInfoPerSpike{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.fieldAmpZ{g} = cat(2,cellprop.fieldAmpZ{g},NaN(size(S.CellInfo.fieldAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfAmpZ{g} = cat(2,cellprop.fieldShfAmpZ{g},NaN(size(S.CellInfo.fieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.rategain{g} = cat(2,cellprop.rategain{g},NaN(size(S.CellInfo.rategain{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.rategainSE{g} = cat(2,cellprop.rategainSE{g},NaN(size(S.CellInfo.rategainSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.rategainZ{g} = cat(2,cellprop.rategainZ{g},NaN(size(S.CellInfo.rategain{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.rategainShfdiffZ{g} = cat(2,cellprop.rategainShfdiffZ{g},NaN(size(S.CellInfo.rategain{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.fieldXgainShift{g} = cat(2,cellprop.fieldXgainShift{g},NaN(size(S.CellInfo.fieldShift{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldXgainShiftSE{g} = cat(2,cellprop.fieldXgainShiftSE{g},NaN(size(S.CellInfo.fieldShiftSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldXgainShiftZ{g} = cat(2,cellprop.fieldXgainShiftZ{g},NaN(size(S.CellInfo.fieldShift{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfXgainShiftZ{g} = cat(2,cellprop.fieldShfXgainShiftZ{g},NaN(size(S.CellInfo.fieldShift{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
%                         cellprop.fieldsize{g} = cat(2,cellprop.fieldsize{g},NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldsize{g} = cat(2,cellprop.fieldsize{g},NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldsizeSE{g} = cat(2,cellprop.fieldsize{g},NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldsize_ahead{g} = cat(2,cellprop.fieldsize_ahead{g},NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldsize_behind{g} = cat(2,cellprop.fieldsize_behind{g},NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field_set1{g} = cat(1,cellprop.field_set1{g}, NaN(size(S.CellInfo.field{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field_set2{g} = cat(1,cellprop.field_set2{g}, NaN(size(S.CellInfo.field{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.fieldunwrap{g} = cat(1,cellprop.fieldunwrap{g},NaN(size(S.CellInfo.unwrapfield{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldunwrapSE{g} = cat(1,cellprop.fieldunwrapSE{g},NaN(size(S.CellInfo.unwrapfieldSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldunwrap_set1{g} = cat(1,cellprop.fieldunwrap_set1{g},NaN(size(squeeze(S.CellInfo.unwrapfield_2fold{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(:,1,:)))));
                        cellprop.fieldunwrap_set2{g} = cat(1,cellprop.fieldunwrap_set2{g},NaN(size(squeeze(S.CellInfo.unwrapfield_2fold{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(:,2,:)))));
                        
                        cellprop.fieldunwraphalfdiff{g} = cat(2,cellprop.fieldunwraphalfdiff{g},NaN(size(S.CellInfo.unwrapfieldhaldiff{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldunwraphalfdiffSE{g} = cat(2,cellprop.fieldunwraphalfdiffSE{g},NaN(size(S.CellInfo.unwrapfieldhaldiffSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldunwraphalfdiffZ{g} = cat(2,cellprop.fieldunwraphalfdiffZ{g},NaN(size(S.CellInfo.unwrapfieldhaldiff{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.fieldShfunwraphalfdiffZ{g} = cat(2,cellprop.fieldShfunwraphalfdiffZ{g},NaN(size(S.CellInfo.unwrapfieldhaldiff{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXPhstheta{g} = cat(1,cellprop.field2dXPhstheta{g},NaN(size(S.CellInfo.field2dXPhstheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXPhsthetaSE{g} = cat(1,cellprop.field2dXPhsthetaSE{g},NaN(size(S.CellInfo.field2dXPhsthetaSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXthetapos{g} = cat(1,cellprop.field2dXthetapos{g},NaN(size(S.CellInfo.field2dXthetapos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXthetaposSE{g} = cat(1,cellprop.field2dXthetaposSE{g},NaN(size(S.CellInfo.field2dXthetaposSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
%                         cellprop.field2dXthetaposNorm{g} = cat(1,cellprop.field2dXthetaposNorm{g},NaN(size(S.CellInfo.field2dXthetaposNorm{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
%                         cellprop.field2dXthetaposSENorm{g} = cat(1,cellprop.field2dXthetaposSENorm{g},NaN(size(S.CellInfo.field2dXthetaposSENorm{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXcorrtheta{g} = cat(1,cellprop.field2dXcorrtheta{g},NaN(size(S.CellInfo.field2dXcorrtheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaSE{g} = cat(1,cellprop.field2dXcorrthetaSE{g},NaN(size(S.CellInfo.field2dXcorrthetaSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamax{g} = cat(1,cellprop.field2dXcorrthetamax{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxSE{g} = cat(1,cellprop.field2dXcorrthetamaxSE{g},NaN(size(S.CellInfo.field2dXcorrthetamaxSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxAmp{g} = cat(2,cellprop.field2dXcorrthetamaxAmp{g},NaN(size(S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxAmpSE{g} = cat(2,cellprop.field2dXcorrthetamaxAmpSE{g},NaN(size(S.CellInfo.field2dXcorrthetamaxAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxsinAmp{g} = cat(2,cellprop.field2dXcorrthetamaxsinAmp{g},NaN(size(S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxsinAmpSE{g} = cat(2,cellprop.field2dXcorrthetamaxsinAmpSE{g},NaN(size(S.CellInfo.field2dXcorrthetamaxsinAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxOffset{g} = cat(2,cellprop.field2dXcorrthetamaxOffset{g},NaN(size(S.CellInfo.field2dXcorrthetamaxOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetamaxOffsetSE{g} = cat(2,cellprop.field2dXcorrthetamaxOffsetSE{g},NaN(size(S.CellInfo.field2dXcorrthetamaxOffsetSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXcorrthetaShf{g} = cat(1,cellprop.field2dXcorrthetaShf{g},NaN(size(S.CellInfo.field2dXcorrthetaShf{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfSE{g} = cat(1,cellprop.field2dXcorrthetaShfSE{g},NaN(size(S.CellInfo.field2dXcorrthetaShfSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmax{g} = cat(1,cellprop.field2dXcorrthetaShfmax{g},NaN(size(S.CellInfo.field2dXcorrthetaShfmax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxSE{g} = cat(1,cellprop.field2dXcorrthetaShfmaxSE{g},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxAmp{g} = cat(2,cellprop.field2dXcorrthetaShfmaxAmp{g},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxAmpSQ{g,1} = cat(2,cellprop.field2dXcorrthetaShfmaxAmpSQ{g,1},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxAmpSQ{g,2} = cat(2,cellprop.field2dXcorrthetaShfmaxAmpSQ{g,2},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxAmpSQ{g,3} = cat(2,cellprop.field2dXcorrthetaShfmaxAmpSQ{g,3},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxsinAmp{g} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmp{g},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,1} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,1},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,2} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,2},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,3} = cat(2,cellprop.field2dXcorrthetaShfmaxsinAmpSQ{g,3},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXcorrthetaShfmaxOffset{g} = cat(2,cellprop.field2dXcorrthetaShfmaxOffset{g},NaN(size(S.CellInfo.field2dXcorrthetaShfmaxOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.ZAmpthetapos{g} = cat(2,cellprop.ZAmpthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.ZsinAmpthetapos{g} = cat(2,cellprop.ZsinAmpthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.ZAmpthetaShfpos{g} = cat(2,cellprop.ZAmpthetaShfpos{g},NaN(size(S.CellInfo.field2dXcorrthetamaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.ZsinAmpthetaShfpos{g} = cat(2,cellprop.ZsinAmpthetaShfpos{g},NaN(size(S.CellInfo.field2dXcorrthetamaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXRefcorrtheta{g} = cat(1,cellprop.field2dXRefcorrtheta{g},NaN(size(S.CellInfo.field2dXRefcorrtheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaSE{g} = cat(1,cellprop.field2dXRefcorrthetaSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetaSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamax{g} = cat(1,cellprop.field2dXRefcorrthetamax{g},NaN(size(S.CellInfo.field2dXRefcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxSE{g} = cat(1,cellprop.field2dXRefcorrthetamaxSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxAmp{g} = cat(2,cellprop.field2dXRefcorrthetamaxAmp{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxAmpSE{g} = cat(2,cellprop.field2dXRefcorrthetamaxAmpSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxsinAmp{g} = cat(2,cellprop.field2dXRefcorrthetamaxsinAmp{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxsinAmpSE{g} = cat(2,cellprop.field2dXRefcorrthetamaxsinAmpSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxsinAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxOffset{g} = cat(2,cellprop.field2dXRefcorrthetamaxOffset{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetamaxOffsetSE{g} = cat(2,cellprop.field2dXRefcorrthetamaxOffsetSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxOffsetSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXRefcorrthetaShf{g} = cat(1,cellprop.field2dXRefcorrthetaShf{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShf{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfSE{g} = cat(1,cellprop.field2dXRefcorrthetaShfSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmax{g} = cat(1,cellprop.field2dXRefcorrthetaShfmax{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxSE{g} = cat(1,cellprop.field2dXRefcorrthetaShfmaxSE{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxAmp{g} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmp{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,1} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,1},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,2} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,2},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,3} = cat(2,cellprop.field2dXRefcorrthetaShfmaxAmpSQ{g,3},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmp{g} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmp{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,1} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,1},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,2} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,2},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,3} = cat(2,cellprop.field2dXRefcorrthetaShfmaxsinAmpSQ{g,3},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dXRefcorrthetaShfmaxOffset{g} = cat(2,cellprop.field2dXRefcorrthetaShfmaxOffset{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.ZAmpRefthetapos{g} = cat(2,cellprop.ZAmpRefthetapos{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.ZsinAmpRefthetapos{g} = cat(2,cellprop.ZsinAmpRefthetapos{g},NaN(size(S.CellInfo.field2dXRefcorrthetamaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.ZAmpRefthetaShfpos{g} = cat(2,cellprop.ZAmpRefthetaShfpos{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.ZsinAmpRefthetaShfpos{g} = cat(2,cellprop.ZsinAmpRefthetaShfpos{g},NaN(size(S.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dslopeXY{g} = cat(2,cellprop.field2dslopeXY{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dphi0XY{g} = cat(2,cellprop.field2dphi0XY{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2drhoXY{g} = cat(2,cellprop.field2drhoXY{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dslopeXYSE{g} = cat(2,cellprop.field2dslopeXYSE{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dphi0XYSE{g} = cat(2,cellprop.field2dphi0XYSE{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2drhoXYSE{g} = cat(2,cellprop.field2drhoXYSE{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dShfslopeXY{g} = cat(2,cellprop.field2dShfslopeXY{g},NaN(size(S.CellInfo.field2dShfslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfphi0XY{g} = cat(2,cellprop.field2dShfphi0XY{g},NaN(size(S.CellInfo.field2dShfslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfrhoXY{g} = cat(2,cellprop.field2dShfrhoXY{g},NaN(size(S.CellInfo.field2dShfslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfslopeXYSQ{g,1} = cat(2,cellprop.field2dShfslopeXYSQ{g,1},NaN(size(S.CellInfo.field2dShfslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfslopeXYSQ{g,2} = cat(2,cellprop.field2dShfslopeXYSQ{g,2},NaN(size(S.CellInfo.field2dShfslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfslopeXYSQ{g,3} = cat(2,cellprop.field2dShfslopeXYSQ{g,3},NaN(size(S.CellInfo.field2dShfslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfphi0XYSQ{g,1} = cat(2,cellprop.field2dShfphi0XYSQ{g,1},NaN(size(S.CellInfo.field2dShfphi0XY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfphi0XYSQ{g,2} = cat(2,cellprop.field2dShfphi0XYSQ{g,2},NaN(size(S.CellInfo.field2dShfphi0XY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfphi0XYSQ{g,3} = cat(2,cellprop.field2dShfphi0XYSQ{g,3},NaN(size(S.CellInfo.field2dShfphi0XY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfrhoXYSQ{g,1} = cat(2,cellprop.field2dShfrhoXYSQ{g,1},NaN(size(S.CellInfo.field2dShfrhoXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfrhoXYSQ{g,2} = cat(2,cellprop.field2dShfrhoXYSQ{g,2},NaN(size(S.CellInfo.field2dShfrhoXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfrhoXYSQ{g,3} = cat(2,cellprop.field2dShfrhoXYSQ{g,3},NaN(size(S.CellInfo.field2dShfrhoXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dslopeXYZ{g} = cat(2,cellprop.field2dslopeXYZ{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2drhoXYZ{g} = cat(2,cellprop.field2drhoXYZ{g},NaN(size(S.CellInfo.field2drhoXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfslopeXYZ{g} = cat(2,cellprop.field2dShfslopeXYZ{g},NaN(size(S.CellInfo.field2dslopeXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dShfrhoXYZ{g} = cat(2,cellprop.field2dShfrhoXYZ{g},NaN(size(S.CellInfo.field2drhoXY{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
%                         cellprop.field2dXcorrtheta_toref{g} = cat(1,cellprop.field2dXcorrtheta_toref{g},NaN(size(S.CellInfo.field2dXcorrtheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
%                         cellprop.field2dXcorrthetamax_toref{g} = cat(1,cellprop.field2dXcorrthetamax_toref{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field2dXcorrthetamax_set1{g} = cat(1,cellprop.field2dXcorrthetamax_set1{g},NaN(size(squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(:,1,:)))));
                        cellprop.field2dXcorrthetamax_set2{g} = cat(1,cellprop.field2dXcorrthetamax_set2{g},NaN(size(squeeze(S.CellInfo.field2dXcorrthetamax_2fold{allcontidx,2,1,S.CellInfo.outcomeVal == 2}(:,2,:)))));
                        
                        cellprop.phsfieldPos{g} = cat(2,cellprop.phsfieldPos{g},NaN(size(S.CellInfo.phsfieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldPosSE{g} = cat(2,cellprop.phsfieldPosSE{g},NaN(size(S.CellInfo.phsfieldPosSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldCOM{g} = cat(2,cellprop.phsfieldCOM{g},NaN(size(S.CellInfo.phsfieldCOM{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldCOMSE{g} = cat(2,cellprop.phsfieldCOMSE{g},NaN(size(S.CellInfo.phsfieldCOMSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfield{g} = cat(1,cellprop.phsfield{g},NaN(size(S.CellInfo.phsfield{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldSE{g} = cat(1,cellprop.phsfieldSE{g},NaN(size(S.CellInfo.phsfieldSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldAmp{g} = cat(2,cellprop.phsfieldAmp{g},NaN(size(S.CellInfo.phsfieldAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldAmpSE{g} = cat(2,cellprop.phsfieldAmpSE{g},NaN(size(S.CellInfo.phsfieldAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldsinAmp{g} = cat(2,cellprop.phsfieldsinAmp{g},NaN(size(S.CellInfo.phsfieldsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldsinAmpSE{g} = cat(2,cellprop.phsfieldsinAmpSE{g},NaN(size(S.CellInfo.phsfieldsinAmpSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldsinOffset{g} = cat(2,cellprop.phsfieldsinOffset{g},NaN(size(S.CellInfo.phsfieldsinOffset{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldsinOffsetSE{g} = cat(2,cellprop.phsfieldsinOffsetSE{g},NaN(size(S.CellInfo.phsfieldsinOffsetSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShf{g} = cat(1,cellprop.phsfieldShf{g},NaN(size(S.CellInfo.phsfieldShf{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfSE{g} = cat(1,cellprop.phsfieldShfSE{g},NaN(size(S.CellInfo.phsfieldShfSE{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfAmp{g} = cat(2,cellprop.phsfieldShfAmp{g},NaN(size(S.CellInfo.phsfieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfAmpSQ{g,1} = cat(2,cellprop.phsfieldShfAmpSQ{g,1},NaN(size(S.CellInfo.phsfieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfAmpSQ{g,2} = cat(2,cellprop.phsfieldShfAmpSQ{g,2},NaN(size(S.CellInfo.phsfieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfAmpSQ{g,3} = cat(2,cellprop.phsfieldShfAmpSQ{g,3},NaN(size(S.CellInfo.phsfieldShfAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfsinAmp{g} = cat(2,cellprop.phsfieldShfsinAmp{g},NaN(size(S.CellInfo.phsfieldShfsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfsinAmpSQ{g,1} = cat(2,cellprop.phsfieldShfsinAmpSQ{g,1},NaN(size(S.CellInfo.phsfieldShfsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfsinAmpSQ{g,2} = cat(2,cellprop.phsfieldShfsinAmpSQ{g,2},NaN(size(S.CellInfo.phsfieldShfsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfsinAmpSQ{g,3} = cat(2,cellprop.phsfieldShfsinAmpSQ{g,3},NaN(size(S.CellInfo.phsfieldShfsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
%                         cellprop.phsfieldZ{g} = cat(2,cellprop.phsfieldZ{g},NaN(size(S.CellInfo.phsfieldZ{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsmodulation{g} = cat(2,cellprop.phsmodulation{g},NaN(size(S.CellInfo.phsmodulation{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.phsfieldAmpZ{g} = cat(2,cellprop.phsfieldAmpZ{g},NaN(size(S.CellInfo.phsfieldAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldsinAmpZ{g} = cat(2,cellprop.phsfieldsinAmpZ{g},NaN(size(S.CellInfo.phsfieldsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfAmpZ{g} = cat(2,cellprop.phsfieldShfAmpZ{g},NaN(size(S.CellInfo.phsfieldAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.phsfieldShfsinAmpZ{g} = cat(2,cellprop.phsfieldShfsinAmpZ{g},NaN(size(S.CellInfo.phsfieldsinAmp{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.ZXmaxthetapos{g} = cat(1,cellprop.ZXmaxthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.ZXminthetapos{g} = cat(1,cellprop.ZXminthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.Xmaxthetapos{g} = cat(1,cellprop.Xmaxthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.Xminthetapos{g} = cat(1,cellprop.Xminthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.ZPhsmaxtheta{g} = cat(1,cellprop.ZPhsmaxtheta{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.ZPhsmintheta{g} = cat(1,cellprop.ZPhsmintheta{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsthetaMod{g} = cat(1,cellprop.PhsthetaMod{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXmaxthetapos{g} = cat(1,cellprop.PhsXmaxthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXminthetapos{g} = cat(1,cellprop.PhsXminthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXaheadthetapos{g} = cat(1,cellprop.PhsXaheadthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        cellprop.PhsXbehindthetapos{g} = cat(1,cellprop.PhsXbehindthetapos{g},NaN(size(S.CellInfo.field2dXcorrthetamax{allcontidx,2,1,S.CellInfo.outcomeVal == 2},1),1));
                        
                        cellprop.PhaseRho{g} = cat(2,cellprop.PhaseRho{g},NaN(size(S.CellInfo.PhaseRho{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.PhasePval{g} = cat(2,cellprop.PhasePval{g},NaN(size(S.CellInfo.PhasePval{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field2dPhstheta{g} = cat(1,cellprop.field2dPhstheta{g},NaN(size(S.CellInfo.field2dPhstheta{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.PhaseRayleighPval{g} = cat(2,cellprop.PhaseRayleighPval{g},NaN(size(S.CellInfo.PhaseRayleighPval{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.PhaseRayleighZ{g} = cat(2,cellprop.PhaseRayleighZ{g},NaN(size(S.CellInfo.PhaseRayleighZ{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        cellprop.field_set1{g} = cat(1,cellprop.field_set1{g}, NaN(size(S.CellInfo.field_set1{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        cellprop.field_set2{g} = cat(1,cellprop.field_set2{g}, NaN(size(S.CellInfo.field_set1{allcontidx,2,1,S.CellInfo.outcomeVal == 2})));
                        
                        reliabilityCorr = NaN(numel(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2}),1);
                        cellprop.reliabilityCorr{g} = cat(1,cellprop.reliabilityCorr{g},reliabilityCorr);
                        cellprop.thetareliabilityCorr{g} = cat(1,cellprop.thetareliabilityCorr{g},reliabilityCorr);
                        
                        maxcorr{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2}));
                        maxcorrstd{g} = NaN(size(S.CellInfo.fieldPos{allcontidx,2,1,S.CellInfo.outcomeVal == 2}));
                        cellprop.fieldXcorrMax{g} = cat(2, cellprop.fieldXcorrMax{g}, maxcorr{g});
                        cellprop.fieldXcorrMaxSE{g} = cat(2, cellprop.fieldXcorrMaxSE{g}, maxcorrstd{g});
                        cellprop.fieldXcorrMax_early{g} = cat(2, cellprop.fieldXcorrMax_early{g}, maxcorr{g});
                        cellprop.fieldXcorrMax_late{g} = cat(2, cellprop.fieldXcorrMax_late{g}, maxcorr{g});
                        
                        if ~isempty(Slfpspkcoherence.resCA1V1)
                            lfpspkcoherence = NaN(size(Slfpspkcoherence.resCA1V1(1).spk_CohSpecChAll'));
                            lfpspkthetaPhscoherence = NaN(size(lfpspkcoherence,1),1);
                            cellprop.lfpSpkCoherence{g} = cat(1,cellprop.lfpSpkCoherence{g},lfpspkcoherence);
                            cellprop.lfpSpkThetaPhsCoherence{g} = cat(1,cellprop.lfpSpkThetaPhsCoherence{g},lfpspkthetaPhscoherence);
                        end
                    end
                end
            end
        end
    end
end
% save(savedfilename_cellprop, 'cellprop','-v7.3');
end