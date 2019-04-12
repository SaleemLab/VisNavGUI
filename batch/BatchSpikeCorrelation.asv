function popresSpkCorr = BatchSpikeCorrelation(batch2p,ErrType)
SetDirs;
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
    datadir = DIRS.data2p;
else
    expt = getExperimentList;
    datadir = 'D:\DATA\batch';%DIRS.multichanspikes;
end

if nargin < 2
    ErrType = 'Error';%'Blank';
end
popresSpkCorr.ErrType = ErrType;
popresSpkCorr.Tsmthwin = 15;%250;%250;%150;%300;%40;%120;%50
popresSpkCorr.Tsmthwin_dec = 250;
popresSpkCorr.nDecbins = 100;
popresSpkCorr.Tcorr = 250;
popresSpkCorr.nbinscorr = 20;%100;%
popresSpkCorr.Xsmthwin = 4;%2%1;%
popresSpkCorr.SpdSmthWin = popresSpkCorr.Tcorr;%popresSpkCorr.Tsmthwin;
popresSpkCorr.SpeedThreshold = 5;
popresSpkCorr.nspeedbins = 5;
popresSpkCorr.neyebins = 1;
popresSpkCorr.nthetaphsbins = 1;%1;%
popresSpkCorr.nphsbins = 1;
popresSpkCorr.cellstr = 'goodonly';%'All_50bins';%'goodonly';%'goodonly_unwrapped';%'goodonly';%'All';%
filesuffix_EXP = ['Twin' num2str(popresSpkCorr.Tsmthwin) '_' 'Xwin' num2str(popresSpkCorr.Xsmthwin) '_' 'spdth' num2str(popresSpkCorr.SpeedThreshold) '_' 'Decwin' num2str(popresSpkCorr.Tsmthwin_dec) '_' 'nDecbins' num2str(popresSpkCorr.nDecbins) '_' num2str(popresSpkCorr.nspeedbins) 'speedbins' '_' num2str(popresSpkCorr.neyebins) 'eyebins' '_' num2str(popresSpkCorr.nphsbins) 'thetabins' '_' popresSpkCorr.cellstr];
disp(filesuffix_EXP);

popresSpkCorr.Pth_field = 0.01;
popresSpkCorr.sampleRate = 60;
popresSpkCorr.nSpdbins = 5;%1;%
popresSpkCorr.nEyebins = 1;%3;%3;
% popresSpkCorr.nXbins = 20;%100;
if strcmp(popresSpkCorr.ErrType,'Blank')
    disp('will look only at blank periods')
    popresSpkCorr.nXbins = 1;%
else
    popresSpkCorr.nXbins = 1;%100;
end

popresSpkCorr.Nshuffle = 100;
popresSpkCorr.PthCrossCA1V1 = 0.05;
lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

contval = [0.1:0.05:0.9];%[0.2 0.3 0.4];%[0.8 0.9];%
outvalcorr = 2;%[0 1 2 3 4 5];%[0 1 2 3 4];%

savedfilename_popresCorr = ['D:\DATA\batch\All\Correlations\popresSpkCorr_' popresSpkCorr.ErrType '_Twin' num2str(popresSpkCorr.Tsmthwin) '_Xwin' num2str(popresSpkCorr.Xsmthwin) '_spdth' num2str(popresSpkCorr.SpeedThreshold)...
                     '_Decwin' num2str(popresSpkCorr.Tsmthwin_dec) '_nDecbins' num2str(popresSpkCorr.nDecbins) '_' num2str(popresSpkCorr.nspeedbins) 'speedbins_' num2str(popresSpkCorr.neyebins) 'eyebins_' num2str(popresSpkCorr.nphsbins) 'thetabins_' popresSpkCorr.cellstr '.mat'];


for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
        if exist([dDIRname filesep 'EXP_' filesuffix_EXP '.mat'],'file')
            if (expt(ianimal).goodCA1{iseries} == 1) && (expt(ianimal).goodV1{iseries} == 1)
                S = load([dDIRname filesep 'EXP_' filesuffix_EXP '.mat']);
                EXP = TVRData;
                EXP.Copyobj(S.EXP);
                S = [];
                
                cont_list = find(ismember(EXP.SubsetVal.contrast, contval));
                RL_list = find(ismember(EXP.SubsetVal.roomlength, [1]));
                outcome_list = find(ismember(EXP.SubsetVal.outcome, outvalcorr));
                Xrange = max(floor(EXP.Bayes.X));%1;%
                speeds = NaN(size(EXP.data.es.ballspeed));
                speeds(~isnan(EXP.data.es.ballspeed)) = smthInTime(EXP.data.es.ballspeed(~isnan(EXP.data.es.ballspeed)), popresSpkCorr.sampleRate, popresSpkCorr.SpdSmthWin, 'same', [], 'boxcar_centered');
                eyeX = NaN(size(EXP.data.es.eyeXpos));
                eyeX(~isnan(EXP.data.es.eyeXpos)) = smthInTime(EXP.data.es.eyeXpos(~isnan(EXP.data.es.eyeXpos)), popresSpkCorr.sampleRate, popresSpkCorr.SpdSmthWin, 'same', [], 'boxcar_centered');
                
                Speedbinned = cell(1,4);
                Eyebinned = cell(1,4);
                
                gooddecCA1 = false(1,3);
                gooddecV1 = false(1,3);
                for g = [2 1 3]
                    gooddecCA1(g) = expt(ianimal).goodCA1dec{g}{iseries};
                    gooddecV1(g) = expt(ianimal).goodV1dec{g}{iseries};
                end
                gooddecCA1(4) = gooddecCA1(2);
                gooddecV1(4) = gooddecV1(2);
                
                cbase = numel(EXP.SubsetVal.contrast) + 1;%find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
                gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
                rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
                obase = find(EXP.SubsetVal.outcome == 2);
                Amp = abs(EXP.CellInfo.fieldAmp{cbase, gbase, rbase, obase});
                ShfAmp = abs(EXP.CellInfo.fieldShfAmpiter{cbase, gbase, rbase, obase})';
                nShf = size(ShfAmp,1);
                ShfAmpZ = sum(repmat(Amp,[nShf 1])< ShfAmp,1)/nShf;
                signiAmpfields = ShfAmpZ <= popresSpkCorr.Pth_field;
                
                min2maxSpkwf = NaN(size(EXP.CellInfo.min2maxSpkwf));
                for icell = 1:size(EXP.CellInfo.waveform,1)
                    [pks,minidx] = findpeaks(-EXP.CellInfo.waveform(icell,:));%find(S.CellInfo.waveform(icell,:) == min(S.CellInfo.waveform(icell,:)));
                    pks = pks(minidx >= 10 & minidx <= size(EXP.CellInfo.waveform,2)-10);
                    minidx = minidx(minidx >= 10 & minidx <= size(EXP.CellInfo.waveform,2)-10);
                    minidx = minidx(pks == max(pks));
                    if minidx >= 10 && minidx <= size(EXP.CellInfo.waveform,2)-10
                        minidx = minidx(1);
                        [~,maxidx] = max(EXP.CellInfo.waveform(icell,(minidx+1):end));
                        min2maxSpkwf(icell) = maxidx;
                    end
                end
                
                goodcells = EXP.CellInfo.Goodcluster & ~EXP.CellInfo.Finterneuron & ~isnan(EXP.CellInfo.fieldPos{cbase,gbase,rbase,obase}) & signiAmpfields;
                ProbeID = EXP.CellInfo.Probe(goodcells);
                min2maxSpkwf = min2maxSpkwf(goodcells);
                fieldPosAll = EXP.CellInfo.fieldPos{cbase,gbase,rbase,obase}(goodcells);
                field = EXP.CellInfo.field{cbase,gbase,rbase,obase}(goodcells,:);
                
                popresSpkCorr.V1min2maxSpkwf{ianimal,iseries} = min2maxSpkwf(ProbeID == 2);
                popresSpkCorr.CA1min2maxSpkwf{ianimal,iseries} = min2maxSpkwf(ProbeID == 1);
                
                ncells = sum(goodcells);
                PosError0{1} = EXP.data.es.spikeTrain(:,goodcells);%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                sampleRate = mean(1./EXP.data.es.sampleSize);
                for icell = 1:ncells
                    PosError0{1}(:,icell) = smthInTime(PosError0{1}(:,icell), sampleRate, popresSpkCorr.Tcorr, 'same', [], 'boxcarsum_centered');%popresSpkCorr.Tsmthwin, 'same', [], 'boxcarsum_centered');
                end
                PosError0{2} = PosError0{1};
                
%                 tvec = 1:size(PosError0{2},1);
%                 for xx = 1:size(PosError0{2},2)
%                     PosError0{2}(:,xx) = interp1(tvec(~isnan(PosError0{2}(:,xx)) & sum(PosError0{2},2)>0),PosError0{2}(~isnan(PosError0{2}(:,xx)) & sum(PosError0{2},2)>0,xx),1:size(PosError0{2}(:,xx),1));
%                 end
%                 for xx = 1:size(PosError0{1},2)
%                     PosError0{1}(:,xx) = interp1(tvec(~isnan(PosError0{1}(:,xx)) & sum(PosError0{1},2)>0),PosError0{1}(~isnan(PosError0{1}(:,xx)) & sum(PosError0{1},2)>0,xx),1:size(PosError0{1}(:,xx),1));
%                 end
                
                for g = [2 1 3]%2%
                    if gooddecCA1(g) && gooddecV1(g)
                        tidx = false(size(EXP.Bayes.X));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    if ~isempty(EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                        tidx = tidx | EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                    end
                                end
                            end
                        end
                        if strcmp(popresSpkCorr.ErrType,'Blank')
                            blanks = ~EXP.Subset.noBlanks;
                            %blanks = blanks & EXP.Bayes.X == mode(EXP.Bayes.X(blanks));
                            for tshift = -30:30
                                blanks = blanks & circshift(~EXP.Subset.noBlanks,tshift);
                            end
                            for tshift = -15:15
                                blanks = blanks & circshift(~EXP.data.es.lick,tshift);
                            end
                            if g==2
                                idxblank_onset = find(diff(blanks)>0);
                                blank_duration = zeros(1,numel(idxblank_onset));
                                for k = 1:numel(idxblank_onset)
                                    blank_duration(k) = find(blanks(idxblank_onset(k)+1:end)<0.5,1,'first');
                                end
                                popresSpkCorr.blank_duration{ianimal,iseries} = blank_duration/60+1;
                            end
                            tidx = blanks & EXP.data.es.smthBallSpd > EXP.Bayes.speed_th & ~isnan(EXP.data.es.smthBallSpd) & ~isnan(EXP.data.es.smthTrajSpd);
                        end                        
                        
                        idxref = tidx;
                        
                        if ~strcmp(popresSpkCorr.ErrType,'Blank')
                            Xbinned = round(EXP.Bayes.X);%floor(((EXP.Bayes.X-1)/Xrange)*popresSpkCorr.nXbins)+1;%EXP.Bayes.X;%
                        else
                            Xbinned = ones(size(EXP.Bayes.X));
                        end
                        Xrange = round(max(Xbinned));
                        
                        itraj = Xbinned;
                        ntrajbins = max(itraj);
                        spdquantilelim = zeros(ntrajbins,2);
                        Speedbinned{g} = NaN(size(EXP.data.es.ballspeed));
                        if popresSpkCorr.nSpdbins > 1
                            for spd = 1:popresSpkCorr.nSpdbins
                                for xx = 1:ntrajbins
                                    spdquantilelim(xx,1) = quantile(speeds(idxref & itraj == xx),max(0,(spd-1)/popresSpkCorr.nSpdbins));
                                    spdquantilelim(xx,2) = quantile(speeds(idxref & itraj == xx),min(1,(spd)/popresSpkCorr.nSpdbins));
                                end
                                Speedbinned{g}(speeds >= spdquantilelim(itraj,1) & speeds < spdquantilelim(itraj,2)) = spd;
                                if spd == 1
                                    Speedbinned{g}(speeds <= spdquantilelim(itraj,1)) = spd;
                                end
                                if spd == popresSpkCorr.nSpdbins
                                    Speedbinned{g}(speeds >= spdquantilelim(itraj,2)) = spd;
                                end
                            end
                        else
                            Speedbinned{g} = ones(size(EXP.data.es.ballspeed));
                        end
                        
%                         if popresSpkCorr.nSpdbins >= 3
%                             Speedbinned{g} = NaN(size(EXP.data.es.ballspeed));
%                             speedprofile = NaN(EXP.Bayes.numBins,1);
%                             for xx = 1:EXP.Bayes.numBins
%                                 speedprofile(xx) = median(speeds(idxref & EXP.Bayes.X==xx));
%                             end
%                             %speedprofile = fast1Dmap(res.X{ianimal,iseries}(idxref),res.runSpeed{ianimal,iseries}(idxref),1,1,1/(EXP.Bayes.Xsmth_win/EXP.Bayes.numBins),EXP.data.es.CircularMaze);
%                             speeds = (speeds - speedprofile(EXP.Bayes.X));%./speedprofile(res.X{ianimal,iseries});
%                             speedrange = 10;%from -5 cm/s to +5 cm/s around the median speed
%                             if popresSpkCorr.nSpdbins > 3
%                                 binsize = speedrange/(popresSpkCorr.nSpdbins-2 - 1);
%                             else
%                                 binsize = 0;
%                             end
%                             minbinval = -floor(speedrange/2) - binsize/2;
%                             maxbinval = floor(speedrange/2) + binsize/2;
%                             speeds(speeds < minbinval) = minbinval-1;
%                             speeds(speeds > maxbinval) = maxbinval+1;
%                             [Speedbinned{g}(speeds >= minbinval & speeds <= maxbinval), ~] = normalise1var(speeds(speeds >= minbinval & speeds <= maxbinval), popresSpkCorr.nSpdbins-2,[],[minbinval maxbinval]);
%                             Speedbinned{g}(speeds >= minbinval & speeds <= maxbinval) = Speedbinned{g}(speeds >= minbinval & speeds <= maxbinval) + 1;
%                             Speedbinned{g}(speeds < minbinval) = 1;
%                             Speedbinned{g}(speeds > maxbinval) = popresSpkCorr.nSpdbins;
%                         else
%                             Speedbinned{g} = ones(size(EXP.data.es.ballspeed));
%                             popresSpkCorr.nSpdbins = 1;
%                         end
                        
                        itraj = EXP.Bayes.X;
                        ntrajbins = max(itraj);
                        eyequantilelim = zeros(ntrajbins,2);
                        Eyebinned{g} = NaN(size(EXP.data.es.eyeXpos));
                        if popresSpkCorr.nEyebins > 1
                            for ieye = 1:popresSpkCorr.nEyebins
                                for xx = 1:ntrajbins
                                    eyequantilelim(xx,1) = quantile(eyeX(idxref & itraj == xx),max(0,(ieye-1)/popresSpkCorr.nEyebins));
                                    eyequantilelim(xx,2) = quantile(eyeX(idxref & itraj == xx),min(1,(ieye)/popresSpkCorr.nEyebins));
                                end
                                Eyebinned{g}(eyeX >= eyequantilelim(itraj,1) & eyeX < eyequantilelim(itraj,2)) = ieye;
                                if ieye == 1
                                    Eyebinned{g}(eyeX <= eyequantilelim(itraj,1)) = ieye;
                                end
                                if ieye == popresSpkCorr.nEyebins
                                    Eyebinned{g}(eyeX >= eyequantilelim(itraj,2)) = ieye;
                                end
                            end
                        else
                            Eyebinned{g} = ones(size(EXP.data.es.eyeXpos));
                        end
                        
                        if ismember(g,[2 1 3])
                            for ieye = 1:popresSpkCorr.nEyebins
                                for spd = 1:popresSpkCorr.nSpdbins
                                    for xx = 1:Xrange
                                        tidx_corr = tidx & Xbinned == xx & Speedbinned{g} == spd & Eyebinned{g} == ieye & ~isnan(sum(PosError0{1},2))  & ~isnan(sum(PosError0{2},2));
                                        idx_corrCA1 = find(tidx_corr);
                                        idx_corrV1 = find(tidx_corr);
                                        PosError0{1}(idx_corrCA1,:) = PosError0{1}(idx_corrCA1,:) - repmat(nanmean(PosError0{1}(idx_corrCA1,:),1),[numel(idx_corrCA1) 1]);
                                        PosError0{2}(idx_corrV1,:) = PosError0{2}(idx_corrV1,:) - repmat(nanmean(PosError0{2}(idx_corrV1,:),1),[numel(idx_corrV1) 1]);
                                    end
                                end
                            end
                        end
%                         PosError0{1}(tidx,:) = PosError0{1}(tidx,:)./sqrt(repmat(sum(PosError0{1}(tidx,:).^2,1),[sum(tidx) 1]));
%                         PosError0{2}(tidx,:) = PosError0{2}(tidx,:)./sqrt(repmat(sum(PosError0{2}(tidx,:).^2,1),[sum(tidx) 1]));
                    end
                end
                    
                for g = [2 1 3]%2%
                    if gooddecCA1(g) && gooddecV1(g)
                        tidx = false(size(EXP.Bayes.X));
                        for cont = 1:numel(cont_list)
                            for r = 1:numel(RL_list)
                                for o = 1:numel(outcome_list)
                                    if g < 4
                                        if ~isempty(EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                            tidx = tidx | EXP.getSubsets(cont_list(cont),g,RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                        end
                                    else
                                        if ~isempty(EXP.getSubsets(cont_list(cont),[1 2 3],RL_list(r),outcome_list(o), EXP.Bayes.speed_th))
                                            tidx = tidx | EXP.getSubsets(cont_list(cont),[1 2 3],RL_list(r),outcome_list(o), EXP.Bayes.speed_th);
                                        end
                                    end
                                end
                            end
                        end
                        if strcmp(popresSpkCorr.ErrType,'Blank')
                            blanks = ~EXP.Subset.noBlanks;
                            for tshift = -30:30
                                blanks = blanks & circshift(~EXP.Subset.noBlanks,tshift);
                            end
                            for tshift = -15:15
                                blanks = blanks & circshift(~EXP.data.es.lick,tshift);
                            end
                            tidx = blanks & EXP.data.es.smthBallSpd > EXP.Bayes.speed_th & ~isnan(EXP.data.es.smthBallSpd) & ~isnan(EXP.data.es.smthTrajSpd);
%                             tidx = tidx & EXP.Bayes.X == mode(EXP.Bayes.X(tidx));
                        end
                        
                        Prange = size(PosError0{1},2);
                        
                        nXbins = popresSpkCorr.nXbins;
                        nSpdbins = 1;%popresSpkCorr.nSpdbins;
                        nEyebins = 1;%popresSpkCorr.nEyebins;
                        s_CA1V1Cross = zeros(Prange,Prange);
                        s_CA1V1Cross_shuffled = zeros(Prange,Prange,popresSpkCorr.Nshuffle);
                        s_CA1marginal = zeros(1,Prange);
                        s_V1marginal = zeros(1,Prange);
                        pnXbins = popresSpkCorr.nXbins;
                        pNshuffle = popresSpkCorr.Nshuffle;
                        for ieye = 1:nEyebins
                            for spd = 1:nSpdbins
                                for xx = 1:nXbins
                                    tidx_corr = tidx;% & (Xbinned >= (xx-1)*100/pnXbins) & (Xbinned <= xx*100/pnXbins) & Speedbinned{g} == spd & Eyebinned{g} == ieye;
                                    idx_corrCA1 = find(tidx_corr);
                                    idx_corrV1 = find(tidx_corr);
                                    s_CA1V1Cross = s_CA1V1Cross + PosError0{2}(idx_corrV1,:)'*PosError0{1}(idx_corrCA1,:);
                                    for ishf = 1:pNshuffle
                                        shuffledidxV1 = circshift(idx_corrV1,round(numel(idx_corrV1)/2) + randi(round(numel(idx_corrV1)/2)) - round(numel(idx_corrV1)/4));
                                        s_CA1V1Cross_shuffled(:,:,ishf) = s_CA1V1Cross_shuffled(:,:,ishf) + PosError0{2}(shuffledidxV1,:)'*PosError0{1}(idx_corrCA1,:);
%                                         shuffledidxV1 = randperm(numel(idx_corrV1));
%                                         s_CA1V1Cross_shuffled(:,:,ishf) = s_CA1V1Cross_shuffled(:,:,ishf) + PosError0{2}(idx_corrV1(shuffledidxV1),:)'*PosError0{1}(idx_corrCA1,:);
                                    end
                                    s_CA1marginal = s_CA1marginal + sum(PosError0{1}(idx_corrCA1,:).^2,1);
                                    s_V1marginal = s_V1marginal + sum(PosError0{2}(idx_corrV1,:).^2,1);
                                end
                            end
                        end
%                         s_CA1V1Cross = squeeze(nansum(nansum(nansum(s_CA1V1Cross,1),2),3));
%                         s_CA1V1Cross_shuffled = squeeze(nansum(nansum(nansum(s_CA1V1Cross_shuffled,1),2),3));
%                         s_CA1marginal = squeeze(nansum(nansum(nansum(s_CA1marginal,1),2),3));
%                         s_V1marginal = squeeze(nansum(nansum(nansum(s_V1marginal,1),2),3));
                        
                        popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g} = s_CA1V1Cross./sqrt(s_V1marginal'*s_CA1marginal);
                        popresSpkCorr.s_CA1V1Cross_shuffled{ianimal,iseries,g} = s_CA1V1Cross_shuffled./repmat(sqrt(s_V1marginal'*s_CA1marginal),[1 1 popresSpkCorr.Nshuffle]);
                        
                        CA1cells = ProbeID == 1;
                        V1cells = ProbeID == 2;
                        [fposCA1,isortCA1] = sort(fieldPosAll(CA1cells));
                        [fposV1,isortV1] = sort(fieldPosAll(V1cells));
                        
                        popresSpkCorr.s_fposCA1{ianimal,iseries,g} = fposCA1;
                        popresSpkCorr.s_fposV1{ianimal,iseries,g} = fposV1;
                        fieldnorm = field;
                        fieldnorm = fieldnorm-repmat(nanmean(fieldnorm,2),[1 size(fieldnorm,2)]);
                        fieldnorm = fieldnorm./sqrt(repmat(nanmean(fieldnorm.^2,2),[1 size(fieldnorm,2)]));
                        fieldnorm_centered = fieldnorm;
                        for icell = 1:size(fieldnorm,1)
                            if sum(isnan(fieldnorm_centered(icell,:))) == 0
                                fieldnorm_centered(icell,:) = circshift(fieldnorm_centered(icell,:),round(floor(size(fieldnorm,2)/2) - fieldPosAll(icell)));
                            end
                        end
                        fielddist = NaN(size(fieldnorm,1));
                        try
                        for icell = 1:size(fieldnorm,1)
                            fielddist(icell,:) = size(fieldnorm,2)/(2*pi)*circ_dist(2*pi/size(fieldnorm,2)*repmat(fieldPosAll(icell),[1 size(fieldnorm,1)]),2*pi/size(fieldnorm,2)*fieldPosAll);
                        end
                        catch
                            keyboard
                        end
                        popresSpkCorr.s_CA1V1fielddistCA1V1{ianimal,iseries,g} = fielddist(CA1cells,V1cells);
                        popresSpkCorr.s_CA1V1fielddistCA1V1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1fielddistCA1V1{ianimal,iseries,g}(isortCA1,isortV1);
                         
                        s_CA1V1fieldXcorrCA1 = NaN(sum(CA1cells),sum(CA1cells),size(fieldnorm,2));
                        s_CA1V1fieldXcorrV1 = NaN(sum(V1cells),sum(V1cells),size(fieldnorm,2));
                        s_CA1V1fieldXcorrCA1V1 = NaN(sum(CA1cells),sum(V1cells),size(fieldnorm,2));
                        ishift = 0;
                        for xshift = (-floor(size(fieldnorm,2)/2)+1):floor(size(fieldnorm,2)/2)
                            ishift = ishift + 1;
                            s_CA1V1fieldXcorrCA1(:,:,ishift) = fieldnorm(CA1cells,:)*(circshift(fieldnorm(CA1cells,:),xshift,2)')/size(fieldnorm,2);
                            s_CA1V1fieldXcorrCA1(:,:,ishift) = s_CA1V1fieldXcorrCA1(isortCA1,isortCA1,ishift);
                            s_CA1V1fieldXcorrV1(:,:,ishift) = fieldnorm(V1cells,:)*(circshift(fieldnorm(V1cells,:),xshift,2)')/size(fieldnorm,2);
                            s_CA1V1fieldXcorrV1(:,:,ishift) = s_CA1V1fieldXcorrV1(isortV1,isortV1,ishift);
                            s_CA1V1fieldXcorrCA1V1(:,:,ishift) = fieldnorm(CA1cells,:)*(circshift(fieldnorm(V1cells,:),xshift,2)')/size(fieldnorm,2);
                            s_CA1V1fieldXcorrCA1V1(:,:,ishift) = s_CA1V1fieldXcorrCA1V1(isortCA1,isortV1,ishift);
                        end
                        [popresSpkCorr.s_CA1V1fieldXcorrmaxCA1{ianimal,iseries,g},popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1{ianimal,iseries,g}] = max(s_CA1V1fieldXcorrCA1,[],3);
                        [popresSpkCorr.s_CA1V1fieldXcorrmaxV1{ianimal,iseries,g},popresSpkCorr.s_CA1V1fieldXcorrmaxPosV1{ianimal,iseries,g}] = max(s_CA1V1fieldXcorrV1,[],3);
                        [popresSpkCorr.s_CA1V1fieldXcorrmaxCA1V1{ianimal,iseries,g},popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1V1{ianimal,iseries,g}] = max(s_CA1V1fieldXcorrCA1V1,[],3);
                        
                        popresSpkCorr.s_CA1V1fieldcorrCA1{ianimal,iseries,g} = fieldnorm(CA1cells,:)*(fieldnorm(CA1cells,:)')/size(fieldnorm,2);
                        popresSpkCorr.s_CA1V1fieldcorrCA1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1fieldcorrCA1{ianimal,iseries,g}(isortCA1,isortCA1);
                        popresSpkCorr.s_CA1V1fieldcorrV1{ianimal,iseries,g} = fieldnorm(V1cells,:)*(fieldnorm(V1cells,:)')/size(fieldnorm,2);
                        popresSpkCorr.s_CA1V1fieldcorrV1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1fieldcorrV1{ianimal,iseries,g}(isortV1,isortV1);
                        popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g} = fieldnorm(CA1cells,:)*(fieldnorm(V1cells,:)')/size(fieldnorm,2);
                        popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(isortCA1,isortV1);
                                                
                        popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g}(CA1cells,CA1cells);
                        popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}(isortCA1,isortCA1);
                        popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}(logical(eye(size(popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g})))) = NaN;
                        popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g}(V1cells,V1cells);
                        popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}(isortV1,isortV1);
                        popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}(logical(eye(size(popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g})))) = NaN;
                        popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g}(CA1cells,V1cells);
                        popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(isortCA1,isortV1);
%                         popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(logical(eye(size(popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g})))) = NaN;
                        
                        popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1Cross_shuffled{ianimal,iseries,g}(CA1cells,CA1cells,:);
                        popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g}(isortCA1,isortCA1,:);
                        eye3D = [];
                        for ishf = 1:popresSpkCorr.Nshuffle
                            eye3D = cat(3,eye3D,logical(eye(size(popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}))));
                        end
                        popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g}(logical(eye3D)) = NaN;
                        
                        popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1Cross_shuffled{ianimal,iseries,g}(V1cells,V1cells,:);
                        popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g}(isortV1,isortV1,:);
                        eye3D = [];
                        for ishf = 1:popresSpkCorr.Nshuffle
                            eye3D = cat(3,eye3D,logical(eye(size(popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}))));
                        end
                        popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g}(logical(eye3D)) = NaN;
                        popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1Cross_shuffled{ianimal,iseries,g}(CA1cells,V1cells,:);
                        popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}(isortCA1,isortV1,:);
%                         popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}(logical(eye(size(popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g})))) = NaN;
                        
                        popresSpkCorr.s_CA1V1CrossPvalCA1{ianimal,iseries,g} = NaN(size(popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}));
                        popresSpkCorr.s_CA1V1CrossPvalV1{ianimal,iseries,g} = NaN(size(popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}));
                        popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,g} = NaN(size(popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}));
                        popresSpkCorr.s_CA1V1CrossPvalCA1{ianimal,iseries,g} = sum(repmat(abs(popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}),[1 1 popresSpkCorr.Nshuffle]) < abs(popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g}), 3)/popresSpkCorr.Nshuffle;
                        popresSpkCorr.s_CA1V1CrossPvalCA1{ianimal,iseries,g}(logical(eye(size(popresSpkCorr.s_CA1V1CrossPvalCA1{ianimal,iseries,g})))) = NaN;
                        popresSpkCorr.s_CA1V1CrossPvalV1{ianimal,iseries,g} = sum(repmat(abs(popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}),[1 1 popresSpkCorr.Nshuffle]) < abs(popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g}), 3)/popresSpkCorr.Nshuffle;
                        popresSpkCorr.s_CA1V1CrossPvalV1{ianimal,iseries,g}(logical(eye(size(popresSpkCorr.s_CA1V1CrossPvalV1{ianimal,iseries,g})))) = NaN;
                        popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,g} = sum(repmat(abs(popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}),[1 1 popresSpkCorr.Nshuffle]) < abs(popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}), 3)/popresSpkCorr.Nshuffle;
                        
                        poscorr = popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(:) > 0;
                        signi = popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,g}(:) <= popresSpkCorr.PthCrossCA1V1;
                        valid = ~isnan(popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(:)) & ~isnan(popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(:));
                        fieldcorr = popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(valid & signi & poscorr);%popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1V1{ianimal,iseries,g}(valid)+50;%popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1V1{ianimal,iseries,g}(valid);%popresSpkCorr.s_CA1V1fieldXcorrmaxCA1V1{ianimal,iseries,g}(valid);%
%                         fieldcorr(fieldcorr==min(fieldcorr(:))) = -1;%1;
%                         fieldcorr(fieldcorr==max(fieldcorr(:))) = 1;%size(fieldnorm,2);
                        fieldcorr = normalise1var(fieldcorr, popresSpkCorr.nbinscorr, [], [-1 1]);%normalise1var(fieldcorrPos, popresSpkCorr.nbinscorr, [], [1 size(fieldnorm,2)]);%
                        fieldcorr = fieldcorr - 1;
                        
                        fieldcorrPos = popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1V1{ianimal,iseries,g}(valid & signi & poscorr);%popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(valid);%popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1V1{ianimal,iseries,g}(valid)+50;%popresSpkCorr.s_CA1V1fieldXcorrmaxCA1V1{ianimal,iseries,g}(valid);%
%                         fieldcorrPos(fieldcorrPos==min(fieldcorrPos(:))) = 1;
%                         fieldcorrPos(fieldcorrPos==max(fieldcorrPos(:))) = size(fieldnorm,2);
                        fieldcorrPos = normalise1var(fieldcorrPos, popresSpkCorr.nbinscorr, [], [1 size(fieldnorm,2)]);
                        fieldcorrPos = fieldcorrPos - 1;
                        
                        fieldPos = popresSpkCorr.s_CA1V1fielddistCA1V1{ianimal,iseries,g}(valid & signi & poscorr)+50;%popresSpkCorr.s_CA1V1fieldXcorrmaxCA1V1{ianimal,iseries,g}(valid);%
%                         fieldPos(fieldPos==min(fieldPos(:))) = 0;
%                         fieldPos(fieldPos==max(fieldPos(:))) = size(fieldnorm,2);
                        fieldPos = normalise1var(fieldPos, popresSpkCorr.nbinscorr, [], [0 size(fieldnorm,2)]);
                        fieldPos = fieldPos - 1;
                        
                        noisecorr = popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(valid & signi & poscorr);
                        respcorr = popresSpkCorr.s_CA1V1fieldXcorrmaxCA1V1{ianimal,iseries,g}(valid & signi & poscorr);%popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(valid);%
                        popresSpkCorr.s_CA1V1CrossAveCA1V1_Corr{ianimal,iseries,g} = fast1Dmap(fieldcorr,noisecorr,1,1,popresSpkCorr.nbinscorr/1,false,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1CrossAveCA1V1_dCorr{ianimal,iseries,g} = fast1Dmap(fieldcorrPos,noisecorr,1,1,popresSpkCorr.nbinscorr/1,true,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1CrossAveCA1V1_dPeak{ianimal,iseries,g} = fast1Dmap(fieldPos,noisecorr,1,1,popresSpkCorr.nbinscorr/1,false,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1fieldXcorrAveCA1V1_Corr{ianimal,iseries,g} = fast1Dmap(fieldcorr,respcorr,1,1,popresSpkCorr.nbinscorr/1,false,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1fieldXcorrAveCA1V1_dCorr{ianimal,iseries,g} = fast1Dmap(fieldcorrPos,respcorr,1,1,popresSpkCorr.nbinscorr/1,true,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1fieldXcorrAveCA1V1_dPeak{ianimal,iseries,g} = fast1Dmap(fieldPos,respcorr,1,1,popresSpkCorr.nbinscorr/1,false,popresSpkCorr.nbinscorr);%true);
                        
                        
                        crossCA1V1_shuffled = popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}(:,:,1);
                        noisecorr_shuffled = crossCA1V1_shuffled(valid & signi & poscorr);
                        popresSpkCorr.s_CA1V1CrossAveCA1V1_Corr_shuffled{ianimal,iseries,g} = fast1Dmap(fieldcorr,noisecorr_shuffled,1,1,popresSpkCorr.nbinscorr/1,false,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1CrossAveCA1V1_dCorr_shuffled{ianimal,iseries,g} = fast1Dmap(fieldcorrPos,noisecorr_shuffled,1,1,popresSpkCorr.nbinscorr/1,true,popresSpkCorr.nbinscorr);%true);
                        popresSpkCorr.s_CA1V1CrossAveCA1V1_dPeak_shuffled{ianimal,iseries,g} = fast1Dmap(fieldPos,noisecorr_shuffled,1,1,popresSpkCorr.nbinscorr/1,false,popresSpkCorr.nbinscorr);%true);
                        
                        popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g}(:,:,1);
                        popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g}(:,:,1);
                        popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g} = popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}(:,:,1);
                    end
                end
            end
        end
    end
end

%now we do the grand averages by averaging across sessions
nbXbins = popresSpkCorr.nbinscorr;
for g = [2 1 3]
    popresSpkCorr.CA1V1CrossCA1V1{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1CrossCA1{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1CrossV1{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1CrossCA1V1_shuffled{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1CrossCA1_shuffled{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1CrossV1_shuffled{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1fieldcorrCA1{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1fieldcorrV1{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1fieldcorrCA1V1{g} = cell(nbXbins,nbXbins);
    popresSpkCorr.CA1V1CrossCA1_All{g} = [];
    popresSpkCorr.CA1V1CrossPvalCA1_All{g} = [];
    popresSpkCorr.CA1V1CrossCA1_shuffled_All{g} = [];
    popresSpkCorr.CA1V1CrossV1_All{g} = [];
    popresSpkCorr.CA1V1CrossPvalV1_All{g} = [];
    popresSpkCorr.CA1V1CrossV1_shuffled_All{g} = [];
    popresSpkCorr.CA1V1CrossCA1V1_All{g} = [];
    popresSpkCorr.CA1V1CrossPvalCA1V1_All{g} = [];
    popresSpkCorr.CA1V1CrossCA1V1_shuffled_All{g} = [];
    popresSpkCorr.CA1V1fieldcorrCA1V1_All{g} = [];
    popresSpkCorr.CA1V1fieldXcorrmaxPosCA1V1_All{g} = [];
    popresSpkCorr.CA1V1fieldXcorrmaxCA1V1_All{g} = [];
    popresSpkCorr.CA1V1fielddistCA1V1_All{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled{g} = [];
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr_SE{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled_SE{g} = [];
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr_SE{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g} = [];
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_SE{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled_SE{g} = [];
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr_SE{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g} = [];
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_SE{g} = [];
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled_SE{g} = [];
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak_SE{g} = [];
    for ianimal = 1:size(popresSpkCorr.s_CA1V1Cross,1)
        for iseries = 1:size(popresSpkCorr.s_CA1V1Cross,2)
            if ~isempty(popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g})
               if sum(isnan(popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g}(:))) ~= numel(popresSpkCorr.s_CA1V1Cross{ianimal,iseries,g})
                   fposCA1 = popresSpkCorr.s_fposCA1{ianimal,iseries,g};
                   fposV1 = popresSpkCorr.s_fposV1{ianimal,iseries,g};
                   for i=1:nbXbins
                       for j=1:nbXbins
                           goodCA1cells_i = double(fposCA1>=(i-1)*100/popresSpkCorr.nbinscorr & fposCA1<=i*100/popresSpkCorr.nbinscorr);
                           goodV1cells_i = double(fposV1>=(i-1)*100/popresSpkCorr.nbinscorr & fposV1<=i*100/popresSpkCorr.nbinscorr);
                           goodCA1cells_j = double(fposCA1>=(j-1)*100/popresSpkCorr.nbinscorr & fposCA1<=j*100/popresSpkCorr.nbinscorr);
                           goodV1cells_j = double(fposV1>=(j-1)*100/popresSpkCorr.nbinscorr & fposV1<=j*100/popresSpkCorr.nbinscorr);
                           CA1poscorr = popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g} > 0;
                           CA1signi = popresSpkCorr.s_CA1V1CrossPvalCA1{ianimal,iseries,g}<=popresSpkCorr.PthCrossCA1V1;
                           V1poscorr = popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g} > 0;
                           V1signi = popresSpkCorr.s_CA1V1CrossPvalV1{ianimal,iseries,g}<=popresSpkCorr.PthCrossCA1V1;
                           CA1V1poscorr = popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g} > 0;
                           CA1V1signi = popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,g}<=popresSpkCorr.PthCrossCA1V1;
                           
                           mtemp = popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}(goodCA1cells_i'*goodCA1cells_j & CA1poscorr & CA1signi);
                           popresSpkCorr.CA1V1CrossCA1{g}{i,j} = cat(1,popresSpkCorr.CA1V1CrossCA1{g}{i,j},mtemp(:));
                           mtemp = popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}(goodV1cells_i'*goodV1cells_j & V1poscorr & V1signi);
                           popresSpkCorr.CA1V1CrossV1{g}{i,j} = cat(1,popresSpkCorr.CA1V1CrossV1{g}{i,j},mtemp(:));
                           mtemp = popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(goodCA1cells_i'*goodV1cells_j & CA1V1poscorr & CA1V1signi);
                           popresSpkCorr.CA1V1CrossCA1V1{g}{i,j} = cat(1,popresSpkCorr.CA1V1CrossCA1V1{g}{i,j},mtemp(:));
                           
                           mtemp = popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g}(goodCA1cells_i'*goodCA1cells_j & CA1poscorr & CA1signi);
                           popresSpkCorr.CA1V1CrossCA1_shuffled{g}{i,j} = cat(1,popresSpkCorr.CA1V1CrossCA1_shuffled{g}{i,j},mtemp(:));
                           mtemp = popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g}(goodV1cells_i'*goodV1cells_j & V1poscorr & V1signi);
                           popresSpkCorr.CA1V1CrossV1_shuffled{g}{i,j} = cat(1,popresSpkCorr.CA1V1CrossV1_shuffled{g}{i,j},mtemp(:));
                           mtemp = popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}(goodCA1cells_i'*goodV1cells_j & CA1V1poscorr & CA1V1signi);
                           popresSpkCorr.CA1V1CrossCA1V1_shuffled{g}{i,j} = cat(1,popresSpkCorr.CA1V1CrossCA1V1_shuffled{g}{i,j},mtemp(:));
                           
                           mtemp = popresSpkCorr.s_CA1V1fieldcorrCA1{ianimal,iseries,g}(goodCA1cells_i'*goodCA1cells_j & CA1poscorr & CA1signi);
                           popresSpkCorr.CA1V1fieldcorrCA1{g}{i,j} = cat(1,popresSpkCorr.CA1V1fieldcorrCA1{g}{i,j},mtemp(:));
                           mtemp = popresSpkCorr.s_CA1V1fieldcorrV1{ianimal,iseries,g}(goodV1cells_i'*goodV1cells_j & V1poscorr & V1signi);
                           popresSpkCorr.CA1V1fieldcorrV1{g}{i,j} = cat(1,popresSpkCorr.CA1V1fieldcorrV1{g}{i,j},mtemp(:));
                           mtemp = popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(goodCA1cells_i'*goodV1cells_j & CA1V1poscorr & CA1V1signi);
                           popresSpkCorr.CA1V1fieldcorrCA1V1{g}{i,j} = cat(1,popresSpkCorr.CA1V1fieldcorrCA1V1{g}{i,j},mtemp(:));
                       end
                   end
                   popresSpkCorr.CA1V1CrossCA1_All{g} = cat(1,popresSpkCorr.CA1V1CrossCA1_All{g},popresSpkCorr.s_CA1V1CrossCA1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossPvalCA1_All{g} = cat(1,popresSpkCorr.CA1V1CrossPvalCA1_All{g},popresSpkCorr.s_CA1V1CrossPvalCA1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossCA1_shuffled_All{g} = cat(1,popresSpkCorr.CA1V1CrossCA1_shuffled_All{g},popresSpkCorr.s_CA1V1CrossCA1_shuffled{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossV1_All{g} = cat(1,popresSpkCorr.CA1V1CrossV1_All{g},popresSpkCorr.s_CA1V1CrossV1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossPvalV1_All{g} = cat(1,popresSpkCorr.CA1V1CrossPvalV1_All{g},popresSpkCorr.s_CA1V1CrossPvalV1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossV1_shuffled_All{g} = cat(1,popresSpkCorr.CA1V1CrossV1_shuffled_All{g},popresSpkCorr.s_CA1V1CrossV1_shuffled{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossCA1V1_All{g} = cat(1,popresSpkCorr.CA1V1CrossCA1V1_All{g},popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossPvalCA1V1_All{g} = cat(1,popresSpkCorr.CA1V1CrossPvalCA1V1_All{g},popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1CrossCA1V1_shuffled_All{g} = cat(1,popresSpkCorr.CA1V1CrossCA1V1_shuffled_All{g},popresSpkCorr.s_CA1V1CrossCA1V1_shuffled{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1fieldcorrCA1V1_All{g} = cat(1,popresSpkCorr.CA1V1fieldcorrCA1V1_All{g},popresSpkCorr.s_CA1V1fieldcorrCA1V1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1fieldXcorrmaxPosCA1V1_All{g} = cat(1,popresSpkCorr.CA1V1fieldXcorrmaxPosCA1V1_All{g},popresSpkCorr.s_CA1V1fieldXcorrmaxPosCA1V1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1fieldXcorrmaxCA1V1_All{g} = cat(1,popresSpkCorr.CA1V1fieldXcorrmaxCA1V1_All{g},popresSpkCorr.s_CA1V1fieldXcorrmaxCA1V1{ianimal,iseries,g}(:));
                   popresSpkCorr.CA1V1fielddistCA1V1_All{g} = cat(1,popresSpkCorr.CA1V1fielddistCA1V1_All{g},popresSpkCorr.s_CA1V1fielddistCA1V1{ianimal,iseries,g}(:)+50);
                   
                   popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g} = cat(2,popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g},popresSpkCorr.s_CA1V1CrossAveCA1V1_Corr{ianimal,iseries,g});
                   popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled{g} = cat(2,popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled{g},popresSpkCorr.s_CA1V1CrossAveCA1V1_Corr_shuffled{ianimal,iseries,g});
                   popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr{g} = cat(2,popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr{g},popresSpkCorr.s_CA1V1fieldXcorrAveCA1V1_Corr{ianimal,iseries,g});
                   
                   popresSpkCorr.CA1V1CrossAveCA1V1_dCorr{g} = cat(2,popresSpkCorr.CA1V1CrossAveCA1V1_dCorr{g},popresSpkCorr.s_CA1V1CrossAveCA1V1_dCorr{ianimal,iseries,g});
                   popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g} = cat(2,popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g},popresSpkCorr.s_CA1V1CrossAveCA1V1_dCorr_shuffled{ianimal,iseries,g});
                   popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr{g} = cat(2,popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr{g},popresSpkCorr.s_CA1V1fieldXcorrAveCA1V1_dCorr{ianimal,iseries,g});
                   
                   popresSpkCorr.CA1V1CrossAveCA1V1_dPeak{g} = cat(2,popresSpkCorr.CA1V1CrossAveCA1V1_dPeak{g},popresSpkCorr.s_CA1V1CrossAveCA1V1_dPeak{ianimal,iseries,g});
                   popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g} = cat(2,popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g},popresSpkCorr.s_CA1V1CrossAveCA1V1_dPeak_shuffled{ianimal,iseries,g});
                   popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak{g} = cat(2,popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak{g},popresSpkCorr.s_CA1V1fieldXcorrAveCA1V1_dPeak{ianimal,iseries,g});
                end
            end
        end
    end
    for i=1:nbXbins
        for j=1:nbXbins
            popresSpkCorr.CA1V1CrossCA1{g}{i,j} = nanmean(popresSpkCorr.CA1V1CrossCA1{g}{i,j},1);
            popresSpkCorr.CA1V1CrossV1{g}{i,j} = nanmean(popresSpkCorr.CA1V1CrossV1{g}{i,j},1);
            popresSpkCorr.CA1V1CrossCA1V1{g}{i,j} = nanmean(popresSpkCorr.CA1V1CrossCA1V1{g}{i,j},1);
            
            popresSpkCorr.CA1V1CrossCA1_shuffled{g}{i,j} = nanmean(popresSpkCorr.CA1V1CrossCA1_shuffled{g}{i,j},1);
            popresSpkCorr.CA1V1CrossV1_shuffled{g}{i,j} = nanmean(popresSpkCorr.CA1V1CrossV1_shuffled{g}{i,j},1);
            popresSpkCorr.CA1V1CrossCA1V1_shuffled{g}{i,j} = nanmean(popresSpkCorr.CA1V1CrossCA1V1_shuffled{g}{i,j},1);
            
            popresSpkCorr.CA1V1fieldcorrCA1{g}{i,j} = nanmean(popresSpkCorr.CA1V1fieldcorrCA1{g}{i,j},1);
            popresSpkCorr.CA1V1fieldcorrV1{g}{i,j} = nanmean(popresSpkCorr.CA1V1fieldcorrV1{g}{i,j},1);
            popresSpkCorr.CA1V1fieldcorrCA1V1{g}{i,j} = nanmean(popresSpkCorr.CA1V1fieldcorrCA1V1{g}{i,j},1);
        end
    end
    popresSpkCorr.CA1V1CrossCA1{g} = cell2mat(popresSpkCorr.CA1V1CrossCA1{g});
    popresSpkCorr.CA1V1CrossV1{g} = cell2mat(popresSpkCorr.CA1V1CrossV1{g});
    popresSpkCorr.CA1V1CrossCA1V1{g} = cell2mat(popresSpkCorr.CA1V1CrossCA1V1{g});
    
    popresSpkCorr.CA1V1CrossCA1_shuffled{g} = cell2mat(popresSpkCorr.CA1V1CrossCA1_shuffled{g});
    popresSpkCorr.CA1V1CrossV1_shuffled{g} = cell2mat(popresSpkCorr.CA1V1CrossV1_shuffled{g});
    popresSpkCorr.CA1V1CrossCA1V1_shuffled{g} = cell2mat(popresSpkCorr.CA1V1CrossCA1V1_shuffled{g});
    
    popresSpkCorr.CA1V1fieldcorrCA1{g} = cell2mat(popresSpkCorr.CA1V1fieldcorrCA1{g});
    popresSpkCorr.CA1V1fieldcorrV1{g} = cell2mat(popresSpkCorr.CA1V1fieldcorrV1{g});
    popresSpkCorr.CA1V1fieldcorrCA1V1{g} = cell2mat(popresSpkCorr.CA1V1fieldcorrCA1V1{g});
    
    CrossAve_Corr = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g},2);
    CrossAve_Corr_shuffled = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled{g},2);
    fieldCrossAve_Corr = nanmean(popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr{g},2);
    CrossAve_Corr_std = 0;
    CrossAve_Corr_shuffled_std = 0;
    fieldCrossAve_Corr_std = 0;
    CrossAve_dCorr = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dCorr{g},2);
    CrossAve_dCorr_shuffled = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g},2);
    fieldCrossAve_dCorr = nanmean(popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr{g},2);
    CrossAve_dCorr_std = 0;
    CrossAve_dCorr_shuffled_std = 0;
    fieldCrossAve_dCorr_std = 0;
    CrossAve_dPeak = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dPeak{g},2);
    CrossAve_dPeak_shuffled = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g},2);
    fieldCrossAve_dPeak = nanmean(popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak{g},2);
    CrossAve_dPeak_std = 0;
    CrossAve_dPeak_shuffled_std = 0;
    fieldCrossAve_dPeak_std = 0;
    kfold = size(popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g},2);
    allsessions = 1:size(popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g},2);
    for k = 1:kfold
        CrossAve_Corr_iter = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g}(:,~ismember(allsessions,k)),2);
        CrossAve_Corr_shuffled_iter = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled{g}(:,~ismember(allsessions,k)),2);
        fieldCrossAve_Corr_iter = nanmean(popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr{g}(:,~ismember(allsessions,k)),2);
        CrossAve_Corr_std = CrossAve_Corr_std + (kfold - 1)/kfold*(CrossAve_Corr_iter - CrossAve_Corr).^2;
        CrossAve_Corr_shuffled_std = CrossAve_Corr_shuffled_std + (kfold - 1)/kfold*(CrossAve_Corr_shuffled_iter - CrossAve_Corr_shuffled).^2;
        fieldCrossAve_Corr_std = fieldCrossAve_Corr_std + (kfold - 1)/kfold*(fieldCrossAve_Corr_iter - fieldCrossAve_Corr).^2;
        
        CrossAve_dCorr_iter = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dCorr{g}(:,~ismember(allsessions,k)),2);
        CrossAve_dCorr_shuffled_iter = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g}(:,~ismember(allsessions,k)),2);
        fieldCrossAve_dCorr_iter = nanmean(popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr{g}(:,~ismember(allsessions,k)),2);
        CrossAve_dCorr_std = CrossAve_dCorr_std + (kfold - 1)/kfold*(CrossAve_dCorr_iter - CrossAve_dCorr).^2;
        CrossAve_dCorr_shuffled_std = CrossAve_dCorr_shuffled_std + (kfold - 1)/kfold*(CrossAve_dCorr_shuffled_iter - CrossAve_dCorr_shuffled).^2;
        fieldCrossAve_dCorr_std = fieldCrossAve_dCorr_std + (kfold - 1)/kfold*(fieldCrossAve_dCorr_iter - fieldCrossAve_dCorr).^2;
        
        CrossAve_dPeak_iter = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dPeak{g}(:,~ismember(allsessions,k)),2);
        CrossAve_dPeak_shuffled_iter = nanmean(popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g}(:,~ismember(allsessions,k)),2);
        fieldCrossAve_dPeak_iter = nanmean(popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak{g}(:,~ismember(allsessions,k)),2);
        CrossAve_dPeak_std = CrossAve_dPeak_std + (kfold - 1)/kfold*(CrossAve_dPeak_iter - CrossAve_dPeak).^2;
        CrossAve_dPeak_shuffled_std = CrossAve_dPeak_shuffled_std + (kfold - 1)/kfold*(CrossAve_dPeak_shuffled_iter - CrossAve_dPeak_shuffled).^2;
        fieldCrossAve_dPeak_std = fieldCrossAve_dPeak_std + (kfold - 1)/kfold*(fieldCrossAve_dPeak_iter - fieldCrossAve_dPeak).^2;
    end
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr_SE{g} = sqrt(CrossAve_Corr_std);
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled_SE{g} = sqrt(CrossAve_Corr_shuffled_std);
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr_SE{g} = sqrt(fieldCrossAve_Corr_std);
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr{g} = CrossAve_Corr;
    popresSpkCorr.CA1V1CrossAveCA1V1_Corr_shuffled{g} = CrossAve_Corr_shuffled;
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_Corr{g} = fieldCrossAve_Corr;
    
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_SE{g} = sqrt(CrossAve_dCorr_std);
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled_SE{g} = sqrt(CrossAve_dCorr_shuffled_std);
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr_SE{g} = sqrt(fieldCrossAve_dCorr_std);
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr{g} = CrossAve_dCorr;
    popresSpkCorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g} = CrossAve_dCorr_shuffled;
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dCorr{g} = fieldCrossAve_dCorr;
    
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_SE{g} = sqrt(CrossAve_dPeak_std);
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled_SE{g} = sqrt(CrossAve_dPeak_shuffled_std);
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak_SE{g} = sqrt(fieldCrossAve_dPeak_std);
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak{g} = CrossAve_dPeak;
    popresSpkCorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g} = CrossAve_dPeak_shuffled;
    popresSpkCorr.CA1V1fieldXcorrAveCA1V1_dPeak{g} = fieldCrossAve_dPeak;
end
save(savedfilename_popresCorr, 'popresSpkCorr','-v7.3');
end


%do this before plotting distri of corrcoeff in order to remove duplicates
% nSignicorrV1 = [];
% FinterneuronV1 = [];
% corrAll = [];
% pvalAll = [];
% for ianimal = 1:size(popresSpkCorr.s_CA1V1CrossPvalCA1V1,1)
% for iseries = 1:size(popresSpkCorr.s_CA1V1CrossPvalCA1V1,2)
% if ~isempty(popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,2})
% CrossCorr = popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,2} + tril(inf(size(popresSpkCorr.s_CA1V1CrossCA1V1{ianimal,iseries,2})));
% Pval = popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,2} + tril(NaN(size(popresSpkCorr.s_CA1V1CrossPvalCA1V1{ianimal,iseries,2})));
% corrAll = [corrAll CrossCorr(~isnan(Pval))'];
% pvalAll = [pvalAll Pval(~isnan(Pval))'];
% nSignicorrV1 = [nSignicorrV1 nansum(Pval<=0.05,1)];
% FinterneuronV1 = [FinterneuronV1 popresSpkCorr.V1min2maxSpkwf{ianimal,iseries}(:)'<18];
% end
% end
% end
% nSignicorrV1_blank = [];
% FinterneuronV1_blank = [];
% corrAll_blank = [];
% pvalAll_blank = [];
% for ianimal = 1:size(popresSpkCorr_blank.s_CA1V1CrossPvalCA1V1,1)
% for iseries = 1:size(popresSpkCorr_blank.s_CA1V1CrossPvalCA1V1,2)
% if ~isempty(popresSpkCorr_blank.s_CA1V1CrossPvalCA1V1{ianimal,iseries,2})
% CrossCorr_blank = popresSpkCorr_blank.s_CA1V1CrossCA1V1{ianimal,iseries,2} + tril(inf(size(popresSpkCorr_blank.s_CA1V1CrossCA1V1{ianimal,iseries,2})));
% Pval_blank = popresSpkCorr_blank.s_CA1V1CrossPvalCA1V1{ianimal,iseries,2} + tril(NaN(size(popresSpkCorr_blank.s_CA1V1CrossPvalCA1V1{ianimal,iseries,2})));
% corrAll_blank = [corrAll_blank CrossCorr_blank(~isnan(Pval_blank))'];
% pvalAll_blank = [pvalAll_blank Pval_blank(~isnan(Pval_blank))'];
% nSignicorrV1_blank = [nSignicorrV1_blank nansum(Pval_blank<=0.05,1)];
% FinterneuronV1_blank = [FinterneuronV1_blank popresSpkCorr_blank.V1min2maxSpkwf{ianimal,iseries}(:)'<18];
% end
% end
% end
% figure;
% subplot(1,3,1)
% scatter(corrAll(~isnan(corrAll_blank)),corrAll_blank(~isnan(corrAll_blank)),'MarkerEdgeAlpha',0.6,'MarkerEdgeColor','k')
% set(gca,'PlotBoxAspectRatio', [1 1 1])
% title('Blank vs. Corridor')
% subplot(1,3,2)
% histogram(corrAll(~isnan(corrAll_blank)),-0.4:0.01:0.4,'EdgeColor','none','FaceColor',[0.5 0.5 0.5]);
% hold on;histogram(corrAll(~isnan(corrAll_blank) & pvalAll<=0.05),-0.4:0.01:0.4,'EdgeColor','none','FaceColor','k');
% set(gca,'PlotBoxAspectRatio', [1 1 1])
% title('Corridor')
% subplot(1,3,3)
% histogram(corrAll_blank(~isnan(corrAll_blank)),-0.4:0.01:0.4,'EdgeColor','none','FaceColor',[0.5 0.5 0.5]);
% hold on;histogram(corrAll_blank(~isnan(corrAll_blank) & pvalAll_blank<=0.05),-0.4:0.01:0.4,'EdgeColor','none','FaceColor','k');
% set(gca,'PlotBoxAspectRatio', [1 1 1])
% title('Blank')

%use the following to visualize the data until there's a GUI or fig
%associated
% for g = [2 1 3]
% for tlag = 1:1
% count = 0;
% CovMat_Ave{tlag,g} = 0;
% CovMat_Shuffled{tlag,g} = 0;
% for ianimal = 1:size(popresSpkCorr.s_CA1V1Cross,1)
% for iseries = 1:size(popresSpkCorr.s_CA1V1Cross,2)
% if ~isempty(popresSpkCorr.s_CA1V1Cross{ianimal,iseries,tlag,g})
% %This will average across sessions than across all trials
% mat = popresSpkCorr.s_CA1V1Joint{ianimal,iseries,tlag,g}./repmat(popresSpkCorr.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% if sum(isnan(mat)) == 0
% CovMat_Ave{tlag,g} = CovMat_Ave{tlag,g} + mat;
% CovMat_Shuffled{tlag,g} = CovMat_Shuffled{tlag,g} + popresSpkCorr.s_CA1V1Joint_shuffled{ianimal,iseries,tlag,g}./repmat(popresSpkCorr.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% count = count+ 1;%popresSpkCorr.s_nDataPoints{ianimal,iseries,tlag,g};
% end
% end
% end
% end
% CovMat_Ave{tlag,g} = CovMat_Ave{tlag,g}/count;%./repmat(count,[1 1 100 100]);
% CovMat_Shuffled{tlag,g} = CovMat_Shuffled{tlag,g}/count;%./repmat(count,[1 1 100 100]);
% end
% end
% 
% figure;
% mat1 = [];
% mat2 = [];
% for g = 1:3
% mat1{g} = 0;
% mat2{g} = 0;
% end
% for xx = 1:100
% for g = 1:3
% mattemp1 = CovMat_Ave{1,g}(xx,:,:,:);
% mattemp2 = CovMat_Shuffled{1,g}(xx,:,:,:);
% mattemp1 = squeeze(nanmean(nanmean(mattemp1,1),2));
% mattemp2 = squeeze(nanmean(nanmean(mattemp2,1),2));
% mat1{g} = (mat1{g}*(xx-1) + mattemp1)/xx;
% mat2{g} = (mat2{g}*(xx-1) + mattemp2)/xx;
% subplot(2,3,g)
% imagesc(mat1{g});
% set(gca,'Clim',[-1 1],'Ydir','normal');
% subplot(2,3,g+3)
% imagesc(mat2{g});
% set(gca,'Clim',[-1 1],'Ydir','normal');
% xlabel(num2str(xx))
% end
% pause
% end


%or alternatively
% for g = [2 1 3]
% count = 0;
% crosscount = 0;
% CovMat_Ave{g} = 0;
% CovLagMat_Ave{g} = 0;
% CrossMat_Ave{g} = 0;
% CovMat_Shuffled{g} = 0;
% CovLagMat_Shuffled{g} = 0;
% CrossMat_Shuffled{g} = 0;
% CovMat_stdCA1{g} = 0;
% CovMat_stdV1{g} = 0;
% for ianimal = 1:size(popresSpkCorr.s_CA1V1Cross,1)
% for iseries = 1:size(popresSpkCorr.s_CA1V1Cov,2)
% if ~isempty(popresSpkCorr.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresSpkCorr.s_CA1V1Joint{ianimal,iseries,g};
% matshf = popresSpkCorr.s_CA1V1Joint_shuffled{ianimal,iseries,g};
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% matcross(xx,spd,:,:) = squeeze(mat(xx,spd,:,:))./(squeeze(popresSpkCorr.s_V1marginal{ianimal,iseries,g}(xx,spd,:))*squeeze(popresSpkCorr.s_CA1marginal{ianimal,iseries,g}(xx,spd,:))').^0.5;
% matcross_shf(xx,spd,:,:) = squeeze(matshf(xx,spd,:,:))./(squeeze(popresSpkCorr.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))*squeeze(popresSpkCorr.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))').^0.5;
% end
% end
% if sum(isnan(mat)) == 0
% CovMat_Ave{g} = CovMat_Ave{g} + mat;
% CrossMat_Ave{g} = CrossMat_Ave{g} + matcross;
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g} + matcross_shf;
% CovLagMat_Ave{g} = CovLagMat_Ave{g} + popresCorr_250.s_CA1V1Jointdiag{ianimal,iseries,g};
% crosscount = crosscount + 1;
% CovMat_stdCA1{g} = CovMat_stdCA1{g} + popresCorr_250.s_CA1marginal{ianimal,iseries,g};
% CovMat_stdV1{g} = CovMat_stdV1{g} + popresCorr_250.s_V1marginal{ianimal,iseries,g};
% CovMat_Shuffled{g} = CovMat_Shuffled{g} + popresCorr_250.s_CA1V1Joint_shuffled{ianimal,iseries,g};
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g} + popresCorr_250.s_CA1V1Jointdiag_shuffled{ianimal,iseries,g};
% count = count+ popresCorr_250.s_nDataPoints{ianimal,iseries,g};
% end
% end
% end
% end
% CovMat_Ave{g} = CovMat_Ave{g}./repmat(count,[1 1 100 100]);
% CovLagMat_Ave{g} = CovLagMat_Ave{g}./repmat(count,[1 1 241 100]);
% CrossMat_Ave{g} = CrossMat_Ave{g}/crosscount;
% CovMat_Shuffled{g} = CovMat_Shuffled{g}./repmat(count,[1 1 100 100]);
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g}./repmat(count,[1 1 241 100]);
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g}/crosscount;
% CovMat_stdCA1{g} = CovMat_stdCA1{g}./repmat(count,[1 1 100]);
% CovMat_stdV1{g} = CovMat_stdV1{g}./repmat(count,[1 1 100]);
% end
% 
% 
% same as above but averaging across sessions
% for g = [2 1 3]
% count = 0;
% crosscount = 0;
% CovMat_Ave{g} = 0;
% CovLagMat_Ave{g} = 0;
% CrossMat_Ave{g} = 0;
% CovMat_Shuffled{g} = 0;
% CovLagMat_Shuffled{g} = 0;
% CrossMat_Shuffled{g} = 0;
% CovMat_stdCA1{g} = 0;
% CovMat_stdV1{g} = 0;
% for ianimal = 1:size(popresSpkCorr.s_CA1V1Cov,1)
% for iseries = 1:size(popresSpkCorr.s_CA1V1Cov,2)
% if ~isempty(popresSpkCorr.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresSpkCorr.s_CA1V1Joint{ianimal,iseries,g}./repmat(popresSpkCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matshf = popresSpkCorr.s_CA1V1Joint_shuffled{ianimal,iseries,g}./repmat(popresSpkCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% for ieye = 1:size(mat,3)
% matcross(xx,spd,ieye,:,:) = squeeze(popresSpkCorr.s_CA1V1Joint{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresSpkCorr.s_V1marginal{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresSpkCorr.s_CA1marginal{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% matcross_shf(xx,spd,ieye,:,:) = squeeze(popresSpkCorr.s_CA1V1Joint_shuffled{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresSpkCorr.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresSpkCorr.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% end
% end
% end
% if sum(isnan(mat)) == 0
% CovMat_Ave{g} = CovMat_Ave{g} + mat;
% CrossMat_Ave{g} = CrossMat_Ave{g} + matcross;
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g} + matcross_shf;
% CovLagMat_Ave{g} = CovLagMat_Ave{g} + popresSpkCorr.s_CA1V1Jointdiag{ianimal,iseries,g}./repmat(popresSpkCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
% crosscount = crosscount + 1;
% CovMat_stdCA1{g} = CovMat_stdCA1{g} + popresSpkCorr.s_CA1marginal{ianimal,iseries,g};
% CovMat_stdV1{g} = CovMat_stdV1{g} + popresSpkCorr.s_V1marginal{ianimal,iseries,g};
% CovMat_Shuffled{g} = CovMat_Shuffled{g} + matshf;
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g} + popresSpkCorr.s_CA1V1Jointdiag_shuffled{ianimal,iseries,g}./repmat(popresSpkCorr.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
% count = count+ 1;%popresCorr_250.s_nDataPoints{ianimal,iseries,g};
% end
% end
% end
% end
% CovMat_Ave{g} = CovMat_Ave{g}./repmat(count,[100 3 3 100 100]);
% CovLagMat_Ave{g} = CovLagMat_Ave{g}./repmat(count,[100 3 3 241 100]);
% CrossMat_Ave{g} = CrossMat_Ave{g}/crosscount;
% CovMat_Shuffled{g} = CovMat_Shuffled{g}./repmat(count,[100 3 3 100 100]);
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g}./repmat(count,[100 3 3 241 100]);
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g}/crosscount;
% CovMat_stdCA1{g} = CovMat_stdCA1{g}./repmat(count,[100 3 3 100]);
% CovMat_stdV1{g} = CovMat_stdV1{g}./repmat(count,[100 3 3 100]);
% end