function popresCorr1D = BatchJointCorrelation3(batch2p,ErrType,AveType,ShiftType)
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
    ErrType = 'Error';%'DecPos';%Bias;%'BlankPos';%'BlankErr';%
end
if nargin < 3
    AveType = 'full';%'Pos';%
end
if nargin < 4
    ShiftType = 'Space';%'Time';%
end
popresCorr1D.ErrType = ErrType;
popresCorr1D.AveType = AveType;
popresCorr1D.ShiftType = ShiftType;
popresCorr1D.Tsmthwin = 15;%250;%250;%150;%300;%40;%120;%50
popresCorr1D.Tsmthwin_dec = 250;%
popresCorr1D.nDecbins = 100;
popresCorr1D.Xsmthwin = 4;%2%1;%
popresCorr1D.SpdSmthWin = popresCorr1D.Tsmthwin_dec;
popresCorr1D.SpeedThreshold = 5;
popresCorr1D.nspeedbins = 5;
popresCorr1D.neyebins = 1;
popresCorr1D.nthetaphsbins = 1;%1;%
popresCorr1D.nphsbins = 1;
popresCorr1D.cellstr = 'goodonly';%'All_50bins';%'goodonly';%'goodonly_unwrapped';%'goodonly';%'All';%
filesuffix_EXP = ['Twin' num2str(popresCorr1D.Tsmthwin) '_' 'Xwin' num2str(popresCorr1D.Xsmthwin) '_' 'spdth' num2str(popresCorr1D.SpeedThreshold) '_' 'Decwin' num2str(popresCorr1D.Tsmthwin_dec) '_' 'nDecbins' num2str(popresCorr1D.nDecbins) '_' num2str(popresCorr1D.nspeedbins) 'speedbins' '_' num2str(popresCorr1D.neyebins) 'eyebins' '_' num2str(popresCorr1D.nthetaphsbins) 'thetabins' '_' popresCorr1D.cellstr];
disp(filesuffix_EXP);
disp(['correlating ' popresCorr1D.ErrType])

popresCorr1D.sampleRate = 60;
popresCorr1D.nSpdbins = 5;%1;%
popresCorr1D.nEyebins = 1;%3;%3;
% popresCorr1D.nXbins = 100;%20;%1;%
if strcmp(popresCorr1D.ErrType,'BlankErr') || strcmp(popresCorr1D.ErrType,'BlankPos')
    popresCorr1D.nXbins = 1;%
else
    popresCorr1D.nXbins = 1;%100;
end

popresCorr1D.Nshuffle = 100;
lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

contval = [0.1:0.05:0.9];%[0.2 0.3 0.4];%[0.8 0.9];%
outvalcorr = 2;%5;%[0 1 2 3 4];%[0 1 2 3 4 5];%
if strcmp(popresCorr1D.ShiftType,'Time')
    savedfilename_popresCorr = ['D:\DATA\batch\All\Correlations\popresCorr1D_' popresCorr1D.ShiftType '_' popresCorr1D.AveType '_' popresCorr1D.ErrType '_Twin' num2str(popresCorr1D.Tsmthwin) '_Xwin' num2str(popresCorr1D.Xsmthwin) '_spdth' num2str(popresCorr1D.SpeedThreshold)...
                                '_Decwin' num2str(popresCorr1D.Tsmthwin_dec) '_nDecbins' num2str(popresCorr1D.nDecbins) '_' num2str(popresCorr1D.nspeedbins) 'speedbins_' num2str(popresCorr1D.neyebins) 'eyebins_' num2str(popresCorr1D.nphsbins) 'thetabins_' popresCorr1D.cellstr '.mat'];
else
    savedfilename_popresCorr = ['D:\DATA\batch\All\Correlations\popresCorr1D_' popresCorr1D.AveType '_' popresCorr1D.ErrType '_Twin' num2str(popresCorr1D.Tsmthwin) '_Xwin' num2str(popresCorr1D.Xsmthwin) '_spdth' num2str(popresCorr1D.SpeedThreshold)...
                                '_Decwin' num2str(popresCorr1D.Tsmthwin_dec) '_nDecbins' num2str(popresCorr1D.nDecbins) '_' num2str(popresCorr1D.nspeedbins) 'speedbins_' num2str(popresCorr1D.neyebins) 'eyebins_' num2str(popresCorr1D.nphsbins) 'thetabins_' popresCorr1D.cellstr '.mat'];
end
for g = [2 1 3]
    popresCorr1D.all_nDataPoints{g} = 0;
    popresCorr1D.all_CA1V1Cross{g} = 0;
    popresCorr1D.all_CA1V1Cross_shuffled{g} = 0;
    popresCorr1D.all_CA1marginal{g} = 0;
    popresCorr1D.all_V1marginal{g} = 0;
end

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
                speeds(~isnan(EXP.data.es.ballspeed)) = smthInTime(EXP.data.es.ballspeed(~isnan(EXP.data.es.ballspeed)), popresCorr1D.sampleRate, popresCorr1D.SpdSmthWin, 'same', [], 'boxcar_centered');
                eyeX = NaN(size(EXP.data.es.eyeXpos));
                eyeX(~isnan(EXP.data.es.eyeXpos)) = smthInTime(EXP.data.es.eyeXpos(~isnan(EXP.data.es.eyeXpos)), popresCorr1D.sampleRate, popresCorr1D.SpdSmthWin, 'same', [], 'boxcar_centered');
                
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
                Prange = size(EXP.Bayes.PosError0{1},2);
                
                if strcmp(popresCorr1D.ErrType,'Error')
                    PosError0{1} = EXP.Bayes.PosError0{1};%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.PosError0{2};%EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                elseif strcmp(popresCorr1D.ErrType,'Bias')
                    PosError0{1} = zeros(size(EXP.Bayes.PosError0{1}));%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = zeros(size(EXP.Bayes.PosError0{2}));%EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                    for tt = 1:size(PosError0{1},1)
                        PosError0{1}(tt,:) = circshift(EXP.Bayes.Posterior0{1}(tt,:),floor(Prange/2)-EXP.Bayes.MaxDecodedPosition0{1}(tt),2);
                        PosError0{2}(tt,:) = circshift(EXP.Bayes.Posterior0{2}(tt,:),floor(Prange/2)-EXP.Bayes.MaxDecodedPosition0{2}(tt),2);
                    end
                elseif strcmp(popresCorr1D.ErrType,'DecPos')
                    PosError0{1} = EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                elseif strcmp(popresCorr1D.ErrType,'BlankErr')
                    PosError0{1} = EXP.Bayes.PosError0{1};%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.PosError0{2};%EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                elseif strcmp(popresCorr1D.ErrType,'BlankPos')
                    PosError0{1} = EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                end
                Prange = size(PosError0{1},2);

                tvec = 1:size(PosError0{2},1);
                for xx = 1:size(PosError0{2},2)
                    PosError0{2}(:,xx) = interp1(tvec(~isnan(PosError0{2}(:,xx)) & sum(PosError0{2},2)>0),PosError0{2}(~isnan(PosError0{2}(:,xx)) & sum(PosError0{2},2)>0,xx),1:size(PosError0{2}(:,xx),1));
                end
                for xx = 1:size(PosError0{1},2)
                    PosError0{1}(:,xx) = interp1(tvec(~isnan(PosError0{1}(:,xx)) & sum(PosError0{1},2)>0),PosError0{1}(~isnan(PosError0{1}(:,xx)) & sum(PosError0{1},2)>0,xx),1:size(PosError0{1}(:,xx),1));
                end
%                 fq = [6 9];
%                 for xdec = 1:size(PosError0{1},2)
%                     vec = PosError0{1}(:,xdec);
%                     vec = LFPfilter(vec, fq(1), fq(2), 60);
%                     PosError0{1}(:,xdec) = vec;%
%                     vec = PosError0{2}(:,xdec);
%                     vec = LFPfilter(vec, fq(1), fq(2), 60);
%                     PosError0{2}(:,xdec) = vec;%
%                 end
                
                goodtidx = cell(1,4);
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
                        if strcmp(popresCorr1D.ErrType,'BlankErr') || strcmp(popresCorr1D.ErrType,'BlankPos')
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
%                         goodtidx{g} = sum(PosError0{1}(:,floor(Prange/4):floor(3*Prange/4)),2)>floor(Prange/2) &...
%                                       sum(PosError0{2}(:,floor(Prange/4):floor(3*Prange/4)),2)>floor(Prange/2);%
%                         if g==2
%                             disp(sum(tidx & goodtidx{g})/sum(tidx));
%                         end
%                         tidx = tidx & goodtidx{g};
                        
                        if ~strcmp(popresCorr1D.ErrType,'BlankErr') && ~strcmp(popresCorr1D.ErrType,'BlankPos')
                            Xbinned = EXP.Bayes.X;%floor(((EXP.Bayes.X-1)/Xrange)*popresCorr1D.nXbins)+1;%EXP.Bayes.X;%
                        else
                            Xbinned = ones(size(EXP.Bayes.X));
                        end
                        Xrange = round(max(Xbinned));
                        
                        idxref = tidx;
                        
                        itraj = Xbinned;
                        ntrajbins = max(itraj);
                        spdquantilelim = zeros(ntrajbins,2);
                        Speedbinned{g} = NaN(size(EXP.data.es.ballspeed));
                        if popresCorr1D.nSpdbins > 1
                            for spd = 1:popresCorr1D.nSpdbins
                                for xx = 1:ntrajbins
                                    spdquantilelim(xx,1) = quantile(speeds(idxref & itraj == xx),max(0,(spd-1)/popresCorr1D.nSpdbins));
                                    spdquantilelim(xx,2) = quantile(speeds(idxref & itraj == xx),min(1,(spd)/popresCorr1D.nSpdbins));
                                end
                                Speedbinned{g}(speeds >= spdquantilelim(itraj,1) & speeds < spdquantilelim(itraj,2)) = spd;
                                if spd == 1
                                    Speedbinned{g}(speeds <= spdquantilelim(itraj,1)) = spd;
                                end
                                if spd == popresCorr1D.nSpdbins
                                    Speedbinned{g}(speeds >= spdquantilelim(itraj,2)) = spd;
                                end
                            end
                        else
                            Speedbinned{g} = ones(size(EXP.data.es.ballspeed));
                        end
                        
                        eyequantilelim = zeros(ntrajbins,2);
                        Eyebinned{g} = NaN(size(EXP.data.es.eyeXpos));
                        if popresCorr1D.nEyebins > 1
                            for ieye = 1:popresCorr1D.nEyebins
                                for xx = 1:ntrajbins
                                    eyequantilelim(xx,1) = quantile(eyeX(idxref & itraj == xx),max(0,(ieye-1)/popresCorr1D.nEyebins));
                                    eyequantilelim(xx,2) = quantile(eyeX(idxref & itraj == xx),min(1,(ieye)/popresCorr1D.nEyebins));
                                end
                                Eyebinned{g}(eyeX >= eyequantilelim(itraj,1) & eyeX < eyequantilelim(itraj,2)) = ieye;
                                if ieye == 1
                                    Eyebinned{g}(eyeX <= eyequantilelim(itraj,1)) = ieye;
                                end
                                if ieye == popresCorr1D.nEyebins
                                    Eyebinned{g}(eyeX >= eyequantilelim(itraj,2)) = ieye;
                                end
                            end
                        else
                            Eyebinned{g} = ones(size(EXP.data.es.eyeXpos));
                        end
                        %
                        
                        %here we remove the average posterior for every
                        %position x spd x eye bin. 
                        %Quick note: one alternative would be
                        %to project out the average: this could account for
                        %coherent fluctuations which might result simply
                        %from change in posterior height, ie in overall
                        %spiking rate. But the problem in that case is that
                        %if the height of the post changes, the width also
                        %changes and projecting out therefore doesn't help
                        %much.
                        if ismember(g,[2 1 3])
                            for ieye = 1:popresCorr1D.nEyebins
                                for spd = 1:popresCorr1D.nSpdbins
                                    for xx = 1:Xrange
                                        tidx_corr = tidx & Xbinned == xx & Speedbinned{g} == spd & Eyebinned{g} == ieye & ~isnan(sum(PosError0{1},2))  & ~isnan(sum(PosError0{2},2));
                                        idx_corrCA1 = find(tidx_corr);%find(tidx & Speedbinned{g} == spd);%
                                        idx_corrV1 = find(tidx_corr);%find(tidx & Speedbinned{g} == spd);%
                                        PosError0{1}(idx_corrCA1,:) = PosError0{1}(idx_corrCA1,:) - repmat(nanmean(PosError0{1}(idx_corrCA1,:),1),[numel(idx_corrCA1) 1]);
                                        PosError0{2}(idx_corrV1,:) = PosError0{2}(idx_corrV1,:) - repmat(nanmean(PosError0{2}(idx_corrV1,:),1),[numel(idx_corrV1) 1]);
                                        
%                                         replace the previous 2 lines with the following if you wanna project out the mean 
%                                         instead of just subtracting it
%                                         PosError0{1}(idx_corrCA1,:) = PosError0{1}(idx_corrCA1,:) - repmat(nanmean(PosError0{1}(idx_corrCA1,:),2),[1 Xrange]);
%                                         PosError0{2}(idx_corrV1,:) = PosError0{2}(idx_corrV1,:) - repmat(nanmean(PosError0{2}(idx_corrV1,:),2),[1 Xrange]);
%                                         CA1PosErrorMean = nanmean(PosError0{1}(idx_corrCA1,:),1);
%                                         V1PosErrorMean = nanmean(PosError0{2}(idx_corrV1,:),1); 
%                                         PosError0{1}(idx_corrCA1,:) = PosError0{1}(idx_corrCA1,:) - (PosError0{1}(idx_corrCA1,:)*((CA1PosErrorMean./sum(CA1PosErrorMean.^2))')).*repmat(CA1PosErrorMean,[numel(idx_corrCA1) 1]);
%                                         PosError0{2}(idx_corrV1,:) = PosError0{2}(idx_corrV1,:) - (PosError0{2}(idx_corrV1,:)*((V1PosErrorMean./sum(V1PosErrorMean.^2))')).*repmat(V1PosErrorMean,[numel(idx_corrV1) 1]);
                                        
%                                         PosError0{1}(idx_corrCA1,:) = PosError0{1}(idx_corrCA1,:)/sqrt(sum(sum(PosError0{1}(idx_corrCA1,:).^2,1),2));
%                                         PosError0{2}(idx_corrV1,:) = PosError0{2}(idx_corrV1,:)/sqrt(sum(sum(PosError0{2}(idx_corrV1,:).^2,1),2));
                                    end
                                end
                            end
                        end                        
                    end
                end
%                 fq = [6 9];
%                 for xdec = 1:size(PosError0{1},2)
%                     vec = PosError0{1}(:,xdec);
%                     vec = LFPfilter(vec, fq(1), fq(2), 60);
%                     PosError0{1}(:,xdec) = vec;%
%                     vec = PosError0{2}(:,xdec);
%                     vec = LFPfilter(vec, fq(1), fq(2), 60);
%                     PosError0{2}(:,xdec) = vec;%
%                 end
                    
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
                        if strcmp(popresCorr1D.ErrType,'BlankErr') || strcmp(popresCorr1D.ErrType,'BlankPos')
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
                        
                        pNshuffle = popresCorr1D.Nshuffle;
                        nXbins = 1;%popresCorr1D.nXbins;
                        nSpdbins = 1;%popresCorr1D.nSpdbins;
                        nEyebins = 1;%popresCorr1D.nEyebins;
                        s_nDataPoints = zeros(nXbins,nSpdbins,nEyebins);
                        s_CA1V1Cross = NaN(nXbins,nSpdbins,nEyebins,Prange);
                        s_CA1V1Cross_shuffled = NaN(nXbins,nSpdbins,nEyebins,Prange,pNshuffle);
                        s_CA1marginal = NaN(nXbins,nSpdbins,nEyebins);
                        s_V1marginal = NaN(nXbins,nSpdbins,nEyebins);
                        
                        ishift = 0;
                        for xshift = -floor(Prange/2)+1:floor(Prange/2)
                            ishift = ishift + 1;
                            if strcmp(popresCorr1D.ShiftType,'Space')
                                PosErrorCA1_shift = circshift(PosError0{1},xshift,2);
                            elseif strcmp(popresCorr1D.ShiftType,'Time')
                                PosErrorCA1_shift = circshift(PosError0{1},xshift,1);
                            end
                            for ieye = 1:nEyebins
                                for spd = 1:nSpdbins
                                    for xx = 1:nXbins
                                        tidx_corr = tidx & ~isnan(sum(PosError0{1},2)) & ~isnan(sum(PosError0{2},2));% & (Xbinned >= (xx-1)*100/popresCorr1D.nXbins) & (Xbinned <= xx*100/popresCorr1D.nXbins) & Speedbinned{g} == spd & Eyebinned{g} == ieye;
                                        idx_corrCA1 = find(tidx_corr);
                                        idx_corrV1 = find(tidx_corr);
                                        s_nDataPoints(xx,spd,ieye) = numel(idx_corrCA1);
                                        s_CA1V1Cross(xx,spd,ieye,ishift) = sum(sum(PosError0{2}(idx_corrV1,:).*PosErrorCA1_shift(idx_corrCA1,:),1),2);%(PosError0{2}(idx_corrV1,:) - repmat(mean(PosError0{2}(idx_corrV1,:),1),[numel(idx_corrV1) 1]))'*(PosError0{1}(idx_corrCA1,:)-repmat(mean(PosError0{1}(idx_corrCA1,:),1),[numel(idx_corrCA1) 1]));%
                                        for ishf = 1:pNshuffle
                                            shuffledidxV1 = circshift(idx_corrV1,round(numel(idx_corrV1)/2) + randi(round(numel(idx_corrV1)/2)) - round(numel(idx_corrV1)/4));
                                            s_CA1V1Cross_shuffled(xx,spd,ieye,ishift,ishf) = sum(sum(PosError0{2}(shuffledidxV1,:).*PosErrorCA1_shift(idx_corrCA1,:),1),2);
                                        end
                                        shuffledidxV1 = randperm(numel(idx_corrV1));
                                        s_CA1V1Cross_shuffled(xx,spd,ieye,ishift) = sum(sum(PosError0{2}(idx_corrV1(shuffledidxV1),:).*PosErrorCA1_shift(idx_corrCA1,:),1),2);
                                        if ishift == 1
                                            s_CA1marginal(xx,spd,ieye) = sum(sum(PosError0{1}(idx_corrCA1,:).^2,1),2);
                                            s_V1marginal(xx,spd,ieye) = sum(sum(PosError0{2}(idx_corrV1,:).^2,1),2);
                                        end
                                    end
                                end
                            end
                        end
                        if strcmp(popresCorr1D.AveType,'full')
                            s_CA1V1Cross = squeeze(nansum(nansum(nansum(s_CA1V1Cross,1),2),3));
                            s_CA1V1Cross_shuffled = squeeze(nansum(nansum(nansum(s_CA1V1Cross_shuffled,1),2),3));
                            s_CA1marginal = squeeze(nansum(nansum(nansum(s_CA1marginal,1),2),3));
                            s_V1marginal = squeeze(nansum(nansum(nansum(s_V1marginal,1),2),3));
                        elseif strcmp(popresCorr1D.AveType,'Pos')
                            s_CA1V1Cross = squeeze(nansum(nansum(s_CA1V1Cross,2),3));
                            s_CA1V1Cross_shuffled = squeeze(nansum(nansum(s_CA1V1Cross_shuffled,2),3));
                            s_CA1marginal = squeeze(nansum(nansum(s_CA1marginal,2),3));
                            s_V1marginal = squeeze(nansum(nansum(s_V1marginal,2),3));
                        end
                        
                        popresCorr1D.all_nDataPoints{g} = popresCorr1D.all_nDataPoints{g} + s_nDataPoints;
                        popresCorr1D.all_CA1V1Cross{g} = popresCorr1D.all_CA1V1Cross{g} + s_CA1V1Cross;
                        popresCorr1D.all_CA1V1Cross_shuffled{g} = popresCorr1D.all_CA1V1Cross_shuffled{g} + s_CA1V1Cross_shuffled(:,1);
                        popresCorr1D.all_CA1marginal{g} = popresCorr1D.all_CA1marginal{g} + s_CA1marginal;
                        popresCorr1D.all_V1marginal{g} = popresCorr1D.all_V1marginal{g} + s_V1marginal;
                        
                        if strcmp(popresCorr1D.AveType,'full')
                            popresCorr1D.s_CA1V1Cross{ianimal,iseries,g} = s_CA1V1Cross./sqrt(s_V1marginal*s_CA1marginal');
                            popresCorr1D.s_CA1V1Cross_shuffled{ianimal,iseries,g} = s_CA1V1Cross_shuffled./sqrt(s_V1marginal*s_CA1marginal');
                        elseif strcmp(popresCorr1D.AveType,'Pos')
                            popresCorr1D.s_CA1V1Cross{ianimal,iseries,g} = NaN(nXbins,Prange);
                            popresCorr1D.s_CA1V1Cross_shuffled{ianimal,iseries,g} = NaN(nXbins,Prange);
                            for xx = 1:nXbins
                                popresCorr1D.s_CA1V1Cross{ianimal,iseries,g}(xx,:) = s_CA1V1Cross(xx,:)./sqrt(s_V1marginal(xx,:)*s_CA1marginal(xx,:)');
                                popresCorr1D.s_CA1V1Cross_shuffled{ianimal,iseries,g}(xx,:) = s_CA1V1Cross_shuffled(xx,:)./sqrt(s_V1marginal(xx,:)*s_CA1marginal(xx,:)');
                            end
                        end
                        popresCorr1D.s_CA1V1CrossPval{ianimal,iseries,g} = sum(repmat(abs(popresCorr1D.s_CA1V1Cross{ianimal,iseries,g}),[1 popresCorr1D.Nshuffle]) < abs(popresCorr1D.s_CA1V1Cross_shuffled{ianimal,iseries,g}), 2);
                    end
                end
            end
        end
    end
end

%now we do the grand averages by averaging across sessions or across trials
if strcmp(popresCorr1D.AveType,'full')
    catdim = 2;
elseif strcmp(popresCorr1D.AveType,'Pos')
    catdim = 3;
end
for g = [2 1 3]
    %averaging acrross sessions
    popresCorr1D.CA1V1Cross{g} = [];
    popresCorr1D.CA1V1Cross_shuffled{g} = [];
    for ianimal = 1:size(popresCorr1D.s_CA1V1Cross,1)
        for iseries = 1:size(popresCorr1D.s_CA1V1Cross,2)
            if ~isempty(popresCorr1D.s_CA1V1Cross{ianimal,iseries,g})
                if sum(isnan(popresCorr1D.s_CA1V1Cross{ianimal,iseries,g}(:))) ~= numel(popresCorr1D.s_CA1V1Cross{ianimal,iseries,g})
                    popresCorr1D.CA1V1Cross{g} = cat(catdim,popresCorr1D.CA1V1Cross{g},popresCorr1D.s_CA1V1Cross{ianimal,iseries,g});
                    popresCorr1D.CA1V1Cross_shuffled{g} = cat(catdim,popresCorr1D.CA1V1Cross_shuffled{g},popresCorr1D.s_CA1V1Cross_shuffled{ianimal,iseries,g});
                end
            end
        end
    end
    
    CA1V1CrossAve = nanmean(popresCorr1D.CA1V1Cross{g},catdim);
    CA1V1CrossAve_shuffled = nanmean(popresCorr1D.CA1V1Cross_shuffled{g},catdim);
    CA1V1CrossAve_std = 0;
    CA1V1CrossAve_shuffled_std = 0;
    kfold = size(popresCorr1D.CA1V1Cross{g},catdim);
    allsessions = 1:size(popresCorr1D.CA1V1Cross{g},catdim);
    for k = 1:kfold
        CA1V1CrossAve_iter = nanmean(popresCorr1D.CA1V1Cross{g}(:,~ismember(allsessions,k)),2);
        CA1V1CrossAve_shuffled_iter = nanmean(popresCorr1D.CA1V1Cross_shuffled{g}(:,~ismember(allsessions,k)),2);
        CA1V1CrossAve_std = CA1V1CrossAve_std + (kfold - 1)/kfold*(CA1V1CrossAve_iter - CA1V1CrossAve).^2;
        CA1V1CrossAve_shuffled_std = CA1V1CrossAve_shuffled_std + (kfold - 1)/kfold*(CA1V1CrossAve_shuffled_iter - CA1V1CrossAve_shuffled).^2;
    end
    popresCorr1D.CA1V1Cross_SE{g} = sqrt(CA1V1CrossAve_std);
    popresCorr1D.CA1V1Cross_shuffled_SE{g} = sqrt(CA1V1CrossAve_shuffled_std);
    popresCorr1D.CA1V1Cross{g} = CA1V1CrossAve;
    popresCorr1D.CA1V1Cross_shuffled{g} = CA1V1CrossAve_shuffled;
    
    %averaging acrross trials
    if strcmp(popresCorr1D.AveType,'full')
        popresCorr1D.all_CA1V1Cross{g} = popresCorr1D.all_CA1V1Cross{g}./sqrt(popresCorr1D.all_V1marginal{g}*popresCorr1D.all_CA1marginal{g}');
        popresCorr1D.all_CA1V1Cross_shuffled{g} = popresCorr1D.all_CA1V1Cross_shuffled{g}./sqrt(popresCorr1D.all_V1marginal{g}*popresCorr1D.all_CA1marginal{g}');
    elseif strcmp(popresCorr1D.AveType,'Pos')
        for xx = 1:nXbins
            popresCorr1D.all_CA1V1Cross{g}(xx,:) = popresCorr1D.all_CA1V1Cross{g}(xx,:)./sqrt(popresCorr1D.all_V1marginal{g}(xx,:)*popresCorr1D.all_CA1marginal{g}(xx,:)');
            popresCorr1D.all_CA1V1Cross_shuffled{g}(xx,:) = popresCorr1D.all_CA1V1Cross_shuffled{g}(xx,:)./sqrt(popresCorr1D.all_V1marginal{g}(xx,:)*popresCorr1D.all_CA1marginal{g}(xx,:)');
        end
    end
end
save(savedfilename_popresCorr, 'popresCorr1D','-v7.3');
end

%to plot correlation coefficient per session
% Pvalcorr = [];
% crosscoeff = [];
% crosscorr = [];
% for ianimal = 1:size(popresCorr1D.s_CA1V1CrossPval,1)
% for iseries = 1:size(popresCorr1D.s_CA1V1CrossPval,2)
% if ~isempty(popresCorr1D.s_CA1V1CrossPval{ianimal,iseries,2})
% [~,imax] = max(popresCorr1D.s_CA1V1Cross{ianimal,iseries,2}(30:70));
% imax = imax+30-1;
% Pvalcorr = [Pvalcorr popresCorr1D.s_CA1V1CrossPval{ianimal,iseries,2}(imax)/100];
% crosscoeff = [crosscoeff popresCorr1D.s_CA1V1Cross{ianimal,iseries,2}(imax)];
% crosscorr = [crosscorr; popresCorr1D.s_CA1V1Cross{ianimal,iseries,2}(:)'];
% end
% end
% end
% Pvalcorr_blank = [];
% crosscoeff_blank = [];
% crosscorr_blank = [];
% for ianimal = 1:size(popresCorr1D_blank.s_CA1V1CrossPval,1)
% for iseries = 1:size(popresCorr1D_blank.s_CA1V1CrossPval,2)
% if ~isempty(popresCorr1D_blank.s_CA1V1CrossPval{ianimal,iseries,2})
% [~,imax] = max(popresCorr1D_blank.s_CA1V1Cross{ianimal,iseries,2}(30:70));
% imax = imax+30-1;
% Pvalcorr_blank = [Pvalcorr_blank popresCorr1D_blank.s_CA1V1CrossPval{ianimal,iseries,2}(imax)/100];
% crosscoeff_blank = [crosscoeff_blank popresCorr1D_blank.s_CA1V1Cross{ianimal,iseries,2}(imax)];
% crosscorr_blank = [crosscorr_blank; popresCorr1D_blank.s_CA1V1Cross{ianimal,iseries,2}(:)'];
% end
% end
% end
% figure;scatter(crosscoeff(Pvalcorr_blank>0.05),crosscoeff_blank(Pvalcorr_blank>0.05))
% hold on;scatter(crosscoeff(Pvalcorr_blank<=0.05),crosscoeff_blank(Pvalcorr_blank<=0.05))



%use the following to visualize the data until there's a GUI or fig
%associated
% for g = [2 1 3]
% for tlag = 1:1
% count = 0;
% CovMat_Ave{tlag,g} = 0;
% CovMat_Shuffled{tlag,g} = 0;
% for ianimal = 1:size(popresCorr1D.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr1D.s_CA1V1Cov,2)
% if ~isempty(popresCorr1D.s_CA1V1Cov{ianimal,iseries,tlag,g})
% %This will average across sessions than across all trials
% mat = popresCorr1D.s_CA1V1Joint{ianimal,iseries,tlag,g}./repmat(popresCorr1D.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% if sum(isnan(mat)) == 0
% CovMat_Ave{tlag,g} = CovMat_Ave{tlag,g} + mat;
% CovMat_Shuffled{tlag,g} = CovMat_Shuffled{tlag,g} + popresCorr1D.s_CA1V1Joint_shuffled{ianimal,iseries,tlag,g}./repmat(popresCorr1D.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% count = count+ 1;%popresCorr1D.s_nDataPoints{ianimal,iseries,tlag,g};
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
% for ianimal = 1:size(popresCorr1D.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr1D.s_CA1V1Cov,2)
% if ~isempty(popresCorr1D.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresCorr1D.s_CA1V1Joint{ianimal,iseries,g};
% matshf = popresCorr1D.s_CA1V1Joint_shuffled{ianimal,iseries,g};
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% matcross(xx,spd,:,:) = squeeze(mat(xx,spd,:,:))./(squeeze(popresCorr1D.s_V1marginal{ianimal,iseries,g}(xx,spd,:))*squeeze(popresCorr1D.s_CA1marginal{ianimal,iseries,g}(xx,spd,:))').^0.5;
% matcross_shf(xx,spd,:,:) = squeeze(matshf(xx,spd,:,:))./(squeeze(popresCorr1D.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))*squeeze(popresCorr1D.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))').^0.5;
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
% for ianimal = 1:size(popresCorr1D.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr1D.s_CA1V1Cov,2)
% if ~isempty(popresCorr1D.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresCorr1D.s_CA1V1Joint{ianimal,iseries,g}./repmat(popresCorr1D.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matshf = popresCorr1D.s_CA1V1Joint_shuffled{ianimal,iseries,g}./repmat(popresCorr1D.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% for ieye = 1:size(mat,3)
% matcross(xx,spd,ieye,:,:) = squeeze(popresCorr1D.s_CA1V1Joint{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresCorr1D.s_V1marginal{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresCorr1D.s_CA1marginal{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% matcross_shf(xx,spd,ieye,:,:) = squeeze(popresCorr1D.s_CA1V1Joint_shuffled{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresCorr1D.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresCorr1D.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% end
% end
% end
% if sum(isnan(mat)) == 0
% CovMat_Ave{g} = CovMat_Ave{g} + mat;
% CrossMat_Ave{g} = CrossMat_Ave{g} + matcross;
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g} + matcross_shf;
% CovLagMat_Ave{g} = CovLagMat_Ave{g} + popresCorr1D.s_CA1V1Jointdiag{ianimal,iseries,g}./repmat(popresCorr1D.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
% crosscount = crosscount + 1;
% CovMat_stdCA1{g} = CovMat_stdCA1{g} + popresCorr1D.s_CA1marginal{ianimal,iseries,g};
% CovMat_stdV1{g} = CovMat_stdV1{g} + popresCorr1D.s_V1marginal{ianimal,iseries,g};
% CovMat_Shuffled{g} = CovMat_Shuffled{g} + matshf;
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g} + popresCorr1D.s_CA1V1Jointdiag_shuffled{ianimal,iseries,g}./repmat(popresCorr1D.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
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


%to plot crosscoeff corridor vs blank session per session with significance
% Pvalcorr_blank = [];
% crosscoeff_blank = [];
% crosscorr_blank = [];
% for ianimal = 1:size(popresCorr1D_blank.s_CA1V1CrossPval,1)
% for iseries = 1:size(popresCorr1D_blank.s_CA1V1CrossPval,2)
% if ~isempty(popresCorr1D_blank.s_CA1V1CrossPval{ianimal,iseries,2})
% [~,imax] = max(popresCorr1D_blank.s_CA1V1Cross{ianimal,iseries,2}(30:70));
% imax = imax+30-1;
% Pvalcorr_blank = [Pvalcorr_blank popresCorr1D_blank.s_CA1V1CrossPval{ianimal,iseries,2}(imax)];
% crosscoeff_blank = [crosscoeff_blank popresCorr1D_blank.s_CA1V1Cross{ianimal,iseries,2}(imax)];
% crosscorr_blank = [crosscorr_blank; popresCorr1D_blank.s_CA1V1Cross{ianimal,iseries,2}(:)'];
% end
% end
% end
%
% Pvalcorr = [];
% crosscoeff = [];
% crosscorr = [];
% for ianimal = 1:size(popresCorr1D.s_CA1V1CrossPval,1)
% for iseries = 1:size(popresCorr1D.s_CA1V1CrossPval,2)
% if ~isempty(popresCorr1D.s_CA1V1CrossPval{ianimal,iseries,2})
% [~,imax] = max(popresCorr1D.s_CA1V1Cross{ianimal,iseries,2}(30:70));
% imax = imax+30-1;
% Pvalcorr = [Pvalcorr popresCorr1D.s_CA1V1CrossPval{ianimal,iseries,2}(imax)];
% crosscoeff = [crosscoeff popresCorr1D.s_CA1V1Cross{ianimal,iseries,2}(imax)];
% crosscorr = [crosscorr; popresCorr1D.s_CA1V1Cross{ianimal,iseries,2}(:)'];
% end
% end
% end
% [rho,pval] = corr(crosscoeff',crosscoeff_blank')