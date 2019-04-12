function popresCorr2D = BatchJointCorrelation2(batch2p,ErrType,AveType)
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
    ErrType = 'DecPos';%'Error';%Bias;%'BlankPos';%'BlankErr';%
end
if nargin < 3
    AveType = 'full';%'Pos';%
end
popresCorr2D.ErrType = ErrType;
popresCorr2D.AveType = AveType;
popresCorr2D.Tsmthwin = 15;%250;%250;%150;%300;%40;%120;%50
popresCorr2D.Tsmthwin_dec = 250;%
popresCorr2D.nDecbins = 100;
popresCorr2D.Xsmthwin = 4;%2%1;%
popresCorr2D.SpdSmthWin = popresCorr2D.Tsmthwin_dec;
popresCorr2D.SpeedThreshold = 5;
popresCorr2D.nspeedbins = 5;
popresCorr2D.neyebins = 1;
popresCorr2D.nthetaphsbins = 1;%1;%
popresCorr2D.nphsbins = 1;
popresCorr2D.cellstr = 'goodonly';%'All_50bins';%'goodonly';%'goodonly_unwrapped';%'goodonly';%'All';%
filesuffix_EXP = ['Twin' num2str(popresCorr2D.Tsmthwin) '_' 'Xwin' num2str(popresCorr2D.Xsmthwin) '_' 'spdth' num2str(popresCorr2D.SpeedThreshold) '_' 'Decwin' num2str(popresCorr2D.Tsmthwin_dec) '_' 'nDecbins' num2str(popresCorr2D.nDecbins) '_' num2str(popresCorr2D.nspeedbins) 'speedbins' '_' num2str(popresCorr2D.neyebins) 'eyebins' '_' num2str(popresCorr2D.nthetaphsbins) 'thetabins' '_' popresCorr2D.cellstr];
disp(filesuffix_EXP);
disp(['correlating ' popresCorr2D.ErrType])

popresCorr2D.sampleRate = 60;
popresCorr2D.nSpdbins = 5;%1;%
popresCorr2D.nEyebins = 1;%3;%3;
if strcmp(popresCorr2D.ErrType,'BlankErr') || strcmp(popresCorr2D.ErrType,'BlankPos')
    popresCorr2D.nXbins = 1;%
    gainlist = 2;
else
    popresCorr2D.nXbins = 1;%100;
    gainlist = [2 1 3];
end

lambdaSmooth = 2;
corrmaxlag = 180;
nanimal = numel(expt);

contval = [0.1:0.05:0.9];%[0.2 0.3 0.4];%[0.8 0.9];%
outvalcorr = 2;%5;%[0 1 2 3 4];%[0 1 2 3 4 5];%


savedfilename_popresCorr = ['D:\DATA\batch\All\Correlations\popresCorr2D_' popresCorr2D.AveType '_' popresCorr2D.ErrType '_Twin' num2str(popresCorr2D.Tsmthwin) '_Xwin' num2str(popresCorr2D.Xsmthwin) '_spdth' num2str(popresCorr2D.SpeedThreshold)...
                     '_Decwin' num2str(popresCorr2D.Tsmthwin_dec) '_nDecbins' num2str(popresCorr2D.nDecbins) '_' num2str(popresCorr2D.nspeedbins) 'speedbins_' num2str(popresCorr2D.neyebins) 'eyebins_' num2str(popresCorr2D.nphsbins) 'thetabins_' popresCorr2D.cellstr '.mat'];

for g = [2 1 3]
    popresCorr2D.all_nDataPoints{g} = 0;
    popresCorr2D.all_CA1V1Cross{g} = 0;
    popresCorr2D.all_CA1V1Cross_shuffled{g} = 0;
    popresCorr2D.all_CA1marginal{g} = 0;
    popresCorr2D.all_V1marginal{g} = 0;
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
                S.EXP = [];
                
                cont_list = find(ismember(EXP.SubsetVal.contrast, contval));
                RL_list = find(ismember(EXP.SubsetVal.roomlength, [1]));
                outcome_list = find(ismember(EXP.SubsetVal.outcome, outvalcorr));
                Xrange = max(floor(EXP.Bayes.X));%1;%
                speeds = NaN(size(EXP.data.es.ballspeed));
                speeds(~isnan(EXP.data.es.ballspeed)) = smthInTime(EXP.data.es.ballspeed(~isnan(EXP.data.es.ballspeed)), popresCorr2D.sampleRate, popresCorr2D.SpdSmthWin, 'same', [], 'boxcar_centered');
                eyeX = NaN(size(EXP.data.es.eyeXpos));
                eyeX(~isnan(EXP.data.es.eyeXpos)) = smthInTime(EXP.data.es.eyeXpos(~isnan(EXP.data.es.eyeXpos)), popresCorr2D.sampleRate, popresCorr2D.SpdSmthWin, 'same', [], 'boxcar_centered');
                
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
                
                if strcmp(popresCorr2D.ErrType,'Error')
                    PosError0{1} = EXP.Bayes.PosError0{1};%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.PosError0{2};%EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                elseif strcmp(popresCorr2D.ErrType,'Bias')
                    PosError0{1} = zeros(size(EXP.Bayes.PosError0{1}));%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = zeros(size(EXP.Bayes.PosError0{2}));%EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                    for tt = 1:size(PosError0{1},1)
                        PosError0{1}(tt,:) = circshift(EXP.Bayes.Posterior0{1}(tt,:),floor(Prange/2)-EXP.Bayes.MaxDecodedPosition0{1}(tt),2);
                        PosError0{2}(tt,:) = circshift(EXP.Bayes.Posterior0{2}(tt,:),floor(Prange/2)-EXP.Bayes.MaxDecodedPosition0{2}(tt),2);
                    end
                elseif strcmp(popresCorr2D.ErrType,'DecPos')
                    PosError0{1} = EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                elseif strcmp(popresCorr2D.ErrType,'BlankErr')
                    PosError0{1} = EXP.Bayes.PosError0{1};%EXP.Bayes.Posterior0{1};%getCircularAverage(EXP.Bayes.PosError0{1}',0,1);%EXP.Bayes.PosError0{1}
                    PosError0{2} = EXP.Bayes.PosError0{2};%EXP.Bayes.Posterior0{2};%getCircularAverage(EXP.Bayes.PosError0{2}',0,1);%EXP.Bayes.PosError0{2}
                elseif strcmp(popresCorr2D.ErrType,'BlankPos')
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
                for g = gainlist%2%
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
                        if strcmp(popresCorr2D.ErrType,'BlankErr') || strcmp(popresCorr2D.ErrType,'BlankPos')
                            blanks = ~EXP.Subset.noBlanks;
                            for tshift = -30:30
                                blanks = blanks & circshift(~EXP.Subset.noBlanks,tshift);
                            end
                            for tshift = -15:15
                                blanks = blanks & circshift(~EXP.data.es.lick,tshift);
                            end
                            tidx = blanks & EXP.data.es.smthBallSpd > EXP.Bayes.speed_th & ~isnan(EXP.data.es.smthBallSpd) & ~isnan(EXP.data.es.smthTrajSpd);
%                             for o = 1:numel(outcome_list)
%                                 tidx = tidx & EXP.data.es.outcome == outcome_list(o);
%                             end
%                             tidx = tidx & EXP.Bayes.X == mode(EXP.Bayes.X(tidx));
                        end
                        if ~strcmp(popresCorr2D.ErrType,'BlankErr') && ~strcmp(popresCorr2D.ErrType,'BlankPos')
                            Xbinned = EXP.Bayes.X;%floor(((EXP.Bayes.X-1)/Xrange)*popresCorr2D.nXbins)+1;%EXP.Bayes.X;%
                        else
                            Xbinned = ones(size(EXP.Bayes.X));
                        end
                        Xrange = round(max(Xbinned));
                        
                        idxref = tidx;
                        
                        itraj = Xbinned;
                        ntrajbins = max(itraj);
                        spdquantilelim = zeros(ntrajbins,2);
                        Speedbinned{g} = NaN(size(EXP.data.es.ballspeed));
                        if popresCorr2D.nSpdbins > 1
                            for spd = 1:popresCorr2D.nSpdbins
                                for xx = 1:ntrajbins
                                    spdquantilelim(xx,1) = quantile(speeds(idxref & itraj == xx),max(0,(spd-1)/popresCorr2D.nSpdbins));
                                    spdquantilelim(xx,2) = quantile(speeds(idxref & itraj == xx),min(1,(spd)/popresCorr2D.nSpdbins));
                                end
                                Speedbinned{g}(speeds >= spdquantilelim(itraj,1) & speeds < spdquantilelim(itraj,2)) = spd;
                                if spd == 1
                                    Speedbinned{g}(speeds <= spdquantilelim(itraj,1)) = spd;
                                end
                                if spd == popresCorr2D.nSpdbins
                                    Speedbinned{g}(speeds >= spdquantilelim(itraj,2)) = spd;
                                end
                            end
                        else
                            Speedbinned{g} = ones(size(EXP.data.es.ballspeed));
                        end
                        
                        eyequantilelim = zeros(ntrajbins,2);
                        Eyebinned{g} = NaN(size(EXP.data.es.eyeXpos));
                        if popresCorr2D.nEyebins > 1
                            for ieye = 1:popresCorr2D.nEyebins
                                for xx = 1:ntrajbins
                                    eyequantilelim(xx,1) = quantile(eyeX(idxref & itraj == xx),max(0,(ieye-1)/popresCorr2D.nEyebins));
                                    eyequantilelim(xx,2) = quantile(eyeX(idxref & itraj == xx),min(1,(ieye)/popresCorr2D.nEyebins));
                                end
                                Eyebinned{g}(eyeX >= eyequantilelim(itraj,1) & eyeX < eyequantilelim(itraj,2)) = ieye;
                                if ieye == 1
                                    Eyebinned{g}(eyeX <= eyequantilelim(itraj,1)) = ieye;
                                end
                                if ieye == popresCorr2D.nEyebins
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
                            for ieye = 1:popresCorr2D.nEyebins
                                for spd = 1:popresCorr2D.nSpdbins
                                    for xx = 1:Xrange
                                        tidx_corr = tidx & Xbinned == xx & Speedbinned{g} == spd & Eyebinned{g} == ieye & ~isnan(sum(PosError0{1},2))  & ~isnan(sum(PosError0{2},2));
                                        idx_corrCA1 = find(tidx_corr);%find(tidx & Speedbinned{g} == spd);%
                                        idx_corrV1 = find(tidx_corr);%find(tidx & Speedbinned{g} == spd);%
                                        PosError0{1}(idx_corrCA1,:) = PosError0{1}(idx_corrCA1,:) - repmat(nanmean(PosError0{1}(idx_corrCA1,:),1),[numel(idx_corrCA1) 1]);
                                        PosError0{2}(idx_corrV1,:) = PosError0{2}(idx_corrV1,:) - repmat(nanmean(PosError0{2}(idx_corrV1,:),1),[numel(idx_corrV1) 1]);
                                    end
                                end
                            end
                        end
%                         PosError0{1}(tidx,:) = PosError0{1}(tidx,:)/repmat(sqrt(sum(PosError0{1}(tidx,:).^2,1)),[sum(tidx) 1]);
%                         PosError0{2}(tidx,:) = PosError0{2}(tidx,:)/repmat(sqrt(sum(PosError0{2}(tidx,:).^2,1)),[sum(tidx) 1]);
                    end
                end
                    
                for g = gainlist%2%
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
                        if strcmp(popresCorr2D.ErrType,'BlankErr') || strcmp(popresCorr2D.ErrType,'BlankPos')
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
                                                
                        nXbins = 1;%popresCorr2D.nXbins;
                        nSpdbins = 1;%popresCorr2D.nSpdbins;
                        nEyebins = 1;%popresCorr2D.nEyebins;
                        s_nDataPoints = zeros(nXbins,nSpdbins,nEyebins);
                        s_CA1V1Cross = NaN(nXbins,nSpdbins,nEyebins,Prange,Prange);
                        s_CA1V1Cross_shuffled = NaN(nXbins,nSpdbins,nEyebins,Prange,Prange);
                        s_CA1marginal = NaN(nXbins,nSpdbins,nEyebins,Prange);
                        s_V1marginal = NaN(nXbins,nSpdbins,nEyebins,Prange);
                        for ieye = 1:nEyebins
                            for spd = 1:nSpdbins
                                for xx = 1:nXbins
                                    tidx_corr = tidx & ~isnan(sum(PosError0{1},2)) & ~isnan(sum(PosError0{2},2));% & (Xbinned >= (xx-1)*100/popresCorr2D.nXbins) & (Xbinned <= xx*100/popresCorr2D.nXbins) & Speedbinned{g} == spd & Eyebinned{g} == ieye;
                                    idx_corrCA1 = find(tidx_corr);
                                    idx_corrV1 = find(tidx_corr);
                                    s_nDataPoints(xx,spd,ieye) = numel(idx_corrCA1);
                                    s_CA1V1Cross(xx,spd,ieye,:,:) = PosError0{2}(idx_corrV1,:)'*PosError0{1}(idx_corrCA1,:);
                                    s_CA1marginal(xx,spd,ieye,:) = sum(PosError0{1}(idx_corrCA1,:).^2,1);
                                    s_V1marginal(xx,spd,ieye,:) = sum(PosError0{2}(idx_corrV1,:).^2,1);
                                    shuffledidxV1 = randperm(numel(idx_corrV1));
                                    s_CA1V1Cross_shuffled(xx,spd,ieye,:,:) = PosError0{2}(idx_corrV1(shuffledidxV1),:)'*PosError0{1}(idx_corrCA1,:);
                                end
                            end
                        end
                        if strcmp(popresCorr2D.AveType,'full')
                            s_CA1V1Cross = squeeze(nansum(nansum(nansum(s_CA1V1Cross,1),2),3));
                            s_CA1V1Cross_shuffled = squeeze(nansum(nansum(nansum(s_CA1V1Cross_shuffled,1),2),3));
                            s_CA1marginal = squeeze(nansum(nansum(nansum(s_CA1marginal,1),2),3));
                            s_V1marginal = squeeze(nansum(nansum(nansum(s_V1marginal,1),2),3));
                        elseif strcmp(popresCorr2D.AveType,'Pos')
                            s_CA1V1Cross = squeeze(nansum(nansum(s_CA1V1Cross,2),3));
                            s_CA1V1Cross_shuffled = squeeze(nansum(nansum(s_CA1V1Cross_shuffled,2),3));
                            s_CA1marginal = squeeze(nansum(nansum(s_CA1marginal,2),3));
                            s_V1marginal = squeeze(nansum(nansum(s_V1marginal,2),3));
                        end
                        
                        popresCorr2D.all_nDataPoints{g} = popresCorr2D.all_nDataPoints{g} + s_nDataPoints;
                        popresCorr2D.all_CA1V1Cross{g} = popresCorr2D.all_CA1V1Cross{g} + s_CA1V1Cross;
                        popresCorr2D.all_CA1V1Cross_shuffled{g} = popresCorr2D.all_CA1V1Cross_shuffled{g} + s_CA1V1Cross_shuffled;
                        popresCorr2D.all_CA1marginal{g} = popresCorr2D.all_CA1marginal{g} + s_CA1marginal;
                        popresCorr2D.all_V1marginal{g} = popresCorr2D.all_V1marginal{g} + s_V1marginal;
                        
                        if strcmp(popresCorr2D.AveType,'full')
                            popresCorr2D.s_CA1V1Cross{ianimal,iseries,g} = s_CA1V1Cross./sqrt(s_V1marginal*s_CA1marginal');
                            popresCorr2D.s_CA1V1Cross_shuffled{ianimal,iseries,g} = s_CA1V1Cross_shuffled./sqrt(s_V1marginal*s_CA1marginal');
                        elseif strcmp(popresCorr2D.AveType,'Pos')
                            popresCorr2D.s_CA1V1Cross{ianimal,iseries,g} = NaN(nXbins,Prange,Prange);
                            popresCorr2D.s_CA1V1Cross_shuffled{ianimal,iseries,g} = NaN(nXbins,Prange,Prange);
                            for xx = 1:nXbins
                                popresCorr2D.s_CA1V1Cross{ianimal,iseries,g}(xx,:,:) = squeeze(s_CA1V1Cross(xx,:,:))./sqrt(s_V1marginal(xx,:)'*s_CA1marginal(xx,:));
                                popresCorr2D.s_CA1V1Cross_shuffled{ianimal,iseries,g}(xx,:,:) = squeeze(s_CA1V1Cross_shuffled(xx,:,:))./sqrt(s_V1marginal(xx,:)'*s_CA1marginal(xx,:));
                            end
                        end
                    end
                end
            end
        end
    end
end

%now we do the grand averages by averaging across sessions or across trials
if strcmp(popresCorr2D.AveType,'full')
    catdim = 3;
elseif strcmp(popresCorr2D.AveType,'Pos')
    catdim = 4;
end
for g = gainlist
    %averaging acrross sessions
    popresCorr2D.CA1V1Cross{g} = [];
    popresCorr2D.CA1V1Cross_shuffled{g} = [];
    for ianimal = 1:size(popresCorr2D.s_CA1V1Cross,1)
        for iseries = 1:size(popresCorr2D.s_CA1V1Cross,2)
            if ~isempty(popresCorr2D.s_CA1V1Cross{ianimal,iseries,g})
               if sum(isnan(popresCorr2D.s_CA1V1Cross{ianimal,iseries,g}(:))) ~= numel(popresCorr2D.s_CA1V1Cross{ianimal,iseries,g})
                    popresCorr2D.CA1V1Cross{g} = cat(catdim,popresCorr2D.CA1V1Cross{g},popresCorr2D.s_CA1V1Cross{ianimal,iseries,g});
                    popresCorr2D.CA1V1Cross_shuffled{g} = cat(catdim,popresCorr2D.CA1V1Cross_shuffled{g},popresCorr2D.s_CA1V1Cross_shuffled{ianimal,iseries,g});
                end
            end
        end
    end
    CA1V1CrossAve = nanmean(popresCorr2D.CA1V1Cross{g},catdim);
    CA1V1CrossAve_shuffled = nanmean(popresCorr2D.CA1V1Cross_shuffled{g},catdim);
    CA1V1CrossAve_std = 0;
    CA1V1CrossAve_shuffled_std = 0;
    kfold = size(popresCorr2D.CA1V1Cross{g},catdim);
    allsessions = 1:size(popresCorr2D.CA1V1Cross{g},catdim);
    for k = 1:kfold
        CA1V1CrossAve_iter = nanmean(popresCorr2D.CA1V1Cross{g}(:,~ismember(allsessions,k)),2);
        CA1V1CrossAve_shuffled_iter = nanmean(popresCorr2D.CA1V1Cross_shuffled{g}(:,~ismember(allsessions,k)),2);
        CA1V1CrossAve_std = CA1V1CrossAve_std + (kfold - 1)/kfold*(CA1V1CrossAve_iter - CA1V1CrossAve).^2;
        CA1V1CrossAve_shuffled_std = CA1V1CrossAve_shuffled_std + (kfold - 1)/kfold*(CA1V1CrossAve_shuffled_iter - CA1V1CrossAve_shuffled).^2;
    end
    popresCorr2D.CA1V1Cross_SE{g} = sqrt(CA1V1CrossAve_std);
    popresCorr2D.CA1V1Cross_shuffled_SE{g} = sqrt(CA1V1CrossAve_shuffled_std);
    popresCorr2D.CA1V1Cross{g} = CA1V1CrossAve;
    popresCorr2D.CA1V1Cross_shuffled{g} = CA1V1CrossAve_shuffled;
    
    %averaging acrross trials
    if strcmp(popresCorr2D.AveType,'full')
        popresCorr2D.all_CA1V1Cross{g} = popresCorr2D.all_CA1V1Cross{g}./sqrt(popresCorr2D.all_V1marginal{g}*popresCorr2D.all_CA1marginal{g}');
        popresCorr2D.all_CA1V1Cross_shuffled{g} = popresCorr2D.all_CA1V1Cross_shuffled{g}./sqrt(popresCorr2D.all_V1marginal{g}*popresCorr2D.all_CA1marginal{g}');
    elseif strcmp(popresCorr2D.AveType,'Pos')
        for xx = 1:nXbins
            popresCorr2D.all_CA1V1Cross{g}(xx,:,:) = squeeze(popresCorr2D.all_CA1V1Cross{g}(xx,:,:))./sqrt(popresCorr2D.all_V1marginal{g}(xx,:)'*popresCorr2D.all_CA1marginal{g}(xx,:));
            popresCorr2D.all_CA1V1Cross_shuffled{g}(xx,:,:) = squeeze(popresCorr2D.all_CA1V1Cross_shuffled{g}(xx,:,:))./sqrt(popresCorr2D.all_V1marginal{g}(xx,:)'*popresCorr2D.all_CA1marginal{g}(xx,:));
        end
    end
end
save(savedfilename_popresCorr, 'popresCorr2D','-v7.3');
end



%use the following to visualize the data until there's a GUI or fig
%associated
% for g = [2 1 3]
% for tlag = 1:1
% count = 0;
% CovMat_Ave{tlag,g} = 0;
% CovMat_Shuffled{tlag,g} = 0;
% for ianimal = 1:size(popresCorr2D.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr2D.s_CA1V1Cov,2)
% if ~isempty(popresCorr2D.s_CA1V1Cov{ianimal,iseries,tlag,g})
% %This will average across sessions than across all trials
% mat = popresCorr2D.s_CA1V1Joint{ianimal,iseries,tlag,g}./repmat(popresCorr2D.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% if sum(isnan(mat)) == 0
% CovMat_Ave{tlag,g} = CovMat_Ave{tlag,g} + mat;
% CovMat_Shuffled{tlag,g} = CovMat_Shuffled{tlag,g} + popresCorr2D.s_CA1V1Joint_shuffled{ianimal,iseries,tlag,g}./repmat(popresCorr2D.s_nDataPoints{ianimal,iseries,tlag,g},[1 1 100 100]);
% count = count+ 1;%popresCorr2D.s_nDataPoints{ianimal,iseries,tlag,g};
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
% for ianimal = 1:size(popresCorr2D.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr2D.s_CA1V1Cov,2)
% if ~isempty(popresCorr2D.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresCorr2D.s_CA1V1Joint{ianimal,iseries,g};
% matshf = popresCorr2D.s_CA1V1Joint_shuffled{ianimal,iseries,g};
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% matcross(xx,spd,:,:) = squeeze(mat(xx,spd,:,:))./(squeeze(popresCorr2D.s_V1marginal{ianimal,iseries,g}(xx,spd,:))*squeeze(popresCorr2D.s_CA1marginal{ianimal,iseries,g}(xx,spd,:))').^0.5;
% matcross_shf(xx,spd,:,:) = squeeze(matshf(xx,spd,:,:))./(squeeze(popresCorr2D.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))*squeeze(popresCorr2D.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,:))').^0.5;
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
% for ianimal = 1:size(popresCorr2D.s_CA1V1Cov,1)
% for iseries = 1:size(popresCorr2D.s_CA1V1Cov,2)
% if ~isempty(popresCorr2D.s_CA1V1Cov{ianimal,iseries,g})
% mat = popresCorr2D.s_CA1V1Joint{ianimal,iseries,g}./repmat(popresCorr2D.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matshf = popresCorr2D.s_CA1V1Joint_shuffled{ianimal,iseries,g}./repmat(popresCorr2D.s_nDataPoints{ianimal,iseries,g},[1 1 1 100 100]);
% matcross = zeros(size(mat));
% matcross_shf = zeros(size(mat));
% for xx = 1:size(mat,1)
% for spd = 1:size(mat,2)
% for ieye = 1:size(mat,3)
% matcross(xx,spd,ieye,:,:) = squeeze(popresCorr2D.s_CA1V1Joint{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresCorr2D.s_V1marginal{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresCorr2D.s_CA1marginal{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% matcross_shf(xx,spd,ieye,:,:) = squeeze(popresCorr2D.s_CA1V1Joint_shuffled{ianimal,iseries,g}(xx,spd,ieye,:,:))./(squeeze(popresCorr2D.s_V1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))*squeeze(popresCorr2D.s_CA1marginal_shuffled{ianimal,iseries,g}(xx,spd,ieye,:))').^0.5;
% end
% end
% end
% if sum(isnan(mat)) == 0
% CovMat_Ave{g} = CovMat_Ave{g} + mat;
% CrossMat_Ave{g} = CrossMat_Ave{g} + matcross;
% CrossMat_Shuffled{g} = CrossMat_Shuffled{g} + matcross_shf;
% CovLagMat_Ave{g} = CovLagMat_Ave{g} + popresCorr2D.s_CA1V1Jointdiag{ianimal,iseries,g}./repmat(popresCorr2D.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
% crosscount = crosscount + 1;
% CovMat_stdCA1{g} = CovMat_stdCA1{g} + popresCorr2D.s_CA1marginal{ianimal,iseries,g};
% CovMat_stdV1{g} = CovMat_stdV1{g} + popresCorr2D.s_V1marginal{ianimal,iseries,g};
% CovMat_Shuffled{g} = CovMat_Shuffled{g} + matshf;
% CovLagMat_Shuffled{g} = CovLagMat_Shuffled{g} + popresCorr2D.s_CA1V1Jointdiag_shuffled{ianimal,iseries,g}./repmat(popresCorr2D.s_nDataPoints{ianimal,iseries,g},[1 1 1 1 100]);
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