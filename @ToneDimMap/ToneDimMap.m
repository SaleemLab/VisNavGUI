classdef ToneDimMap < TspikeMap
    %     Create a spikeMap class object of type '1D'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     [obj, prediction, X] = trainSpikeMap(obj, X, Y, Xsmth_win);
    %     
    % Aman Saleem
    % Jan 2014
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        bins;        % the values of the bin limits
        numBins;     % number of bins by which the variable is discretized
        Fcircular;   % if true, circular smoothing
        qthreshold;
        FcomputePos;
        Fdiscarditer;
        Fshuffle;
        Fgoodcells;
        Fshufflehalf;
    end
    
    methods
        function obj = ToneDimMap(varargin)
            
            pnames = {'dimensionality' 'variable' 'variable_range'...
                'Xsmth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'numBins' 'Fcircular' 'qthreshold' 'FcomputePos' 'Fdiscarditer' 'Fshuffle' 'Fgoodcells' 'Fshufflehalf'};
            dflts  = {'1D' 'P' []...
                1 []...
                20 []...
                [] [] 60 100 true 1 true true false [] false};
            [obj.dimensionality, obj.variable, obj.variable_range, ...
                obj.Xsmth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBins, obj.Fcircular, obj.qthreshold, obj.FcomputePos, obj.Fdiscarditer, obj.Fshuffle, obj.Fgoodcells, obj.Fshufflehalf] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, Prediction, X] = trainSpikeMap(obj, X, Y, T, Xsmth_win, FoptiSmooth)
            if nargin<5
                Xsmth_win = obj.Xsmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
            end
            if nargin < 6
                FoptiSmooth = false;
            end
                        
            if isempty(obj.kfold)
                obj.kfold = 5;
            end
            
            if isempty(obj.train_mean)
                calTrainMean = 1;
            else
                calTrainMean = 0;
            end
            if ~FoptiSmooth
                disp(['fixed spatial window  = ' num2str(Xsmth_win) '%']);
            end
            
            if size(X,1) ~= size(Y,1)
                error('wrong dimension for X in ToneDimMap');
            end
            
            if sum(isnan(X(:)))>0 || sum(isnan(Y(:)))>0
                display('WARNING!!! Nans in the data: making a temp fix');
                t = ones(size(X));
                t(isnan(sum(X,2))) = 0;
                t(isnan(sum(Y,2))) = 0;
                X = X(t>0,:);
                Y = Y(t>0,:);
            end
            
            % if sum(X<0)>0
            %     display('WARNING!!! -ve variable data present: ignorning them');
            %     t = ones(size(X));
            %     t(X<0) = 0;
            %     X = X(t>0);
            %     Y = Y(t>0,:);
            % end
            if isempty(obj.CVO) && ~obj.Fshuffle && ~obj.Fshufflehalf
                obj.CVO = crossValPartition(ones(1,size(X,1)),obj.kfold);
                %     obj.CVO = crossValPartition_MD(obj.CVO, ones(1,length(X)), 1:length(X), obj.kfold);
            end
            
            if obj.kfold == 1 || obj.Fshuffle || obj.Fshufflehalf
                obj.CVO.kfold = 1;
                obj.CVO.train{1}= ones(1,size(X,1));
                obj.CVO.cv{1}   = obj.CVO.train{1};
                obj.CVO.test{1} = obj.CVO.train{1};
            end
            
            obj.performance = zeros(1,obj.kfold);
            Prediction = [];
            
            % Bayes decoder needs discretization
%             [X, obj.bins] = normalise1var(X, obj.numBins);
            
            Ncells = size(Y,2);
            if isempty(obj.Fgoodcells)
                pFgoodcells = true(1,Ncells);
            else
                pFgoodcells = obj.Fgoodcells;
            end
            pnumBins = obj.numBins;
            pbins = obj.bins;
            pCVO = obj.CVO;
            pkfold = obj.kfold;
            psampleRate = obj.sampleRate;
            pXsmth_win = obj.Xsmth_win;
            pFcircular = obj.Fcircular;
            pFcomputePos = obj.FcomputePos;
            pqthreshold = obj.qthreshold;
            pFdiscarditer = obj.Fdiscarditer;
            pFshuffle = obj.Fshuffle;
            pFshufflehalf = obj.Fshufflehalf;
            
            
            for icell = 1:Ncells
                Tuning(icell).respModel = NaN(obj.kfold,obj.numBins);
                
                Tuning(icell).meanrespModel = NaN(1,obj.numBins);
                Tuning(icell).SErespModel = NaN(1,obj.numBins);
                
                if obj.FcomputePos
                    Tuning(icell).respModelXpos = NaN(obj.kfold,1);
                    Tuning(icell).respModelXmax = NaN(obj.kfold,1);
                    Tuning(icell).respModelXAmp = NaN(obj.kfold,1);
                    Tuning(icell).respModelXsinAmp = NaN(obj.kfold,1);
                    Tuning(icell).respModelXsinOffset = NaN(obj.kfold,1);
                    
                    Tuning(icell).respModelXhalfdiff = NaN(obj.kfold,1);
                    Tuning(icell).SErespModelXhalfdiff = NaN;
                    
                    Tuning(icell).meanrespModelXpos = NaN;
                    Tuning(icell).meanrespModelXmax = NaN;
                    Tuning(icell).meanrespModelXAmp = NaN;
                    Tuning(icell).meanrespModelXsinAmp = NaN;
                    Tuning(icell).meanrespModelXsinOffset = NaN;
                    Tuning(icell).meanrespModelXhalfdiff = NaN;
                    Tuning(icell).SErespModelXpos = NaN;
                    Tuning(icell).SErespModelXmax = NaN;
                    Tuning(icell).SErespModelXAmp = NaN;
                    Tuning(icell).SErespModelXsinAmp = NaN;
                    Tuning(icell).SErespModelXsinOffset = NaN;
                    
                    Tuning(icell).fieldidx = [];
                    Tuning(icell).fieldsize = NaN;
                    Tuning(icell).SEfieldsize = NaN;
                end
                
                if pkfold > 2 && pFdiscarditer
                    if pFcomputePos
                        pTuning(icell).respModelXpos = [];
                        pTuning(icell).respModelXmax = [];
                    end
                end
            end
            
            obj.model.EV = NaN(obj.kfold,Ncells);
            obj.model.L = NaN(obj.kfold,Ncells);
            obj.model.Q = NaN(obj.kfold,Ncells);
            obj.model.skaggs = NaN(obj.kfold,Ncells);
            obj.model.train_mean = NaN(obj.kfold,Ncells);
            obj.model.swin = NaN(obj.kfold,Ncells);
            
            pTuning = Tuning(pFgoodcells);
            goodcellsidx = find(pFgoodcells);
            Ngoodcells = numel(goodcellsidx);
            pEV = NaN(obj.kfold,Ngoodcells);
            pL = NaN(obj.kfold,Ngoodcells);
            pQ = NaN(obj.kfold,Ngoodcells);
            pskaggs = NaN(obj.kfold,Ngoodcells);
            ptrain_mean = NaN(obj.kfold,Ngoodcells);
            pswin = NaN(obj.kfold,Ngoodcells);
            
            if pFshuffle || pFshufflehalf
                spmd
                    rng(0,'combRecursive');
                end
            end
            parfor icell = 1:Ngoodcells
                for iter = 1:pkfold
                    Xcell = X(:,min(end,goodcellsidx(icell)));
                    if pFshuffle
                        Xcell = circshift(Xcell,round(numel(Xcell)/2) + randi(round(numel(Xcell)/2)) - round(numel(Xcell)/4));%Xcell(randperm(numel(Xcell)));
                        kiter = 1;
                    elseif pFshufflehalf
                        half_1 = Xcell <= round(pnumBins/2);
                        half_2 = Xcell > round(pnumBins/2);
                        shift_idx = round(numel(Xcell)/2) + randi(round(numel(Xcell)/2)) - round(numel(Xcell)/4);
                        half_1 = circshift(half_1, shift_idx);
                        half_2 = circshift(half_2, shift_idx);
                        Xcell(half_1) = mod(Xcell(half_1)-1,round(pnumBins/2)) + 1;
                        Xcell(half_2) = mod(Xcell(half_2)-1,round(pnumBins/2)) + 1 + round(pnumBins/2);
%                         for xx = 1:pnumBins
%                             Xidx = find(Xcell == xx | Xcell == xx + round(pnumBins/2) | Xcell == xx - round(pnumBins/2));
%                             Xcell(Xidx) = Xcell(Xidx(randperm(numel(Xidx))));
%                         end
                        kiter = 1;
                    else
                        kiter = iter;
                    end
                    [model, ~] = get1Dmap(Y(:,goodcellsidx(icell)), Xcell', 1./T, pnumBins, pbins, pCVO, kiter, psampleRate, pXsmth_win, FoptiSmooth, pFcircular);
                    pTuning(icell).respModel(iter,:) = model.tuning;
                    pEV(iter,icell) = model.EV;
                    pL(iter,icell) = model.L;
                    pQ(iter,icell) = model.Q;
                    pskaggs(iter,icell) = model.skaggs;
                    ptrain_mean(iter,icell) = model.train_mean;
                    pswin(iter,icell) = model.swin;
                end
            end
            
            Xmesh = 1:obj.numBins;
            parfor icell = 1:Ngoodcells
                if ~pFshuffle && ~pFshufflehalf
                    pTuning(icell).meanrespModel = squeeze(mean(pTuning(icell).respModel,1));
                else
                    pTuning(icell).meanrespModel = squeeze(mean(pTuning(icell).respModel(1,:),1));
                end
                
                if pFcomputePos
                    [~,fieldX] = findfield(pTuning(icell).meanrespModel,pqthreshold);
                    pTuning(icell).fieldidx = fieldX;
                    pTuning(icell).fieldsize = numel(fieldX);
                    outfieldX = find(~ismember(1:pnumBins,fieldX));
                    
                    tuningX0 = squeeze(pTuning(icell).meanrespModel);
                    tuningX = tuningX0;
                    if ~isempty(fieldX)
                        tuningX(outfieldX) = 0;%min(tuningX(fieldX));
                    end
                    if pFcircular
                        pTuning(icell).meanrespModelXpos = getCircularAverage(tuningX',0,1);%mod(circ_mean(pTuning(icell).respModelXpos/pnumBins*2*pi,[],1)+2*pi,2*pi)*pnumBins/(2*pi);
                        pTuning(icell).meanrespModelXmax = getCircularAverage(tuningX0',0,0.1,0.05);%mod(circ_mean(pTuning(icell).respModelXmax/pnumBins*2*pi,[],1)+2*pi,2*pi)*pnumBins/(2*pi);
                        [maxtuningX0, imaxtuningX0] = max(tuningX0);
                        [~, imintuningX0] = min(tuningX0);
                        pTuning(icell).meanrespModelXAmp = tuningX0(max(imaxtuningX0,imintuningX0)) - tuningX0(min(imaxtuningX0,imintuningX0));
                        [pTuning(icell).meanrespModelXsinOffset,pTuning(icell).meanrespModelXsinAmp] = getThetaPrecPhase(tuningX0);
                        
                        imaxtuningX1 = mod(((imaxtuningX0-round(pnumBins/4)):(imaxtuningX0+round(pnumBins/4))) + round(pnumBins/2) - 1, pnumBins)+1;
                        [maxtuningX1, ~] = max(tuningX0(imaxtuningX1));
                        if imaxtuningX0 > round(pnumBins/2)
                            pTuning(icell).meanrespModelXhalfdiff = maxtuningX0 - maxtuningX1;
                        else
                            pTuning(icell).meanrespModelXhalfdiff = maxtuningX1 - maxtuningX0;
                        end
%                         pTuning(icell).meanrespModelXhalfdiff = tuningX0(max(imaxtuningX0,imaxtuningX1)) - tuningX0(min(imaxtuningX0,imaxtuningX1));
                    else
                        pTuning(icell).meanrespModelXpos = sum(tuningX.*Xmesh)./sum(tuningX);%squeeze(nanmean(pTuning(icell).respModelXpos,1));
                        [~, pTuning(icell).meanrespModelXmax] = max(tuningX0);%squeeze(nanmean(pTuning(icell).respModelXmax,1));
                        [maxtuningX0, imaxtuningX0] = max(tuningX0);
                        [~, imintuningX0] = min(tuningX0);
                        pTuning(icell).meanrespModelXAmp = tuningX0(max(imaxtuningX0,imintuningX0)) - tuningX0(min(imaxtuningX0,imintuningX0));
                        pTuning(icell).meanrespModelXsinOffset = 0;
                        pTuning(icell).meanrespModelXsinAmp = 0;
                        
                        imaxtuningX1 = mod(((imaxtuningX0-round(pnumBins/4)):(imaxtuningX0+round(pnumBins/4))) + round(pnumBins/2) - 1, pnumBins)+1;
                        [maxtuningX1, ~] = max(tuningX0(imaxtuningX1));
                        if imaxtuningX0 > round(pnumBins/2)
                            pTuning(icell).meanrespModelXhalfdiff = maxtuningX0 - maxtuningX1;
                        else
                            pTuning(icell).meanrespModelXhalfdiff = maxtuningX1 - maxtuningX0;
                        end                        
%                         imaxtuningX1 = mod(imaxtuningX0 + round(pnumBins/2) - 1, pnumBins)+1;
%                         pTuning(icell).meanrespModelXhalfdiff = tuningX0(max(imaxtuningX0,imaxtuningX1)) - tuningX0(min(imaxtuningX0,imaxtuningX1));
                    end
                    
                    for iter = 1:pkfold
                        tuningX0 = squeeze(pTuning(icell).respModel(iter,:));
                        tuningX = tuningX0;
                        if ~isempty(fieldX)
                            tuningX(outfieldX) = 0;%min(tuningX(fieldX));
                        end
                        if pFcircular
                            pTuning(icell).respModelXpos(iter) = getCircularAverage(tuningX',0,1);
                            pTuning(icell).respModelXmax(iter) = getCircularAverage(tuningX0',0,0.1,0.05);
                            [maxtuningX0, imaxtuningX0] = max(tuningX0);
                            [~, imintuningX0] = min(tuningX0);
                            pTuning(icell).respModelXAmp(iter) = tuningX0(max(imaxtuningX0,imintuningX0)) - tuningX0(min(imaxtuningX0,imintuningX0));
                            [pTuning(icell).respModelXsinOffset(iter),pTuning(icell).respModelXsinAmp(iter)] = getThetaPrecPhase(tuningX0);%getThetaPrecPhase(tuningX0, pTuning(icell).meanrespModelXsinOffset);
                            
                            imaxtuningX1 = mod(((imaxtuningX0-round(pnumBins/4)):(imaxtuningX0+round(pnumBins/4))) + round(pnumBins/2) - 1, pnumBins)+1;
                            [maxtuningX1, ~] = max(tuningX0(imaxtuningX1));
                            if imaxtuningX0 > round(pnumBins/2)
                                pTuning(icell).respModelXhalfdiff(iter) = maxtuningX0 - maxtuningX1;
                            else
                                pTuning(icell).respModelXhalfdiff(iter) = maxtuningX1 - maxtuningX0;
                            end
%                             imaxtuningX1 = mod(imaxtuningX0 + round(pnumBins/2) - 1, pnumBins)+1;
%                             pTuning(icell).respModelXhalfdiff(iter) = tuningX0(max(imaxtuningX0,imaxtuningX1)) - tuningX0(min(imaxtuningX0,imaxtuningX1));
                        else
                            pTuning(icell).respModelXpos(iter) = sum(tuningX.*Xmesh)./sum(tuningX);
                            [~, pTuning(icell).respModelXmax(iter)] = max(tuningX0);
                            [maxtuningX0, imaxtuningX0] = max(tuningX0);
                            [~, imintuningX0] = min(tuningX0);
                            pTuning(icell).respModelXAmp(iter) = tuningX0(max(imaxtuningX0,imintuningX0)) - tuningX0(min(imaxtuningX0,imintuningX0));
                            pTuning(icell).respModelXsinOffset(iter) = 0;
                            pTuning(icell).respModelXsinAmp(iter) = 0;
                            
                            imaxtuningX1 = mod(((imaxtuningX0-round(pnumBins/4)):(imaxtuningX0+round(pnumBins/4))) + round(pnumBins/2) - 1, pnumBins)+1;
                            [maxtuningX1, ~] = max(tuningX0(imaxtuningX1));
                            if imaxtuningX0 > round(pnumBins/2)
                                pTuning(icell).respModelXhalfdiff(iter) = maxtuningX0 - maxtuningX1;
                            else
                                pTuning(icell).respModelXhalfdiff(iter) = maxtuningX1 - maxtuningX0;
                            end                            
%                             imaxtuningX1 = mod(imaxtuningX0 + round(pnumBins/2) - 1, pnumBins)+1;
%                             pTuning(icell).respModelXhalfdiff(iter) = tuningX0(max(imaxtuningX0,imaxtuningX1)) - tuningX0(min(imaxtuningX0,imaxtuningX1));
                        end
                    end
                end
                
                stdresp = 0;
                stdrespXpos = 0;
                stdrespXmax = 0;
                stdrespXAmp = 0;
                stdrespXsinAmp = 0;
                stdrespXsinOffset = 0;
                stdrespSize = 0;
                stdrespXhalfdiff = 0;
                if pFshuffle || pFshufflehalf
                    StdSEfactor = 1/pkfold;
                else
                    StdSEfactor = (pkfold - 1)/pkfold;
                end
                if ~pFshuffle && ~pFshufflehalf
                    for i = 1:pkfold
                        r_iter = pTuning(icell).respModel(i,:);
                        r_ave = pTuning(icell).meanrespModel;
                        stdresp = stdresp + StdSEfactor*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                        if pFcomputePos
                            [~,fieldX_iter] = findfield(r_iter,pqthreshold);
                            stdrespSize = stdrespSize + StdSEfactor*(numel(fieldX_iter) - pTuning(icell).fieldsize).^2;
                            if pFcircular
                                iterdist = circ_dist((pTuning(icell).respModelXpos(i)/pnumBins*2*pi),pTuning(icell).meanrespModelXpos/pnumBins*2*pi)*pnumBins/(2*pi);
                                stdrespXpos = stdrespXpos + StdSEfactor*iterdist.^2;
                                iterdist = circ_dist((pTuning(icell).respModelXmax(i)/pnumBins*2*pi),pTuning(icell).meanrespModelXmax/pnumBins*2*pi)*pnumBins/(2*pi);
                                stdrespXmax = stdrespXmax + StdSEfactor*iterdist.^2;
                                stdrespXAmp = stdrespXAmp + StdSEfactor*(pTuning(icell).respModelXAmp(i) - pTuning(icell).meanrespModelXAmp).^2;
                                stdrespXsinAmp = stdrespXsinAmp + StdSEfactor*(pTuning(icell).respModelXsinAmp(i) - pTuning(icell).meanrespModelXsinAmp).^2;
                                iterdist = circ_dist(pTuning(icell).respModelXsinOffset(i)/360*2*pi,pTuning(icell).meanrespModelXsinOffset/360*2*pi)*360/(2*pi);
                                stdrespXsinOffset = stdrespXsinOffset + StdSEfactor*iterdist.^2;
                                stdrespXhalfdiff = stdrespXhalfdiff + StdSEfactor*(pTuning(icell).respModelXhalfdiff(i) - pTuning(icell).meanrespModelXhalfdiff).^2;
                            else
                                stdrespXpos = stdrespXpos + StdSEfactor*(pTuning(icell).respModelXpos(i) - pTuning(icell).meanrespModelXpos).^2;
                                stdrespXmax = stdrespXmax + StdSEfactor*(pTuning(icell).respModelXmax(i) - pTuning(icell).meanrespModelXmax).^2;
                                stdrespXAmp = stdrespXAmp + StdSEfactor*(pTuning(icell).respModelXAmp(i) - pTuning(icell).meanrespModelXAmp).^2;
                                stdrespXsinAmp = stdrespXsinAmp + StdSEfactor*(pTuning(icell).respModelXsinAmp(i) - pTuning(icell).meanrespModelXsinAmp).^2;
                                iterdist = circ_dist(pTuning(icell).respModelXsinOffset(i)/360*2*pi,pTuning(icell).meanrespModelXsinOffset/360*2*pi)*360/(2*pi);
                                stdrespXsinOffset = stdrespXsinOffset + StdSEfactor*iterdist.^2;
                                stdrespXhalfdiff = stdrespXhalfdiff + StdSEfactor*(pTuning(icell).respModelXhalfdiff(i) - pTuning(icell).meanrespModelXhalfdiff).^2;
                            end
                        end
                    end
                    
                    pTuning(icell).SErespModel = stdresp.^0.5;
                    if pFcomputePos
                        pTuning(icell).SErespModelXpos = stdrespXpos.^0.5;
                        pTuning(icell).SErespModelXmax = stdrespXmax.^0.5;
                        pTuning(icell).SErespModelXAmp = stdrespXAmp.^0.5;
                        pTuning(icell).SErespModelXsinAmp = stdrespXsinAmp.^0.5;
                        pTuning(icell).SErespModelXsinOffset = stdrespXsinOffset.^0.5;
                        pTuning(icell).SEfieldsize = stdrespSize.^0.5;
                        pTuning(icell).SErespModelXhalfdiff = stdrespXhalfdiff.^0.5;
                    end
                end
                if pkfold > 2 && pFdiscarditer
                    pTuning(icell).respModel = [];
                    if pFcomputePos
                        pTuning(icell).respModelXpos = [];
                        pTuning(icell).respModelXmax = [];
%                         pTuning(icell).respModelXAmp = [];
%                         pTuning(icell).respModelXsinAmp = [];
%                         pTuning(icell).respModelXsinOffset = [];
                    end
                end
            end
            
            Tuning(pFgoodcells) = pTuning;
            obj.model.tuning = Tuning;
            obj.model.EV(:,pFgoodcells) = pEV;
            obj.model.L(:,pFgoodcells) = pL;
            obj.model.Q(:,pFgoodcells) = pQ;
            obj.model.skaggs(:,pFgoodcells) = pskaggs;
            obj.model.train_mean(:,pFgoodcells) = ptrain_mean;
            obj.model.swin(:,pFgoodcells) = pswin;
        end
        
        function [obj, Prediction, X] = trainSpikeMap2(obj, X, Y, T, Xsmth_win, FoptiSmooth)
            if nargin<5
                Xsmth_win = obj.Xsmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
            end
            if nargin < 6
                FoptiSmooth = false;
            end
                        
            if isempty(obj.kfold)
                obj.kfold = 5;
            end
            
            if isempty(obj.train_mean)
                calTrainMean = 1;
            else
                calTrainMean = 0;
            end
            
            if ~FoptiSmooth
                disp(['fixed spatial window  = ' num2str(Xsmth_win) '%']);
            end
            
            if size(X,1) ~= size(Y,1)
                error('wrong dimension for X in ToneDimMap');
            end
            
            if sum(isnan(X))>0 || sum(isnan(Y(:)))>0
                display('WARNING!!! Nans in the data: making a temp fix');
                t = ones(size(X));
                t(isnan(X)) = 0;
                t(isnan(sum(Y,2))) = 0;
                X = X(t>0);
                Y = Y(t>0,:);
            end
            
            % if sum(X<0)>0
            %     display('WARNING!!! -ve variable data present: ignorning them');
            %     t = ones(size(X));
            %     t(X<0) = 0;
            %     X = X(t>0);
            %     Y = Y(t>0,:);
            % end
            if isempty(obj.CVO)
                obj.CVO = crossValPartition(ones(1,size(X,1)),obj.kfold);
                %     obj.CVO = crossValPartition_MD(obj.CVO, ones(1,length(X)), 1:length(X), obj.kfold);
            end
            
            if obj.kfold == 1
                obj.CVO.kfold = 1;
                obj.CVO.train{1}= ones(1,size(X,1));
                obj.CVO.cv{1}   = obj.CVO.train{1};
                obj.CVO.test{1} = obj.CVO.train{1};
            end
            
            obj.performance = zeros(1,obj.CVO.kfold);
            Prediction = [];
            
            % Bayes decoder needs discretization
%             [X, obj.bins] = normalise1var(X, obj.numBins);
            
            pred_temp = cell(obj.CVO.kfold,1);
            
            for iter = 1:obj.CVO.kfold
                
                %     Ytrain  = Y(obj.CVO.train{iter},:);
                %     Ycv     = Y(obj.CVO.cv{iter},:);
                %     Ytest   = Y(obj.CVO.test{iter},:);
                %
                %     Vtrain  = X(obj.CVO.train{iter},:);
                %     Vcv     = X(obj.CVO.cv{iter},:);
                %     Vtest   = X(obj.CVO.test{iter},:);
                
                %     if calTrainMean
                %         train_mean = mean(Vtrain,1);
                %     else
                %         train_mean = obj.train_mean;
                %     end
                
                % Getting the 1D map for each neuron (calculating the place fields)
                % ...the main training component
                %     respModel = zeros(size(Y,2),obj.numBins);
                for icell = 1:size(Y,2)
                    [model, pred] = get1Dmap2(Y(:,icell), X', 1./T, obj.numBins, obj.bins, obj.CVO, iter, obj.sampleRate, obj.Xsmth_win, FoptiSmooth, obj.Fcircular);
                    obj.model.tuning(icell).respModel(iter,:,:) = model.tuning;
                    obj.model.EV(iter,icell) = model.EV;
                    obj.model.L(iter,icell) = model.L;
                    obj.model.Q(iter,icell) = model.Q;
                    obj.model.skaggs(iter,icell) = model.skaggs;
                    obj.model.train_mean(iter,icell) = model.train_mean;
                    obj.model.swin(iter,icell) = model.swin;
                    pred_temp{iter}(icell,:) = pred;                                                            
                end
                obj.model.EV(iter,obj.model.EV(iter,:)<0) = 0;                
            end
            for icell = 1:size(Y,2)
                obj.model.tuning(icell).meanrespModel = squeeze(mean(obj.model.tuning(icell).respModel,1));
                stdresp = 0;
                for i = 1:obj.CVO.kfold
%                     stdresp = stdresp + (obj.CVO.kfold - 1)/obj.CVO.kfold*(obj.model.tuning(icell).respModel(i,:) - obj.model.tuning(icell).meanrespModel).^2;
                end
                obj.model.tuning(icell).SErespModel = stdresp.^0.5;
            end
%             for icell = 1:size(Y,2)
%                 fitobject = fit(obj.bins',mean(obj.model.tuning(icell).respModel, 1)',@(A,k,p,x)(A*exp(k*cos(x/100*2*pi - p)-1)), 'StartPoint', [1 1 0] );
%                 vonmisesfit = feval(fitobject,obj.bins);
%                 [~ , obj.model.position(icell)] = max(vonmisesfit);
%             end
            
            if size(pred_temp{1},1)>1
                for iter = 1:size(pred_temp,1)
                    Prediction = [Prediction pred_temp{iter}];
                end
            else
                for iter = 1:size(pred_temp,1)
                    Prediction = [Prediction pred_temp{iter}];
                end
                Prediction = Prediction';
            end
        end
        
        function [Y, EV] = testMap(obj, X, Yorig)
            minX = min(obj.bins);
            maxX = max(obj.bins);
            X(X<minX) = NaN;
            X(X>maxX) = NaN;
            
            input = normalise1var([minX X' maxX]', obj.numBins);
            input = input(2:end-1);
            
            Y = zeros(length(input),length(obj.model.tuning));
            t = ~isnan(input);
            Y(~t) = NaN;
            for icell = 1:length(obj.model.tuning)
               [~,bestModel] = max(obj.model.EV(:,icell));
               
               tuning = obj.model.tuning(icell).respModel(bestModel,:);
               train  = obj.model.train_mean(bestModel,icell);
               
               Y(t,icell) = tuning(input(t));
               
               EV(icell)  = calCrossValExpVar(train, Yorig(t,icell), Y(t,icell), Yorig(t,icell), Y(t,icell));
            end
        end
    end
end

