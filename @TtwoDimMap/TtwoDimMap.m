classdef TtwoDimMap < TspikeMap
    %     Create a spikeMap class object of type '2D'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     [obj, prediction, X] = trainSpikeMap(obj, X, Y, Z, Xsmth_win, Ysmth_win);
    %     
    % Aman Saleem
    % Jan 2014
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        binsX;        % the values of the bin limits
        binsY;        % the values of the bin limits
        numBinsX;     % number of bins by which the variable is discretized
        numBinsY;     % number of bins by which the variable is discretized
        FcircularX;   % if true, circular smoothing
        FcircularY;   % if true, circular smoothing
        qthreshold;
        FcomputePos;
        FcomputeMarg;
        FcomputeCorr;
        Fdiscarditer;
        FshuffleXY;
        tuningRef;
        Fgoodcells;
        refYindex;
    end
    
    methods
        function obj = TtwoDimMap(varargin)
            pnames = {'dimensionality' 'variable' 'variable_range'...
                'Xsmth_win' 'Ysmth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'XnumBins' 'FcircularX' 'YnumBins' 'FcircularY' 'qthreshold' 'FcomputePos' 'FcomputeMarg' 'FcomputeCorr' 'Fdiscarditer' 'tuningRef' 'FshuffleXY' 'Fgoodcells' 'refYindex'};
            dflts  = {'2D' 'PxS' []...
                1 1 []...
                20 []...
                [] [] 60 100 true 18 false 1 true true true true [] false [] []};
            [obj.dimensionality, obj.variable, obj.variable_range, ...
                obj.Xsmth_win, obj.Ysmth_win,obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate, obj.numBinsX, obj.FcircularX, obj.numBinsY, obj.FcircularY, obj.qthreshold,...
                obj.FcomputePos, obj.FcomputeMarg, obj.FcomputeCorr, obj.Fdiscarditer, obj.tuningRef, obj.FshuffleXY, obj.Fgoodcells, obj.refYindex] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, Prediction, X] = trainSpikeMap(obj, X, Y, Z, T, Xsmth_win, Ysmth_win, FoptiSmooth)
            if nargin<6
                Xsmth_win = obj.Xsmth_win;
                Ysmth_win = obj.Ysmth_win;
            else
                obj.Xsmth_win = Xsmth_win;
                obj.Ysmth_win = Ysmth_win;
            end
            if nargin < 8
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
            
            if size(X,1) ~= size(Z,1)
                error('wrong dimension for X in TtwoDimMap');
            end
            if size(Y,1) ~= size(Z,1)
                error('wrong dimension for Y in TtwoDimMap');
            end
            
            if sum(isnan(X(:)))>0 || sum(isnan(Y(:)))>0
                display('WARNING!!! Nans in the data: making a temp fix');
                t = ones(size(X));
                t(isnan(sum(X,2))) = 0;
                t(isnan(sum(Y,2))) = 0;
                X = X(t>0);
                Y = Y(t>0);
                Z = Z(t>0,:);
            end
            
            % if sum(X<0)>0
            %     display('WARNING!!! -ve variable data present: ignorning them');
            %     t = ones(size(X));
            %     t(X<0) = 0;
            %     X = X(t>0);
            %     Y = Y(t>0,:);
            % end
            if isempty(obj.CVO) && ~obj.FshuffleXY
                obj.CVO = crossValPartition(ones(1,size(X,1)),obj.kfold);
                %     obj.CVO = crossValPartition_MD(obj.CVO, ones(1,length(X)), 1:length(X), obj.kfold);
            end
            
            if obj.kfold == 1 || obj.FshuffleXY
                obj.CVO.kfold = 1;
                obj.CVO.train{1}= ones(1,size(X,1));
                obj.CVO.cv{1}   = obj.CVO.train{1};
                obj.CVO.test{1} = obj.CVO.train{1};
            end
            
            obj.performance = zeros(1,obj.CVO.kfold);
            Prediction = [];
            
            % Bayes decoder needs discretization
%             [X, obj.binsX] = normalise1var(X, obj.numBinsX);
%             [Y, obj.binsY] = normalise1var(Y, obj.numBinsY);
            
            Ncells = size(Z,2);
            if isempty(obj.Fgoodcells)
                pFgoodcells = true(1,Ncells);
            else
                pFgoodcells = obj.Fgoodcells;
            end
            pnumBinsX = obj.numBinsX;
            pbinsX = obj.binsX;
            pnumBinsY = obj.numBinsY;
            pbinsY = obj.binsY;
            pCVO = obj.CVO;
            pkfold = obj.kfold;
            psampleRate = obj.sampleRate;
            pXsmth_win = obj.Xsmth_win;
            pYsmth_win = obj.Ysmth_win;
            pFcircularX = obj.FcircularX;
            pFcircularY = obj.FcircularY;
            pFcomputeMarg = obj.FcomputeMarg;
            prefYindex = obj.refYindex;
            pFcomputePos = obj.FcomputePos;
            pFcomputeCorr = obj.FcomputeCorr;
            pqthreshold = obj.qthreshold;
            pFdiscarditer = obj.Fdiscarditer;
            pFshuffleXY = obj.FshuffleXY;
            maxtol = 0.1;
            
            
            
            for icell = 1:Ncells
                Tuning(icell).respModel = NaN(obj.kfold,obj.numBinsY,obj.numBinsX);
                Tuning(icell).meanrespModel = NaN(obj.numBinsY,obj.numBinsX);
                Tuning(icell).SErespModel = NaN(obj.numBinsY,obj.numBinsX);
                pTuning_nosmth(icell).respModel = NaN(obj.kfold,obj.numBinsY,obj.numBinsX);
                pTuning_nosmth(icell).meanrespModel = NaN(obj.numBinsY,obj.numBinsX);                
                if obj.FcomputeMarg
                    Tuning(icell).respModelX = NaN(obj.kfold,obj.numBinsX);
                    Tuning(icell).respModelY = NaN(obj.kfold,obj.numBinsY);
                    
                    Tuning(icell).meanrespModelX = NaN(1,obj.numBinsX);
                    Tuning(icell).meanrespModelY = NaN(obj.numBinsY,1);
                    
                    Tuning(icell).SErespModelX = NaN(obj.numBinsX,1);
                    Tuning(icell).SErespModelY = NaN(obj.numBinsY,1);
                end
                if obj.FcomputeCorr
                    Tuning(icell).corrModelX = NaN(obj.kfold,obj.numBinsY,obj.numBinsX);
                    Tuning(icell).corrModelY = NaN(obj.kfold,obj.numBinsY,obj.numBinsX);
                    Tuning(icell).corrModelXpos = NaN(obj.kfold,obj.numBinsY);
                    Tuning(icell).corrModelXmax = NaN(obj.kfold,obj.numBinsY);
                    Tuning(icell).corrModelXmaxAmp = NaN(obj.kfold,1);
                    Tuning(icell).corrModelXmaxsinAmp = NaN(obj.kfold,1);
                    Tuning(icell).corrModelXmaxOffset = NaN(obj.kfold,1);
                    Tuning(icell).corrModelYpos = NaN(obj.kfold,obj.numBinsX);
                    Tuning(icell).corrModelYmax = NaN(obj.kfold,obj.numBinsX);
                    
                    Tuning(icell).corrModelXRef = NaN(obj.kfold,obj.numBinsY,obj.numBinsX);
                    Tuning(icell).corrModelXRefpos = NaN(obj.kfold,obj.numBinsY);
                    Tuning(icell).corrModelXRefmax = NaN(obj.kfold,obj.numBinsY);
                    Tuning(icell).corrModelXRefmaxAmp = NaN(obj.kfold,1);
                    Tuning(icell).corrModelXRefmaxsinAmp = NaN(obj.kfold,1);
                    Tuning(icell).corrModelXRefmaxOffset = NaN(obj.kfold,1);
                    
                    
                    Tuning(icell).meancorrModelX = NaN(obj.numBinsY,obj.numBinsX);
                    Tuning(icell).meancorrModelY = NaN(obj.numBinsY,obj.numBinsX);
                    Tuning(icell).meancorrModelXpos = NaN(obj.numBinsY,1);
                    Tuning(icell).meancorrModelXmax = NaN(obj.numBinsY,1);
                    Tuning(icell).meancorrModelXmaxAmp = NaN;
                    Tuning(icell).meancorrModelXmaxsinAmp = NaN;
                    Tuning(icell).meancorrModelXmaxOffset = NaN;
                    Tuning(icell).meancorrModelYpos = NaN(obj.numBinsX,1);
                    Tuning(icell).meancorrModelYmax = NaN(obj.numBinsX,1);
                    
                    Tuning(icell).meancorrModelXRef = NaN(obj.numBinsY,obj.numBinsX);
                    Tuning(icell).meancorrModelXRefpos = NaN(obj.numBinsY,1);
                    Tuning(icell).meancorrModelXRefmax = NaN(obj.numBinsY,1);
                    Tuning(icell).meancorrModelXRefmaxAmp = NaN;
                    Tuning(icell).meancorrModelXRefmaxsinAmp = NaN;
                    Tuning(icell).meancorrModelXRefmaxOffset = NaN;
                    
                    Tuning(icell).SEcorrModelX = NaN(obj.numBinsY,obj.numBinsX);
                    Tuning(icell).SEcorrModelY = NaN(obj.numBinsY,obj.numBinsX);
                    Tuning(icell).SEcorrModelXpos = NaN(obj.numBinsY,1);
                    Tuning(icell).SEcorrModelXmax = NaN(obj.numBinsY,1);
                    Tuning(icell).SEcorrModelXmaxAmp = NaN;
                    Tuning(icell).SEcorrModelXmaxsinAmp = NaN;
                    Tuning(icell).SEcorrModelXmaxOffset = NaN;
                    Tuning(icell).SEcorrModelYpos = NaN(obj.numBinsX,1);
                    Tuning(icell).SEcorrModelYmax = NaN(obj.numBinsX,1);
                    
                    Tuning(icell).SEcorrModelXRef = NaN(obj.numBinsY,obj.numBinsX);
                    Tuning(icell).SEcorrModelXRefpos = NaN(obj.numBinsY,1);
                    Tuning(icell).SEcorrModelXRefmax = NaN(obj.numBinsY,1);
                    Tuning(icell).SEcorrModelXRefmaxAmp = NaN;
                    Tuning(icell).SEcorrModelXRefmaxsinAmp = NaN;
                    Tuning(icell).SEcorrModelXRefmaxOffset = NaN;
                    
                    Tuning(icell).corrModelslopeXY = [];
                    Tuning(icell).corrModelphi0XY = [];
                    Tuning(icell).corrModelrhoXY = [];
                    Tuning(icell).meancorrModelslopeXY = [];
                    Tuning(icell).meancorrModelphi0XY = [];
                    Tuning(icell).meancorrModelrhoXY = [];
                    Tuning(icell).SEcorrModelslopeXY = [];
                    Tuning(icell).SEcorrModelphi0XY = [];
                    Tuning(icell).SEcorrModelrhoXY = [];
                end
                if obj.FcomputePos
                    Tuning(icell).respModelXpos = NaN(obj.kfold,obj.numBinsY);
                    Tuning(icell).respModelYpos = NaN(obj.kfold,obj.numBinsX);
                    Tuning(icell).respModelXmax = NaN(obj.kfold,obj.numBinsY);
                    Tuning(icell).respModelYmax = NaN(obj.kfold,obj.numBinsX);
                    
                    Tuning(icell).meanrespModelXpos = NaN(obj.numBinsY,1);
                    Tuning(icell).meanrespModelYpos = NaN(obj.numBinsX,1);
                    Tuning(icell).meanrespModelXmax = NaN(obj.numBinsY,1);
                    Tuning(icell).meanrespModelYmax = NaN(obj.numBinsX,1);
                    
                    Tuning(icell).SErespModelXpos = NaN(obj.numBinsY,1);
                    Tuning(icell).SErespModelYpos = NaN(obj.numBinsX,1);
                    Tuning(icell).SErespModelXmax = NaN(obj.numBinsY,1);
                    Tuning(icell).SErespModelYmax = NaN(obj.numBinsX,1);
                    
                    Tuning(icell).fieldidxX = [];
                    Tuning(icell).fieldidxY = [];
                    
                    Tuning(icell).respModelslopeXY = NaN(1,obj.kfold);
                    Tuning(icell).respModelphi0XY = NaN(1,obj.kfold);
                    Tuning(icell).respModelrhoXY = NaN(1,obj.kfold);
                    Tuning(icell).meanrespModelslopeXY = NaN;
                    Tuning(icell).meanrespModelphi0XY = NaN;
                    Tuning(icell).meanrespModelrhoXY = NaN;
                    Tuning(icell).SErespModelslopeXY = NaN;
                    Tuning(icell).SErespModelphi0XY = NaN;
                    Tuning(icell).SErespModelrhoXY = NaN;
                end
                
                if pkfold ~= 2 && pFdiscarditer
                    if pFcomputeMarg
%                         Tuning(icell).respModel = [];
%                         Tuning(icell).respModelX = [];
%                         Tuning(icell).respModelY = [];
                    end
                    if pFcomputeCorr                    
                        Tuning(icell).corrModelX = [];
                        Tuning(icell).corrModelY = [];
                        Tuning(icell).corrModelXpos = [];
%                         Tuning(icell).corrModelXmax = [];
                        
%                         Tuning(icell).corrModelXmaxAmp = [];
%                         Tuning(icell).corrModelXmaxsinAmp = [];
%                         Tuning(icell).corrModelXmaxOffset = [];
                        Tuning(icell).corrModelYpos = [];
                        Tuning(icell).corrModelYmax = [];
                        
                        Tuning(icell).corrModelXRef = [];
                        Tuning(icell).corrModelXRefpos = [];
%                         Tuning(icell).corrModelXRefmax = [];
%                         Tuning(icell).corrModelXRefmaxAmp = [];
%                         Tuning(icell).corrModelXRefmaxsinAmp = [];
%                         Tuning(icell).corrModelXRefmaxOffset = [];

                        if pFcircularY
                            Tuning(icell).corrModelslopeXY = [];
                            Tuning(icell).corrModelphi0XY = [];
                            Tuning(icell).corrModelrhoXY = [];
                        end
                    end
                    if pFcomputePos
                        Tuning(icell).respModelXpos = [];
%                         Tuning(icell).respModelXmax = [];
                        Tuning(icell).respModelYpos = [];
                        Tuning(icell).respModelYmax = [];
                        
                        if pFcircularY
%                             Tuning(icell).respModelslopeXY = [];
%                             Tuning(icell).respModelphi0XY = [];
%                             Tuning(icell).respModelrhoXY = [];
                        end
                    end
                end
            end
            
            obj.model.EV = NaN(obj.kfold,Ncells);
            obj.model.L = NaN(obj.kfold,Ncells);
            obj.model.Q = NaN(obj.kfold,Ncells);
            obj.model.skaggs = NaN(obj.kfold,Ncells);
            obj.model.train_mean = NaN(obj.kfold,Ncells);
            obj.model.swinX = NaN(obj.kfold,Ncells);
            obj.model.swinY = NaN(obj.kfold,Ncells);
            
            pTuning = Tuning(pFgoodcells);
            pTuning_nosmth = pTuning_nosmth(pFgoodcells);
            if ~isempty(obj.tuningRef)
                TuningRef = obj.tuningRef(pFgoodcells);
            else
                TuningRef = [];
            end
            obj.tuningRef = [];
            
            goodcellsidx = find(pFgoodcells);
            Ngoodcells = numel(goodcellsidx);
            pEV = NaN(obj.kfold,Ngoodcells);
            pL = NaN(obj.kfold,Ngoodcells);
            pQ = NaN(obj.kfold,Ngoodcells);
            pskaggs = NaN(obj.kfold,Ngoodcells);
            ptrain_mean = NaN(obj.kfold,Ngoodcells);
            pswinX = NaN(obj.kfold,Ngoodcells);
            pswinY = NaN(obj.kfold,Ngoodcells);
            
            if pFshuffleXY
                spmd
                    rng(0,'combRecursive');
                end
            end
            parfor icell = 1:Ngoodcells
                if pFshuffleXY
                    stream = RandStream.getGlobalStream();
                    stream.Substream = icell;
                end
                Ycell = Y(:,min(end,goodcellsidx(icell)));
                for iter = 1:pkfold
                    if pFshuffleXY
                        spkidx = find(Z(:,goodcellsidx(icell))>0);
                        nospkidx = find(Z(:,goodcellsidx(icell))==0);
                        if numel(spkidx)>2
                            Ycell(spkidx) = Ycell(spkidx(randperm(numel(spkidx))));%Ycell(spkidx(randperm(s,numel(spkidx))));
                            Ycell(nospkidx) = Ycell(nospkidx(randperm(numel(nospkidx))));
                        end
                        kiter = 1;
                    else
                        kiter = iter;
                    end
                    [model, ~] = get2Dmap(Z(:,goodcellsidx(icell)), X(:,min(end,goodcellsidx(icell)))', Ycell', 1./T, pnumBinsX, pbinsX, pnumBinsY, pbinsY, pCVO, kiter, psampleRate, pXsmth_win, pYsmth_win, FoptiSmooth, pFcircularX, pFcircularY);%, false, X(:,min(end,goodcellsidx(icell)))',Y(:,min(end,goodcellsidx(icell))));
                    pTuning(icell).respModel(iter,:,:) = model.tuning;
                    if pFcomputeMarg
                        if isempty(prefYindex)
                            pTuning(icell).respModelX(iter,:) = nanmean(model.tuning,1);
                        else
                            pTuning(icell).respModelX(iter,:) = nanmean(model.tuning(prefYindex,:),1);
                        end
                        pTuning(icell).respModelY(iter,:) = nanmean(model.tuning,2);
                    end
                    pEV(iter,icell) = model.EV;
                    pL(iter,icell) = model.L;
                    pQ(iter,icell) = model.Q;
                    pskaggs(iter,icell) = model.skaggs;
                    ptrain_mean(iter,icell) = model.train_mean;
                    pswinX(iter,icell) = model.swinX;
                    pswinY(iter,icell) = model.swinY;
                    
                    [model, ~] = get2Dmap(Z(:,goodcellsidx(icell)), X(:,min(end,goodcellsidx(icell)))', Ycell', 1./T, pnumBinsX, pbinsX, pnumBinsY, pbinsY, pCVO, kiter, psampleRate, pXsmth_win, 0, FoptiSmooth, pFcircularX, pFcircularY);%, false, X(:,min(end,goodcellsidx(icell)))',Y(:,min(end,goodcellsidx(icell))));
                    pTuning_nosmth(icell).respModel(iter,:,:) = model.tuning;
                end
            end
            
            Xmesh = repmat(1:obj.numBinsX,[obj.numBinsY 1]);
            Ymesh = repmat((1:obj.numBinsY)',[1 obj.numBinsX]);
            Xmesh_centered = repmat((1:obj.numBinsX)-round(obj.numBinsX/2),[obj.numBinsY 1]);
            Ymesh_centered = repmat(((1:obj.numBinsY)-round(obj.numBinsY/2))',[1 obj.numBinsX]);
            outcorrXrange = [1:floor(obj.numBinsX/4) floor(obj.numBinsX*3/4)+1:obj.numBinsX];
            
            parfor icell = 1:Ngoodcells
                if ~pFshuffleXY
                    pTuning(icell).meanrespModel = squeeze(mean(pTuning(icell).respModel,1));
                    pTuning_nosmth(icell).meanrespModel = squeeze(mean(pTuning_nosmth(icell).respModel,1));
                else
                    pTuning(icell).meanrespModel = squeeze(mean(pTuning(icell).respModel(1,:,:),1));
                    pTuning_nosmth(icell).meanrespModel = squeeze(mean(pTuning_nosmth(icell).respModel(1,:,:),1));
                end
                if pFcomputeMarg || pFcomputePos
                    if isempty(prefYindex)
                        pTuning(icell).meanrespModelX = nanmean(squeeze(pTuning(icell).meanrespModel),1);
                    else
                        pTuning(icell).meanrespModelX = nanmean(squeeze(pTuning(icell).meanrespModel(prefYindex,:)),1);
                    end
                    pTuning(icell).meanrespModelY = nanmean(squeeze(pTuning(icell).meanrespModel),2);
                end
                if pFcomputeCorr
                    pTuning(icell).meancorrModelX = zeros(pnumBinsY,pnumBinsX);
                    pTuning(icell).meancorrModelXRef = zeros(pnumBinsY,pnumBinsX);
                    Xmap_norm = pTuning(icell).meanrespModelX;
                    Xmap_norm = Xmap_norm - nanmean(Xmap_norm);
                    Xmap_norm = Xmap_norm./(nanmean(Xmap_norm.^2)).^0.5;
%                     Xmap_norm = 0*pTuning_nosmth(icell).meanrespModel;
%                     ally = 1:pnumBinsY;
%                     for i_ymean = 1:pnumBinsY
%                         Xmap_norm(i_ymean,:) = nanmean(pTuning_nosmth(icell).meanrespModel(ally(~ismember(ally,i_ymean)),:),1);
%                     end
%                     Xmap_norm = Xmap_norm - repmat(nanmean(Xmap_norm,2),[1 pnumBinsX]);
%                     Xmap_norm = Xmap_norm./(repmat(nanmean(Xmap_norm.^2,2),[1 pnumBinsX])).^0.5;
                    if ~isempty(TuningRef)
                        XmapRef_norm = TuningRef(icell).meanrespModelX;
                        XmapRef_norm = XmapRef_norm - nanmean(XmapRef_norm);
                        XmapRef_norm = XmapRef_norm./(nanmean(XmapRef_norm.^2)).^0.5;
                    else
                        XmapRef_norm = Xmap_norm;
                    end
                    map_normX = pTuning(icell).meanrespModel;
%                     map_normX = pTuning_nosmth(icell).meanrespModel;
                    map_normX = (map_normX - repmat(nanmean(map_normX,2),[1,size(map_normX,2)]));
                    map_normX = map_normX./(repmat(nanmean(map_normX.^2,2),[1,size(map_normX,2)])).^0.5;
                    ishift = 0;
                    for xshift = -floor(pnumBinsX/2)+1:floor(pnumBinsX/2)
                        ishift = ishift + 1;
                        pTuning(icell).meancorrModelX(:,ishift) = (map_normX*circshift(Xmap_norm,xshift,2)')/pnumBinsX;
                        pTuning(icell).meancorrModelXRef(:,ishift) = (map_normX*circshift(XmapRef_norm,xshift,2)')/pnumBinsX;
%                         pTuning(icell).meancorrModelX(:,ishift) = diag(map_normX*circshift(Xmap_norm,xshift,2)')/pnumBinsX;
%                         pTuning(icell).meancorrModelXRef(:,ishift) = diag(map_normX*circshift(XmapRef_norm,xshift,2)')/pnumBinsX;
                    end
                    corrmap = pTuning(icell).meancorrModelX;
                    corrmap(:,outcorrXrange) = 0;
                    corrmapRef = pTuning(icell).meancorrModelXRef;
                    corrmapRef(:,outcorrXrange) = 0;
                    if pFcircularX
                        pTuning(icell).meancorrModelXpos = getCircularAverage(corrmap',0,1);
                        pTuning(icell).meancorrModelXmax = getCircularAverage(corrmap',0,maxtol,0.05);
                        [~,imaxcorrModelXmax] = max(pTuning(icell).meancorrModelXmax);
                        [~,imincorrModelXmax] = min(pTuning(icell).meancorrModelXmax);
                        pTuning(icell).meancorrModelXmaxAmp = pnumBinsX/(2*pi)*circ_dist(pTuning(icell).meancorrModelXmax(max(imaxcorrModelXmax,imincorrModelXmax))/pnumBinsX*2*pi,pTuning(icell).meancorrModelXmax(min(imaxcorrModelXmax,imincorrModelXmax))/pnumBinsX*2*pi);
                        [pTuning(icell).meancorrModelXmaxOffset,pTuning(icell).meancorrModelXmaxsinAmp] = getThetaPrecPhase(pTuning(icell).meancorrModelXmax);
                        if ~isnan(pTuning(icell).meancorrModelXmaxOffset)
                            %                             Xmax_centered = circshift(pTuning(icell).meancorrModelXmax,round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY));
                            %                             pTuning(icell).meancorrModelXmaxsinAmp = mean(Xmax_centered(1:round(pnumBinsY/2)-1)) - mean(Xmax_centered(round(pnumBinsY/2)+1:end));

                            aveAmp = zeros(1,pnumBinsY);
                            for ishift = 1:pnumBinsY
                                Resp_centered = circshift(pTuning_nosmth(icell).meanrespModel,round(pnumBinsY/2) - ishift,1);
                                Xcorr_circ = circXcorr(mean(Resp_centered(1:round(pnumBinsY/2),:),1),mean(Resp_centered(round(pnumBinsY/2)+1:end,:),1));
                                aveAmp(ishift) = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                            end
                            [pTuning(icell).meancorrModelXmaxsinAmp,maxoffset] = max(aveAmp);
                            maxoffset = getCircularAverage(aveAmp',0,maxtol,0.05);
                            pTuning(icell).meancorrModelXmaxOffset = maxoffset/pnumBinsY*360;
                            
%                             Resp_centered = circshift(pTuning(icell).meanrespModel,round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY),1);
%                             Xcorr_circ = circXcorr(mean(Resp_centered(1:round(pnumBinsY/2),:)),mean(Resp_centered(round(pnumBinsY/2)+1:end,:)));
%                             %                             Xcorr_circ(outcorrXrange) = 0;
%                             pTuning(icell).meancorrModelXmaxsinAmp = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                        end
                        pTuning(icell).meancorrModelXRefpos = getCircularAverage(corrmapRef',0,1);
                        pTuning(icell).meancorrModelXRefmax = getCircularAverage(corrmapRef',0,maxtol,0.05);
                        [~,imaxcorrModelXRefmax] = max(pTuning(icell).meancorrModelXRefmax);
                        [~,imincorrModelXRefmax] = min(pTuning(icell).meancorrModelXRefmax);
                        pTuning(icell).meancorrModelXRefmaxAmp = pnumBinsX/(2*pi)*circ_dist(pTuning(icell).meancorrModelXRefmax(max(imaxcorrModelXRefmax,imincorrModelXRefmax))/pnumBinsX*2*pi,pTuning(icell).meancorrModelXRefmax(min(imaxcorrModelXRefmax,imincorrModelXRefmax))/pnumBinsX*2*pi);
                        [pTuning(icell).meancorrModelXRefmaxOffset, pTuning(icell).meancorrModelXRefmaxsinAmp] = getThetaPrecPhase(pTuning(icell).meancorrModelXRefmax);
                    else
                        pTuning(icell).meancorrModelXpos = sum(corrmap.*Xmesh_centered,2)./sum(corrmap,2);
                        [~, pTuning(icell).meancorrModelXmax] = max(corrmap,[],2);
                        [~,imaxcorrModelXmax] = max(pTuning(icell).meancorrModelXmax);
                        [~,imincorrModelXmax] = min(pTuning(icell).meancorrModelXmax);
                        pTuning(icell).meancorrModelXmaxAmp = pTuning(icell).meancorrModelXmax(max(imaxcorrModelXmax,imincorrModelXmax))-pTuning(icell).meancorrModelXmax(min(imaxcorrModelXmax,imincorrModelXmax));
                        [pTuning(icell).meancorrModelXmaxOffset,pTuning(icell).meancorrModelXmaxsinAmp] = getThetaPrecPhase(pTuning(icell).meancorrModelXmax);
                        if ~isnan(pTuning(icell).meancorrModelXmaxOffset)
                            %                             Xmax_centered = circshift(pTuning(icell).meancorrModelXmax,round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY));
                            %                             pTuning(icell).meancorrModelXmaxsinAmp = mean(Xmax_centered(1:round(pnumBinsY/2))-1) - mean(Xmax_centered(round(pnumBinsY/2)+1:end));

                            aveAmp = zeros(1,pnumBinsY);
                            for ishift = 1:pnumBinsY
                                Resp_centered = circshift(pTuning_nosmth(icell).meanrespModel,round(pnumBinsY/2) - ishift,1);
                                Xcorr_circ = xcorr(mean(Resp_centered(1:round(pnumBinsY/2),:),1),mean(Resp_centered(round(pnumBinsY/2)+1:end,:),1));
                                aveAmp(ishift) = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                            end
                            [pTuning(icell).meancorrModelXmaxsinAmp,maxoffset] = max(aveAmp);
                            maxoffset = getCircularAverage(aveAmp',0,maxtol,0.05);
                            pTuning(icell).meancorrModelXmaxOffset = maxoffset/pnumBinsY*360;
                            
%                             Resp_centered = circshift(pTuning(icell).meanrespModel,round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY),1);
%                             Xcorr_circ = xcorr(mean(Resp_centered(1:round(pnumBinsY/2),:)),mean(Resp_centered(round(pnumBinsY/2)+1:end,:)));
%                             %                             Xcorr_circ(outcorrXrange) = 0;
%                             pTuning(icell).meancorrModelXmaxsinAmp = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                        end
                        pTuning(icell).meancorrModelXRefpos = sum(corrmapRef.*Xmesh_centered,2)./sum(corrmapRef,2);
                        [~, pTuning(icell).meancorrModelXRefmax] = max(corrmapRef,[],2);
                        [~,imaxcorrModelXRefmax] = max(pTuning(icell).meancorrModelXRefmax);
                        [~,imincorrModelXRefmax] = min(pTuning(icell).meancorrModelXRefmax);
                        pTuning(icell).meancorrModelXRefmaxAmp = pTuning(icell).meancorrModelXRefmax(max(imaxcorrModelXRefmax,imincorrModelXRefmax))-pTuning(icell).meancorrModelXRefmax(min(imaxcorrModelXRefmax,imincorrModelXRefmax));
                        [pTuning(icell).meancorrModelXRefmaxOffset, pTuning(icell).meancorrModelXRefmaxsinAmp] = getThetaPrecPhase(pTuning(icell).meancorrModelXRefmax);
                    end
                    if pFcircularY
                        %                         [slope,phi0,rho] = circularlinearfit(Ymesh(:)/pnumBinsY*2*pi,Xmesh_centered(:)/pnumBinsX,corrmap(:));
                        % %                         [slope,phi0,rho] = circularlinearfit((1:pnumBinsY)'/pnumBinsY*2*pi,(pTuning(icell).meancorrModelXmax-floor(pnumBinsX/2)),ones(pnumBinsY,1));
                        %                         pTuning(icell).meancorrModelslopeXY = slope;
                        %                         pTuning(icell).meancorrModelphi0XY = phi0*pnumBinsY/(2*pi);
                        %                         pTuning(icell).meancorrModelrhoXY = rho;
                    end

                    pTuning(icell).meancorrModelY = zeros(pnumBinsY,pnumBinsX);
                    Ymap_norm = pTuning(icell).meanrespModelY;
                    Ymap_norm = Ymap_norm - nanmean(Ymap_norm);
                    Ymap_norm = Ymap_norm./(nanmean(Ymap_norm.^2)).^0.5;
                    map_normY = pTuning(icell).meanrespModel;
                    map_normY = (map_normY - repmat(nanmean(map_normY,1),[size(map_normY,1) 1]));
                    map_normY = map_normY./(repmat(nanmean(map_normY.^2,1),[size(map_normY,1) 1])).^0.5;
                    ishift = 0;
                    for yshift = -round(pnumBinsY/2)+1:floor(pnumBinsY/2)
                        ishift = ishift + 1;
                        pTuning(icell).meancorrModelY(ishift,:) = (map_normY'*circshift(Ymap_norm,yshift,1))/pnumBinsY;
                    end
                    corrmap = pTuning(icell).meancorrModelY;
                    corrmap(:,outcorrXrange) = 0;
                    if pFcircularY
                        pTuning(icell).meancorrModelYpos = getCircularAverage(corrmap,0,1);
                        pTuning(icell).meancorrModelYmax = getCircularAverage(corrmap,0,maxtol,0.05);
                    else
                        pTuning(icell).meancorrModelYpos = sum(corrmap.*Ymesh_centered,1)./sum(corrmap,1);
                        [~, pTuning(icell).meancorrModelYmax] = max(corrmap,[],1);
                    end
                    
                    for iter = 1:pkfold
                        Xmap_norm = pTuning(icell).respModelX(iter,:);
                        Xmap_norm = Xmap_norm - nanmean(Xmap_norm);
                        Xmap_norm = Xmap_norm./(nanmean(Xmap_norm.^2)).^0.5;
                        if ~isempty(TuningRef)
                            XmapRef_norm = TuningRef(icell).meanrespModelX;
                            XmapRef_norm = XmapRef_norm - nanmean(XmapRef_norm);
                            XmapRef_norm = XmapRef_norm./(nanmean(XmapRef_norm.^2)).^0.5;
                        else
                            XmapRef_norm = Xmap_norm;
                        end
                        map_normX = squeeze(pTuning(icell).respModel(iter,:,:));
                        map_normX = (map_normX - repmat(nanmean(map_normX,2),[1,size(map_normX,2)]));
                        map_normX = map_normX./(repmat(nanmean(map_normX.^2,2),[1,size(map_normX,2)])).^0.5;
                        ishift = 0;
                        for xshift = -floor(pnumBinsX/2)+1:floor(pnumBinsX/2)
                            ishift = ishift + 1;
                            pTuning(icell).corrModelX(iter,:,ishift) = (map_normX*circshift(Xmap_norm,xshift,2)')/pnumBinsX;
                            pTuning(icell).corrModelXRef(iter,:,ishift) = (map_normX*circshift(XmapRef_norm,xshift,2)')/pnumBinsX;
                        end
                        corrmap = squeeze(pTuning(icell).corrModelX(iter,:,:));
                        corrmap(:,outcorrXrange) = 0;
                        corrmapRef = squeeze(pTuning(icell).corrModelXRef(iter,:,:));
                        corrmapRef(:,outcorrXrange) = 0;
                        if pFcircularX
                            pTuning(icell).corrModelXpos(iter,:) = getCircularAverage(squeeze(corrmap)',0,1);
                            pTuning(icell).corrModelXmax(iter,:) = getCircularAverage(squeeze(corrmap)',0,maxtol,0.05);
                            [~,imaxcorrModelXmax] = max(pTuning(icell).corrModelXmax(iter,:));
                            [~,imincorrModelXmax] = min(pTuning(icell).corrModelXmax(iter,:));
                            pTuning(icell).corrModelXmaxAmp(iter) = pnumBinsX/(2*pi)*circ_dist(pTuning(icell).corrModelXmax(iter,max(imaxcorrModelXmax,imincorrModelXmax))/pnumBinsX*2*pi,pTuning(icell).corrModelXmax(iter,min(imaxcorrModelXmax,imincorrModelXmax))/pnumBinsX*2*pi);
                            [pTuning(icell).corrModelXmaxOffset(iter), pTuning(icell).corrModelXmaxsinAmp(iter)] = getThetaPrecPhase(pTuning(icell).corrModelXmax(iter,:),pTuning(icell).meancorrModelXmaxOffset);
                            if ~isnan(pTuning(icell).corrModelXmaxOffset(iter))
                                %                                 Xmax_centered = circshift(pTuning(icell).corrModelXmax(iter,:),round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY));
                                %                                 pTuning(icell).corrModelXmaxsinAmp(iter) = mean(Xmax_centered(1:round(pnumBinsY/2))) - mean(Xmax_centered(round(pnumBinsY/2):end));

                                %                                 Resp_centered = circshift(squeeze(pTuning(icell).respModel(iter,:,:)),round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY),1);

                                aveAmp = zeros(1,pnumBinsY);
                                for ishift = 1:pnumBinsY
                                    Resp_centered = circshift(squeeze(pTuning_nosmth(icell).respModel(iter,:,:)),round(pnumBinsY/2) - ishift,1);
                                    Xcorr_circ = circXcorr(mean(Resp_centered(1:round(pnumBinsY/2),:),1),mean(Resp_centered(round(pnumBinsY/2)+1:end,:),1));
                                    aveAmp(ishift) = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                                end
                                [pTuning(icell).corrModelXmaxsinAmp(iter),maxoffset] = max(aveAmp);
                                maxoffset = getCircularAverage(aveAmp',0,maxtol,0.05);
                                pTuning(icell).corrModelXmaxOffset(iter) = maxoffset/pnumBinsY*360;

%                                 Resp_centered = circshift(squeeze(pTuning(icell).respModel(iter,:,:)),round(pnumBinsY/2) - round(pTuning(icell).corrModelXmaxOffset(iter)/360*pnumBinsY),1);
%                                 Xcorr_circ = circXcorr(mean(Resp_centered(1:round(pnumBinsY/2),:)),mean(Resp_centered(round(pnumBinsY/2)+1:end,:)));
%                                 %                                 Xcorr_circ(outcorrXrange) = 0;
%                                 pTuning(icell).corrModelXmaxsinAmp(iter) = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                            end
                            pTuning(icell).corrModelXRefpos(iter,:) = getCircularAverage(squeeze(corrmapRef)',0,1);
                            pTuning(icell).corrModelXRefmax(iter,:) = getCircularAverage(squeeze(corrmapRef)',0,maxtol,0.05);
                            [~,imaxcorrModelXRefmax] = max(pTuning(icell).corrModelXRefmax(iter,:));
                            [~,imincorrModelXRefmax] = min(pTuning(icell).corrModelXRefmax(iter,:));
                            pTuning(icell).corrModelXRefmaxAmp(iter) = pnumBinsX/(2*pi)*circ_dist(pTuning(icell).corrModelXRefmax(iter,max(imaxcorrModelXRefmax,imincorrModelXRefmax))/pnumBinsX*2*pi,pTuning(icell).corrModelXRefmax(iter,min(imaxcorrModelXRefmax,imincorrModelXRefmax))/pnumBinsX*2*pi);
                            [pTuning(icell).corrModelXRefmaxOffset(iter), pTuning(icell).corrModelXRefmaxsinAmp(iter)] = getThetaPrecPhase(pTuning(icell).corrModelXRefmax(iter,:));%,pTuning(icell).meancorrModelXRefmaxOffset);
                        else
                            pTuning(icell).corrModelXpos(iter,:) = sum(squeeze(corrmap).*Xmesh_centered,2)./sum(squeeze(corrmap),2);
                            [~, pTuning(icell).corrModelXmax(iter,:)] = max(squeeze(corrmap),[],2);
                            [~,imaxcorrModelXmax] = max(pTuning(icell).corrModelXmax(iter,:));
                            [~,imincorrModelXmax] = min(pTuning(icell).corrModelXmax(iter,:));
                            pTuning(icell).corrModelXmaxAmp(iter) = pTuning(icell).corrModelXmax(iter,max(imaxcorrModelXmax,imincorrModelXmax)) - pTuning(icell).corrModelXmax(iter,min(imaxcorrModelXmax,imincorrModelXmax));
                            [pTuning(icell).corrModelXmaxOffset(iter), pTuning(icell).corrModelXmaxsinAmp(iter)] = getThetaPrecPhase(pTuning(icell).corrModelXmax(iter,:),pTuning(icell).meancorrModelXmaxOffset);
                            if ~isnan(pTuning(icell).corrModelXmaxOffset(iter))
                                %                                 Xmax_centered = circshift(pTuning(icell).corrModelXmax(iter,:),round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY));
                                %                                 pTuning(icell).corrModelXmaxsinAmp(iter) = mean(Xmax_centered(1:round(pnumBinsY/2))) - mean(Xmax_centered(round(pnumBinsY/2):end));

                                %                                 Resp_centered = circshift(squeeze(pTuning(icell).respModel(iter,:,:)),round(pnumBinsY/2) - round(pTuning(icell).meancorrModelXmaxOffset/360*pnumBinsY),1);

                                aveAmp = zeros(1,pnumBinsY);
                                for ishift = 1:pnumBinsY
                                    Resp_centered = circshift(squeeze(pTuning_nosmth(icell).respModel(iter,:,:)),round(pnumBinsY/2) - ishift,1);
                                    Xcorr_circ = xcorr(mean(Resp_centered(1:round(pnumBinsY/2),:)),mean(Resp_centered(round(pnumBinsY/2)+1:end,:)));
                                    aveAmp(ishift) = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                                end
                                [pTuning(icell).corrModelXmaxsinAmp(iter),maxoffset] = max(aveAmp);
                                maxoffset = getCircularAverage(aveAmp',0,maxtol,0.05);
                                pTuning(icell).corrModelXmaxOffset(iter) = maxoffset/pnumBinsY*360;

%                                 Resp_centered = circshift(squeeze(pTuning(icell).respModel(iter,:,:)),round(pnumBinsY/2) - round(pTuning(icell).corrModelXmaxOffset(iter)/360*pnumBinsY),1);
%                                 Xcorr_circ = xcorr(mean(Resp_centered(1:round(pnumBinsY/2),:)),mean(Resp_centered(round(pnumBinsY/2)+1:end,:)));
%                                 %                                 Xcorr_circ(outcorrXrange) = 0;
%                                 pTuning(icell).corrModelXmaxsinAmp(iter) = getCircularAverage(Xcorr_circ(:),0,maxtol,0.05) - pnumBinsX/2;
                            end
                            pTuning(icell).corrModelXRefpos(iter,:) = sum(squeeze(corrmapRef).*Xmesh_centered,2)./sum(squeeze(corrmapRef),2);
                            [~, pTuning(icell).corrModelXRefmax(iter,:)] = max(squeeze(corrmapRef),[],2);
                            [~,imaxcorrModelXRefmax] = max(pTuning(icell).corrModelXRefmax(iter,:));
                            [~,imincorrModelXRefmax] = min(pTuning(icell).corrModelXRefmax(iter,:));
                            pTuning(icell).corrModelXRefmaxAmp(iter) = pTuning(icell).corrModelXRefmax(iter,max(imaxcorrModelXRefmax,imincorrModelXRefmax)) - pTuning(icell).corrModelXRefmax(iter,min(imaxcorrModelXRefmax,imincorrModelXRefmax));
                            [pTuning(icell).corrModelXRefmaxOffset(iter), pTuning(icell).corrModelXRefmaxsinAmp(iter)] = getThetaPrecPhase(pTuning(icell).corrModelXRefmax(iter,:));
                        end
                        if pFcircularY
                            %                             [slope,phi0,rho] = circularlinearfit(Ymesh(:)/pnumBinsY*2*pi,Xmesh_centered(:)/pnumBinsX,corrmap(:));
                            % %                             [slope,phi0,rho] = circularlinearfit((1:pnumBinsY)'/pnumBinsY*2*pi,(pTuning(icell).corrModelXmax(iter,:)-floor(pnumBinsX/2))',ones(pnumBinsY,1));
                            %                             pTuning(icell).corrModelslopeXY(iter) = slope;
                            %                             pTuning(icell).corrModelphi0XY(iter) = phi0*pnumBinsY/(2*pi);
                            %                             pTuning(icell).corrModelrhoXY(iter) = rho;
                        end

                        Ymap_norm = (pTuning(icell).respModelY(iter,:))';
                        Ymap_norm = Ymap_norm - nanmean(Ymap_norm);
                        Ymap_norm = Ymap_norm./(nanmean(Ymap_norm.^2)).^0.5;
                        map_normY = squeeze(pTuning(icell).respModel(iter,:,:));
                        map_normY = (map_normY - repmat(nanmean(map_normY,2),[1,size(map_normY,2)]));
                        map_normY = map_normY./(repmat(nanmean(map_normY.^2,2),[1,size(map_normY,2)])).^0.5;
                        ishift = 0;
                        for yshift = -round(pnumBinsY/2)+1:floor(pnumBinsY/2)
                            ishift = ishift + 1;
                            pTuning(icell).corrModelY(iter,ishift,:) = (map_normY'*circshift(Ymap_norm,yshift,1))/pnumBinsX;
                        end
                        corrmap = squeeze(pTuning(icell).corrModelY(iter,:,:));
                        corrmap(:,outcorrXrange) = 0;
                        if pFcircularY
                            pTuning(icell).corrModelYpos(iter,:) = getCircularAverage(squeeze(corrmap),0,1);
                            pTuning(icell).corrModelYmax(iter,:) = getCircularAverage(squeeze(corrmap),0,maxtol,0.05);
                        else
                            pTuning(icell).corrModelYpos(iter,:) = sum(squeeze(corrmap).*Ymesh_centered,1)./sum(squeeze(corrmap),1);
                            [~, pTuning(icell).corrModelYmax(iter,:)] = max(squeeze(corrmap),[],1);
                        end

                    end
                end
                if pFcomputePos
                    [q_idxX,fieldX] = findfield(pTuning(icell).meanrespModelX,pqthreshold);
                    pTuning(icell).fieldidxX = fieldX;
                    outfieldX = find(~ismember(1:pnumBinsX,fieldX));
                    [q_idxY,fieldY] = findfield(pTuning(icell).meanrespModelY,pqthreshold);
                    pTuning(icell).fieldidxY = fieldY;
                    outfieldY = find(~ismember(1:pnumBinsY,fieldY));
                    [~,fieldXmax] = max(pTuning(icell).meanrespModelX);

                    tuningX0 = squeeze(pTuning(icell).meanrespModel);
                    tuningX = tuningX0;
                    if ~isempty(fieldX)
                        tuningX(:,outfieldX) = 0;%repmat(min(tuningX(:,fieldX),[],2),[1 numel(outfieldX)]);
                    end
                    if pFcircularX
                        pTuning(icell).meanrespModelXpos = getCircularAverage(tuningX',0,1);
                        pTuning(icell).meanrespModelXmax = getCircularAverage(tuningX0',0,maxtol,0.05);
                    else
                        pTuning(icell).meanrespModelXpos = sum(tuningX.*Xmesh,2)./sum(tuningX,2);
                        [~, pTuning(icell).meanrespModelXmax] = max(tuningX0,[],2);
                    end

                    tuningY0 = squeeze(pTuning(icell).meanrespModel);
                    tuningY = tuningY0;
                    if ~isempty(fieldX)
                        tuningY(:,outfieldX) = 0;%repmat(min(tuningY(fieldY,:),[],1),[numel(outfieldY) 1]);
                    end
                    if pFcircularY
                        pTuning(icell).meanrespModelYpos = getCircularAverage(tuningY,0,1);
                        pTuning(icell).meanrespModelYmax = getCircularAverage(tuningY0,0,maxtol,0.05);

                        tuningY_centered = circshift(tuningY,-fieldXmax+floor(pnumBinsX/2),2);
                        [slope,phi0,rho] = circularlinearfit(Ymesh(:)/pnumBinsY*2*pi,Xmesh_centered(:)/pnumBinsX,tuningY_centered(:));
                        pTuning(icell).meanrespModelslopeXY = slope;
                        pTuning(icell).meanrespModelphi0XY = phi0*pnumBinsY/(2*pi);
                        pTuning(icell).meanrespModelrhoXY = rho;
                    else
                        pTuning(icell).meanrespModelYpos = sum(tuningY.*Ymesh,1)./sum(tuningY,1);
                        [~, pTuning(icell).meanrespModelYmax] = max(tuningY0,[],1);
                    end

                    for iter = 1:pkfold
                        tuningX0 = squeeze(pTuning(icell).respModel(iter,:,:));
                        tuningX = tuningX0;
                        if ~isempty(fieldX)
                            tuningX(:,outfieldX) = 0;%repmat(min(tuningX(:,fieldX),[],2),[1 numel(outfieldX)]);
                        end

                        if pFcircularX
                            pTuning(icell).respModelXpos(iter,:) = getCircularAverage(tuningX',0,1);
                            pTuning(icell).respModelXmax(iter,:) = getCircularAverage(tuningX0',0,maxtol,0.05);
                        else
                            pTuning(icell).respModelXpos(iter,:) = sum(tuningX.*Xmesh,1)./sum(tuningX,2);
                            [~, pTuning(icell).respModelXmax(iter,:)] = max(tuningX0,[],2);
                        end

                        tuningY0 = squeeze(pTuning(icell).respModel(iter,:,:));
                        tuningY = tuningY0;
                        if ~isempty(fieldX)
                            tuningY(:,outfieldX) = 0;%repmat(min(tuningY(fieldY,:),[],1),[numel(outfieldY) 1]);
                        end
                        if pFcircularY
                            pTuning(icell).respModelYpos(iter,:) = getCircularAverage(tuningY,0,1);
                            pTuning(icell).respModelYmax(iter,:) = getCircularAverage(tuningY0,0,maxtol,0.05);

                            tuningY_centered = circshift(tuningY,-fieldXmax+floor(pnumBinsX/2),2);
                            [slope,phi0,rho] = circularlinearfit(Ymesh(:)/pnumBinsY*2*pi,Xmesh_centered(:)/pnumBinsX,tuningY_centered(:));
                            pTuning(icell).respModelslopeXY(iter) = slope;
                            pTuning(icell).respModelphi0XY(iter) = phi0*pnumBinsY/(2*pi);
                            pTuning(icell).respModelrhoXY(iter) = rho;
                        else
                            pTuning(icell).respModelYpos(iter,:) = sum(tuningY.*Ymesh,1)./sum(tuningY,1);
                            [~, pTuning(icell).respModelYmax(iter,:)] = max(tuningY0,[],1);
                        end
                    end
                end
                stdresp = 0;
                stdrespX = 0;
                stdrespY = 0;
                stdrespXpos = 0;
                stdrespXmax = 0;
                stdrespYpos = 0;
                stdrespYmax = 0;
                stdslopeXY = 0;
                stdphi0XY = 0;
                stdrhoXY = 0;
                stdcorrX = 0;
                stdcorrXpos = 0;
                stdcorrXmax = 0;
                stdcorrXmaxAmp = 0;
                stdcorrXmaxsinAmp = 0;
                stdcorrXmaxOffset = 0;
                stdcorrY = 0;
                stdcorrYpos = 0;
                stdcorrYmax = 0;
                stdcorrXRef = 0;
                stdcorrXRefpos = 0;
                stdcorrXRefmax = 0;
                stdcorrXRefmaxAmp = 0;
                stdcorrXRefmaxsinAmp = 0;
                stdcorrXRefmaxOffset = 0;
                stdcorrslopeXY = 0;
                stdcorrphi0XY = 0;
                stdcorrrhoXY = 0;
                if pFshuffleXY
                    StdSEfactor = 1/pkfold;
                else
                    StdSEfactor = (pkfold - 1)/pkfold;
                end
                if ~pFshuffleXY && pkfold > 2
                    for i = 1:pkfold
                        r_iter = squeeze(pTuning(icell).respModel(i,:,:));
                        r_ave = pTuning(icell).meanrespModel;
                        stdresp = stdresp + StdSEfactor*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                        if pFcomputeMarg
                            r_iter = squeeze(pTuning(icell).respModelX(i,:));
                            r_ave = pTuning(icell).meanrespModelX;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespX = stdrespX + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                            r_iter = squeeze(pTuning(icell).respModelY(i,:));
                            r_ave = pTuning(icell).meanrespModelY;
                            r_iter = r_iter(:);
                            r_ave = r_ave(:);
                            stdrespY = stdrespY + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                        end
                        if pFcomputeCorr
                            r_iter = squeeze(pTuning(icell).corrModelX(i,:,:));
                            r_ave = pTuning(icell).meancorrModelX;
                            stdcorrX = stdcorrX + StdSEfactor*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                            r_iter = squeeze(pTuning(icell).corrModelXRef(i,:,:));
                            r_ave = pTuning(icell).meancorrModelXRef;
                            stdcorrXRef = stdcorrXRef + StdSEfactor*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                            r_iter = squeeze(pTuning(icell).corrModelY(i,:,:));
                            r_ave = pTuning(icell).meancorrModelY;
                            stdcorrY = stdcorrY + StdSEfactor*((r_iter - mean(r_iter(:))) - (r_ave - mean(r_ave(:)))).^2;
                            if pFcircularX
                                r_iter = (pTuning(icell).corrModelXpos(i,:)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXpos/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(pTuning(icell).respModelXpos(i,:)' - pTuning(icell).meanrespModelXpos);
                                stdcorrXpos = stdcorrXpos + StdSEfactor*iterdist.^2;
                                r_iter = (pTuning(icell).corrModelXmax(i,:)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXmax/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(pTuning(icell).respModelXmax(i,:)' - pTuning(icell).meanrespModelXmax);
                                stdcorrXmax = stdcorrXmax + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXmaxAmp(i)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXmaxAmp/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = abs(circ_dist(r_iter,r_ave)*pnumBinsX/(2*pi));
                                stdcorrXmaxAmp = stdcorrXmaxAmp + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXmaxsinAmp(i)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXmaxsinAmp/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = abs(circ_dist(r_iter,r_ave)*pnumBinsX/(2*pi));
                                stdcorrXmaxsinAmp = stdcorrXmaxsinAmp + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXmaxOffset(i)/360*2*pi);
                                r_ave = pTuning(icell).meancorrModelXmaxOffset/360*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter,r_ave)*360/(2*pi);
                                stdcorrXmaxOffset = stdcorrXmaxOffset + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXRefpos(i,:)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXRefpos/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(pTuning(icell).respModelXpos(i,:)' - pTuning(icell).meanrespModelXpos);
                                stdcorrXRefpos = stdcorrXRefpos + StdSEfactor*iterdist.^2;
                                r_iter = (pTuning(icell).corrModelXRefmax(i,:)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXRefmax/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(pTuning(icell).respModelXmax(i,:)' - pTuning(icell).meanrespModelXmax);
                                stdcorrXRefmax = stdcorrXRefmax + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXRefmaxAmp(i)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXRefmaxAmp/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = abs(circ_dist(r_iter,r_ave)*pnumBinsX/(2*pi));
                                stdcorrXRefmaxAmp = stdcorrXRefmaxAmp + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXRefmaxsinAmp(i)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meancorrModelXRefmaxsinAmp/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = abs(circ_dist(r_iter,r_ave)*pnumBinsX/(2*pi));
                                stdcorrXRefmaxsinAmp = stdcorrXRefmaxsinAmp + StdSEfactor*iterdist.^2;

                                r_iter = (pTuning(icell).corrModelXRefmaxOffset(i)/360*2*pi);
                                r_ave = pTuning(icell).meancorrModelXRefmaxOffset/360*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter,r_ave)*360/(2*pi);
                                stdcorrXRefmaxOffset = stdcorrXRefmaxOffset + StdSEfactor*iterdist.^2;
                            else
                                r_iter = pTuning(icell).corrModelXpos(i,:)';
                                r_ave = pTuning(icell).meancorrModelXpos;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXpos = stdcorrXpos + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                                r_iter = pTuning(icell).corrModelXmax(i,:)';
                                r_ave = pTuning(icell).meancorrModelXmax;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXmax = stdcorrXmax + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                                r_iter = pTuning(icell).corrModelXmaxAmp(i);
                                r_ave = pTuning(icell).meancorrModelXmaxAmp;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXmaxAmp = stdcorrXmaxAmp + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                                r_iter = pTuning(icell).corrModelXmaxsinAmp(i);
                                r_ave = pTuning(icell).meancorrModelXmaxsinAmp;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXmaxsinAmp = stdcorrXmaxsinAmp + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                                r_iter = (pTuning(icell).corrModelXmaxOffset(i)/360*2*pi);
                                r_ave = pTuning(icell).meancorrModelXmaxOffset/360*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter,r_ave)*360/(2*pi);
                                stdcorrXmaxOffset = stdcorrXmaxOffset + StdSEfactor*iterdist.^2;

                                r_iter = pTuning(icell).corrModelXRefpos(i,:)';
                                r_ave = pTuning(icell).meancorrModelXRefpos;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXRefpos = stdcorrXRefpos + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                                r_iter = pTuning(icell).corrModelXRefmax(i,:)';
                                r_ave = pTuning(icell).meancorrModelXRefmax;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXRefmax = stdcorrXRefmax + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                                r_iter = pTuning(icell).corrModelXRefmaxAmp(i);
                                r_ave = pTuning(icell).meancorrModelXRefmaxAmp;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXRefmaxAmp = stdcorrXRefmaxAmp + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                                r_iter = pTuning(icell).corrModelXRefmaxsinAmp(i);
                                r_ave = pTuning(icell).meancorrModelXRefmaxsinAmp;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrXRefmaxsinAmp = stdcorrXRefmaxsinAmp + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                                r_iter = (pTuning(icell).corrModelXRefmaxOffset(i)/360*2*pi);
                                r_ave = pTuning(icell).meancorrModelXRefmaxOffset/360*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter,r_ave)*360/(2*pi);
                                stdcorrXRefmaxOffset = stdcorrXRefmaxOffset + StdSEfactor*iterdist.^2;

                            end
                            %                         if pFcircularY
                            %                             stdcorrslopeXY = stdcorrslopeXY + StdSEfactor*(pTuning(icell).corrModelslopeXY(i) - pTuning(icell).meancorrModelslopeXY).^2;
                            %                             stdcorrphi0XY = stdcorrphi0XY + StdSEfactor*(pTuning(icell).corrModelphi0XY(i) - pTuning(icell).meancorrModelphi0XY).^2;
                            %                             stdcorrrhoXY = stdcorrrhoXY + StdSEfactor*(pTuning(icell).corrModelrhoXY(i) - pTuning(icell).meancorrModelrhoXY).^2;
                            %                         end
                            if pFcircularY
                                r_iter = (pTuning(icell).corrModelYpos(i,:)/pnumBinsY*2*pi);
                                r_ave = pTuning(icell).meancorrModelYpos/pnumBinsY*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(pTuning(icell).respModelXpos(i,:)' - pTuning(icell).meanrespModelXpos);
                                stdcorrYpos = stdcorrYpos + StdSEfactor*iterdist.^2;
                                r_iter = (pTuning(icell).corrModelYmax(i,:)/pnumBinsY*2*pi);
                                r_ave = pTuning(icell).meancorrModelYmax/pnumBinsY*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(pTuning(icell).respModelXmax(i,:)' - pTuning(icell).meanrespModelXmax);
                                stdcorrYmax = stdcorrYmax + StdSEfactor*iterdist.^2;

                            else
                                r_iter = pTuning(icell).corrModelYpos(i,:)';
                                r_ave = pTuning(icell).meancorrModelYpos;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrYpos = stdcorrYpos + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                                r_iter = pTuning(icell).corrModelYmax(i,:)';
                                r_ave = pTuning(icell).meancorrModelYmax;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdcorrYmax = stdcorrYmax + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                            end
                        end
                        if pFcomputePos
                            if pFcircularX
                                r_iter = (pTuning(icell).respModelXpos(i,:)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meanrespModelXpos/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(pTuning(icell).respModelXpos(i,:)' - pTuning(icell).meanrespModelXpos);
                                stdrespXpos = stdrespXpos + StdSEfactor*iterdist.^2;
                                r_iter = (pTuning(icell).respModelXmax(i,:)/pnumBinsX*2*pi);
                                r_ave = pTuning(icell).meanrespModelXmax/pnumBinsX*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsX/(2*pi);%abs(pTuning(icell).respModelXmax(i,:)' - pTuning(icell).meanrespModelXmax);
                                stdrespXmax = stdrespXmax + StdSEfactor*iterdist.^2;

                            else
                                r_iter = pTuning(icell).respModelXpos(i,:)';
                                r_ave = pTuning(icell).meanrespModelXpos;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdrespXpos = stdrespXpos + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                                r_iter = pTuning(icell).respModelXmax(i,:)';
                                r_ave = pTuning(icell).meanrespModelXmax;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdrespXmax = stdrespXmax + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                            end
                            if pFcircularY
                                r_iter = (pTuning(icell).respModelYpos(i,:)/pnumBinsY*2*pi);
                                r_ave = pTuning(icell).meanrespModelYpos/pnumBinsY*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(pTuning(icell).respModelXpos(i,:)' - pTuning(icell).meanrespModelXpos);
                                stdrespYpos = stdrespYpos + StdSEfactor*iterdist.^2;
                                r_iter = (pTuning(icell).respModelYmax(i,:)/pnumBinsY*2*pi);
                                r_ave = pTuning(icell).meanrespModelYmax/pnumBinsY*2*pi;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                iterdist = circ_dist(r_iter - circ_mean(r_iter,[],1),r_ave - circ_mean(r_ave,[],1))*pnumBinsY/(2*pi);%abs(pTuning(icell).respModelXmax(i,:)' - pTuning(icell).meanrespModelXmax);
                                stdrespYmax = stdrespYmax + StdSEfactor*iterdist.^2;


                                stdslopeXY = stdslopeXY + StdSEfactor*(pTuning(icell).respModelslopeXY(i) - pTuning(icell).meanrespModelslopeXY).^2;
                                stdphi0XY = stdphi0XY + StdSEfactor*(pTuning(icell).respModelphi0XY(i) - pTuning(icell).meanrespModelphi0XY).^2;
                                stdrhoXY = stdrhoXY + StdSEfactor*(pTuning(icell).respModelrhoXY(i) - pTuning(icell).meanrespModelrhoXY).^2;
                            else
                                r_iter = pTuning(icell).respModelYpos(i,:)';
                                r_ave = pTuning(icell).meanrespModelYpos;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdrespYpos = stdrespYpos + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;
                                r_iter = pTuning(icell).respModelYmax(i,:)';
                                r_ave = pTuning(icell).meanrespModelYmax;
                                r_iter = r_iter(:);
                                r_ave = r_ave(:);
                                stdrespYmax = stdrespYmax + StdSEfactor*((r_iter - mean(r_iter)) - (r_ave - mean(r_ave))).^2;

                            end
                        end
                    end
                    
                    pTuning(icell).SErespModel = stdresp.^0.5;
                    if pFcomputeMarg
                        pTuning(icell).SErespModelX = stdrespX.^0.5;
                        pTuning(icell).SErespModelY = stdrespY.^0.5;
                    end
                    if pFcomputeCorr
                        pTuning(icell).SEcorrModelX = stdcorrX.^0.5;
                        pTuning(icell).SEcorrModelY = stdcorrY.^0.5;
                        pTuning(icell).SEcorrModelXpos = stdcorrXpos.^0.5;
                        pTuning(icell).SEcorrModelYpos = stdcorrYpos.^0.5;
                        pTuning(icell).SEcorrModelXmax = stdcorrXmax.^0.5;
                        pTuning(icell).SEcorrModelXmaxAmp = stdcorrXmaxAmp.^0.5;
                        pTuning(icell).SEcorrModelXmaxsinAmp = stdcorrXmaxsinAmp.^0.5;
                        pTuning(icell).SEcorrModelXmaxOffset = stdcorrXmaxOffset.^0.5;
                        pTuning(icell).SEcorrModelYmax = stdcorrYmax.^0.5;
                        
                        pTuning(icell).SEcorrModelXRef = stdcorrXRef.^0.5;
                        pTuning(icell).SEcorrModelXRefpos = stdcorrXRefpos.^0.5;
                        pTuning(icell).SEcorrModelXRefmax = stdcorrXRefmax.^0.5;
                        pTuning(icell).SEcorrModelXRefmaxAmp = stdcorrXRefmaxAmp.^0.5;
                        pTuning(icell).SEcorrModelXRefmaxsinAmp = stdcorrXRefmaxsinAmp.^0.5;
                        pTuning(icell).SEcorrModelXRefmaxOffset = stdcorrXRefmaxOffset.^0.5;
                        
                        if pFcircularY
                            %                         pTuning(icell).SEcorrModelslopeXY = stdcorrslopeXY.^0.5;
                            %                         pTuning(icell).SEcorrModelphi0XY = stdcorrphi0XY.^0.5;
                            %                         pTuning(icell).SEcorrModelrhoXY = stdcorrrhoXY.^0.5;
                        end
                    end
                    
                    if pFcomputePos
                        pTuning(icell).SErespModelXpos = stdrespXpos.^0.5;
                        pTuning(icell).SErespModelYpos = stdrespYpos.^0.5;
                        pTuning(icell).SErespModelXmax = stdrespXmax.^0.5;
                        pTuning(icell).SErespModelYmax = stdrespYmax.^0.5;
                        
                        
                        if pFcircularY
                            pTuning(icell).SErespModelslopeXY = stdslopeXY.^0.5;
                            pTuning(icell).SErespModelphi0XY = stdphi0XY.^0.5;
                            pTuning(icell).SErespModelrhoXY = stdrhoXY.^0.5;
                        end
                    end
                end

                %we empty those to save memory. if necessary, SE of these 
                %components can be recovered from the estimates of the 
                %full map
                if pkfold ~= 2 && pFdiscarditer
                    pTuning(icell).respModel = [];
                    if pFcomputeMarg
%                         pTuning(icell).respModelX = [];
%                         pTuning(icell).respModelY = [];
                    end
                    if pFcomputeCorr
                        pTuning(icell).corrModelX = [];
                        pTuning(icell).corrModelY = [];
                        pTuning(icell).corrModelXpos = [];
%                         pTuning(icell).corrModelXmax = [];
                        
%                         pTuning(icell).corrModelXmaxAmp = [];
%                         pTuning(icell).corrModelXmaxsinAmp = [];
%                         pTuning(icell).corrModelXmaxOffset = [];
                        pTuning(icell).corrModelYpos = [];
                        pTuning(icell).corrModelYmax = [];
                        
                        pTuning(icell).corrModelXRef = [];
                        pTuning(icell).corrModelXRefpos = [];
%                         pTuning(icell).corrModelXRefmax = [];
%                         pTuning(icell).corrModelXRefmaxAmp = [];
%                         pTuning(icell).corrModelXRefmaxsinAmp = [];
%                         pTuning(icell).corrModelXRefmaxOffset = [];
                        
                        if pFcircularY
                            pTuning(icell).corrModelslopeXY = [];
                            pTuning(icell).corrModelphi0XY = [];
                            pTuning(icell).corrModelrhoXY = [];
                        end
                    end
                    if pFcomputePos
                        pTuning(icell).respModelXpos = [];
%                         pTuning(icell).respModelXmax = [];
                        pTuning(icell).respModelYpos = [];
                        pTuning(icell).respModelYmax = [];
                        
                        if pFcircularY
%                             pTuning(icell).respModelslopeXY = [];
%                             pTuning(icell).respModelphi0XY = [];
%                             pTuning(icell).respModelrhoXY = [];
                        end
                    end
                end
            end
            Tuning(pFgoodcells) = pTuning;
            obj.tuningRef = [];
            obj.model.tuning = Tuning;
            obj.model.EV(:,pFgoodcells) = pEV;
            obj.model.L(:,pFgoodcells) = pL;
            obj.model.Q(:,pFgoodcells) = pQ;
            obj.model.skaggs(:,pFgoodcells) = pskaggs;
            obj.model.train_mean(:,pFgoodcells) = ptrain_mean;
            obj.model.swinX(:,pFgoodcells) = pswinX;
            obj.model.swinY(:,pFgoodcells) = pswinY;
            
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

