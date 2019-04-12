function UpdateplotBehavSpks(PlotVar,EXP,Layout)
PlotVar.Plots = [];

page = PlotVar.Page;
win = PlotVar.Window;
Layout.DivideWindow(page, win, size(PlotVar.ChosenObj,1)*size(PlotVar.ChosenCell,2)*size(PlotVar.ChosenContrast,2), size(PlotVar.ChosenGain,2)*size(PlotVar.ChosenRoomlength,2)*size(PlotVar.ChosenOutcome,2));

nplot =  0;
nplotfast1Dmap = 0;
maxfast1Dmap = 0;

trialID = EXP.Nav.trialID;
for c = 1:size(PlotVar.ChosenContrast,2)
    for g = 1:size(PlotVar.ChosenGain,2)
        for r = 1:size(PlotVar.ChosenRoomlength,2)
            for o = 1:size(PlotVar.ChosenOutcome,2)
                EXP.Nav.trialID = trialID;
                Fnoblanks = true;
                Fnoafterblanks = true;
                tidx = GetSubsets(EXP.Nav,PlotVar.ChosenContrast(:,c),PlotVar.ChosenGain(:,g),PlotVar.ChosenRoomlength(:,r),PlotVar.ChosenOutcome(:,o),PlotVar.speed_th, Fnoblanks, Fnoafterblanks);
                tidxAllphs = tidx;
                trialIDnum = unique(EXP.Nav.trialID(tidx));
                trialID = EXP.Nav.trialID;
                trialID(~ismember(trialID,trialIDnum)) = 0;
                for itr = 1:numel(trialIDnum)
                    trialID(trialID == trialIDnum(itr)) = itr;
                end
                
                for cell = 1:size(PlotVar.ChosenCell,2)                    
                    icell = PlotVar.ChosenCell(:,cell);
                    spktrain = EXP.Spk.spikeTrain(:,icell);
                    spktrain = sum(spktrain,2);
                    
                    Fdispmat = PlotVar.Fdispmat;
                    if Fdispmat
                        spkidx = true(size(spktrain));
                    else
                        spkidx = spktrain>0;
                    end
                    
                    try
                        varX = EXP.Nav.(PlotVar.ChosenVarX{1})(:);
                        varY = EXP.Nav.(PlotVar.ChosenVarY{1})(:);
                    catch
                        error('selected variable not a field of the EXP.data.es structure')
                    end
                    
                    for k = 1:size(PlotVar.ChosenObj,1)
                        wintitle = ['Cell: ' num2str(PlotVar.ChosenCell(:,cell)') ' Contrast: ' num2str(EXP.Nav.contrastValues(PlotVar.ChosenContrast(:,c)')) ' Gain: ' num2str(EXP.Nav.gainValues(PlotVar.ChosenGain(:,g)'))];
                        Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}.Title = wintitle;
                                                
                        switch PlotVar.ChosenObj{k}                            
                            case 'Raster'
                                nplot = nplot + 1;
                                PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                varxx = varX(tidx & EXP.Nav.outcome == 2);
                                varyy = varY(tidx & EXP.Nav.outcome == 2);
                                
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & EXP.Nav.outcome == 2, varX, varY, spktrain, Fdispmat, [0 0 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & EXP.Nav.outcome == 3, varX, varY, spktrain, Fdispmat, [1 0 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & EXP.Nav.outcome == 4, varX, varY, spktrain, Fdispmat, [0 1 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & EXP.Nav.outcome == 0, varX, varY, spktrain, Fdispmat, [1 0 0]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & EXP.Nav.outcome == 1, varX, varY, spktrain, Fdispmat, [0 0 1]);
                                PlotBehavSpks(PlotVar.Plots{nplot}, spkidx & tidx & EXP.Nav.outcome == 5, varX, varY, spktrain, Fdispmat, [0.5 0.5 0.5]);
                            case 'fast1Dmap' 
                                varX = EXP.Nav.(PlotVar.ChosenVarX{1})(:);
                                varX = normalise1var(varX,100,[],[0 100]);
                                dx = 1;
                                samplerate = 1./EXP.Nav.sampleSize(tidx);
                                nbXbinsmth = round(1/(4/100));
                                spktrain(~isnan(spktrain)) = smthInTime(spktrain(~isnan(spktrain)), mean(samplerate), 15, 'same', [], 'boxcar_centered');
                                
                                [map,~] = fast1Dmap(varX(tidx), spktrain(tidx), dx, samplerate,nbXbinsmth,EXP.Nav.CircularMaze);
                                maxfast1Dmap = max(max(map),maxfast1Dmap);
                                x = 1:numel(map);
                                nplot = nplot + 1;
                                if nplotfast1Dmap == 0
                                    nplotfast1Dmap = nplot;
                                    PlotVar.Plots{nplotfast1Dmap} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                end                                
%                                 PlotVar.Plots{nplotfast1Dmap}.PlotVector(x, map, [], [], true, 'Color', 'k', 'Ylim', [0 1.1*max(map)+1*(max(map)==0)]);

                                PlotVar.Plots{nplotfast1Dmap}.PlotVector(x, map, [], [], true, 'Ylim', [0 1.1*maxfast1Dmap+1*(maxfast1Dmap==0)]);
                                
                            case '1Dmap'
                                if isprop(EXP, 'Spkmaps_trajPercent')
                                    nplot = nplot + 1;
                                    contidx = numel(EXP.Nav.contrastValues) + 1;%PlotVar.ChosenContrast(:,c);
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                    for cc = 1:size(PlotVar.ChosenContrast,1)
                                        for gg = 1:size(PlotVar.ChosenGain,1)
                                            for rr = 1:size(PlotVar.ChosenRoomlength,1)
                                                for oo = 1:size(PlotVar.ChosenOutcome,1)
                                                    if sum(isnan(EXP.Spkmaps_trajPercent.map{contidx,PlotVar.ChosenGain(gg,g),PlotVar.ChosenRoomlength(rr,r),PlotVar.ChosenOutcome(oo,o)}.model.tuning(icell).meanrespModel)) == 0
                                                        Plot1DmapSpks(PlotVar.Plots{nplot}, EXP.Spkmaps_trajPercent.map{contidx,PlotVar.ChosenGain(gg,g),PlotVar.ChosenRoomlength(rr,r),PlotVar.ChosenOutcome(oo,o)}, icell);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            case '2DmapXSpd'
                                if isprop(EXP, 'Spkmaps_trajPercent_smthBallSpd')
                                    cont = numel(EXP.Nav.contrastValues) + 1;
                                    map = EXP.Spkmaps_trajPercent_smthBallSpd.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
                                    vecmean = EXP.Spkmaps_trajPercent_smthBallSpd.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXpos;
                                    vecSE = EXP.Spkmaps_trajPercent_smthBallSpd.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXpos;
                                    
                                    x = 1:size(map,2);
                                    y = 1:size(map,1);
                                    
                                    nplot = nplot + 1;
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                    PlotVar.Plots{nplot}.palette = 'parula';%'hot';%
                                    if max(map(:)) > 0
                                        PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], true,'Clim',[0 max(map(:))]);
                                        PlotVar.Plots{nplot}.PlotVector(vecmean, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mod(vecmean+vecSE,size(map,2)), y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mod(vecmean-vecSE,size(map,2)), y, [], [], true);
                                    end
                                end
                            case '2DmapXTheta'
                                if isprop(EXP, 'Spkmaps_trajPercent_LFPphase2')
                                    cont = numel(EXP.Nav.contrastValues) + 1;
                                    map = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
                                    map = (map - mean(map(:)))/mean(map(:));
                                    map = circshift(map,round(size(map,1)/2));
                                    map = repmat(map,[2 1]);                                    
                                    vecmean = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXmax;                                    
                                    vecmean = circshift(vecmean,round(numel(vecmean)/2));
                                    vecmean = repmat(vecmean,[2 1]);
                                    vecSE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXmax;                                    
                                    vecSE = circshift(vecSE,round(numel(vecSE)/2));
                                    vecSE = repmat(vecSE,[2 1]);
                                    [~,imaxprec] = max(abs(vecmean-mean(vecmean)));
                                    try
                                        zCOM = (vecmean(imaxprec) - mean(vecmean))/vecSE(imaxprec);
                                    catch
                                        zCOM = NaN;
                                    end
                                    x = 1:size(map,2);
                                    
                                    y = linspace(0,720,size(map,1));
                                    
                                    nplot = nplot + 1;
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                    PlotVar.Plots{nplot}.palette = 'dense';%'parula';%'hot';%
                                    if max(map(:)) > 0
                                        PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], true,'Clim',[min(map(:)) max(map(:))]);
                                        PlotVar.Plots{nplot}.PlotVector(vecmean, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mod(vecmean+vecSE,max(x)), y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mod(vecmean-vecSE,max(x)), y, [], [], true);
                                    end
                                    slope = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelslopeXY;
                                    phi0 = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelphi0XY;
                                    rho = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelrhoXY;
                                    slopeSE = quantile(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelslopeXY,0.95);
                                    phi0SE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelphi0XY;
                                    rhoSE = quantile(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelrhoXY,0.95);
                                    nShf = size(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelrhoXY,2);
                                    rhoZpos = sum(rho > EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelrhoXY)/nShf;
                                    rhoZneg = sum(rho < EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelrhoXY)/nShf;
                                    rhoZ = min(rhoZpos,rhoZneg);
                                    slopeZpos = sum(slope > EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelslopeXY)/nShf;
                                    slopeZneg = sum(slope < EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelslopeXY)/nShf;
                                    slopeZ = min(slopeZpos,slopeZneg);
%                                     [map1d,x] = fast1Dmap(theta(tidxAllphs),spktrain(tidxAllphs),dy,samplerate,nbYbinsmth,false);
                                    map1d = mean(map,2);
                                    map_phs = map1d/sum(map1d)*sum(spktrain(tidxAllphs));
                                    text(1,1,['rho = ' num2str(rho) ';z=' num2str(rhoZ) '  slope = ' num2str(slope) ';z=' num2str(slopeZ) '  offset = ' num2str(phi0) ';SE=' num2str(phi0SE) '     zCOM = ' num2str(zCOM)]);
                                end
                            case '2DmapXTheta-Correlation'
                                if isprop(EXP, 'Spkmaps_trajPercent_LFPphase2')
                                    cont = numel(EXP.Nav.contrastValues) + 1;
                                    xcorrtheta = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelX;
                                    Xpos = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelXmax;
                                    XposSE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelXmax;
                                    XposAmp = (EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelXRefmaxsinAmp);
                                    XposAmpShf = 0;%EXP.Spkmaps2d.trajPercent_LFPphase2_Shf{PlotVar.ChosenContrast(1,c),PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelXRefmaxsinAmp;
                                    XposAmpSE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelXmaxsinAmp;
                                    XposAmpShfSE = quantile(abs(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).corrModelXmaxsinAmp),0.95);
                                    XposOffset = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelXRefmaxOffset;
                                    XposOffsetSE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelXmaxOffset;
                                    
                                    zAmp = XposAmp/XposAmpSE;
%                                     zAmpShf = XposAmp/XposAmpShfSE;%
                                    zAmpShf = sum(abs(XposAmp) < abs(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).corrModelXRefmaxsinAmp))/numel(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).corrModelXRefmaxsinAmp);%(XposAmp-XposAmpShf)/XposAmpShfSE;
                                    zAmpShf2 = sum(abs(EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelXRefmaxAmp) < abs(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).corrModelXRefmaxAmp))/numel(EXP.Spkmaps_trajPercent_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).corrModelXRefmaxsinAmp);%(XposAmp-XposAmpShf)/XposAmpShfSE;
                                    zAmpShf = zAmpShf;%min(zAmpShf,zAmpShf2);
                                    x = 1:size(xcorrtheta,2);
                                    y = linspace(0,360,size(xcorrtheta,1));
                                    
                                    nplot = nplot + 1;
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g}); 
                                    PlotVar.Plots{nplot}.palette = 'parula';%'hot';%
                                    if max(map(:)) > 0
                                        PlotVar.Plots{nplot}.PlotMatrix(x, (y), xcorrtheta, [], true,'Clim',[0.5 1]);
%                                         PlotVar.Plots{nplot}.PlotMatrix(x, (y), map, [], true,'Clim',[0 max(map(:))]);
                                        PlotVar.Plots{nplot}.PlotVector(Xpos, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(mean(Xpos)*ones(1,numel(y)), y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(Xpos+XposSE, y, [], [], true);
                                        PlotVar.Plots{nplot}.PlotVector(Xpos-XposSE, y, [], [], true,'Xlim',[40 60]);
                                        
%                                         PlotVar.Plots{nplot}.PlotVector(XposShf, y, [], [], true);
%                                         PlotVar.Plots{nplot}.PlotVector(XposShf+XposShfSE, y, [], [], true);
%                                         PlotVar.Plots{nplot}.PlotVector(XposShf-XposShfSE, y, [], [], true);
                                    end
                                    slope = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelslopeXY;
                                    phi0 = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelphi0XY;
                                    rho = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meancorrModelrhoXY;
                                    slopeSE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelslopeXY;
                                    phi0SE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelphi0XY;
                                    rhoSE = EXP.Spkmaps_trajPercent_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SEcorrModelrhoXY;
                                    
                                    text(50,1,['Phsoffset = ' num2str(XposOffset) ';PhsoffsetSE=' num2str(XposOffsetSE) '  slope = ' num2str(slope) ';z=' num2str(slope/slopeSE) '  offset = ' num2str(phi0) ';SE=' num2str(phi0SE) '     zAmp = ' num2str(zAmp) '   zAmpshf = '  num2str(zAmpShf)]);
                                end
                            case 'Theta Mod.'
                                if isprop(EXP, 'Spkmaps_LFPphase2')
                                    cont = numel(EXP.SubsetVal.contrast) + 1;
                                    vecmean = EXP.Spkmaps_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModel;
                                    vecSE = EXP.Spkmaps_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModel;
                                    nplot = nplot + 1;
                                    vecmean = [vecmean vecmean(1)];
                                    vecSE = [vecSE vecSE(1)];
                                    x = linspace(0,360,numel(vecmean));
                                    PlotVar.Plots{nplot} = TDataPlot(Layout.subwindow{page, win, ((cell-1)*size(PlotVar.ChosenContrast,2) + (c-1))*size(PlotVar.ChosenObj,1) + k, (o-1)*size(PlotVar.ChosenGain,2) + g});
                                    PlotVar.Plots{nplot}.PlotVector(x, vecmean, vecSE, [], true);
                                    
                                    Amp = EXP.Spkmaps_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXsinAmp;
                                    AmpShf = EXP.Spkmaps_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXsinAmp;
                                    AmpSE = EXP.Spkmaps_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXsinAmp;
                                    Offset = EXP.Spkmaps_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).meanrespModelXsinOffset;
                                    OffsetSE = EXP.Spkmaps_LFPphase2.map{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).SErespModelXsinOffset;
                                    zAmp = (Amp-AmpShf)/AmpSE;
                                    nShf = numel(EXP.Spkmaps_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelXsinAmp);
                                    zAmpShf = sum(abs(Amp) <= abs(EXP.Spkmaps_LFPphase2.map_Shf{cont,PlotVar.ChosenGain(1,g),PlotVar.ChosenRoomlength(1,r),PlotVar.ChosenOutcome(1,o)}.model.tuning(icell).respModelXsinAmp))/nShf;
                                    text(1,mean(vecmean),['zAmp = ' num2str(zAmp) '   zAmpShf = ' num2str(zAmpShf) '   PhsOffset = ' num2str(Offset) '   PhsOffsetSE = ' num2str(OffsetSE)]);
                                end
                        end
                    end
                end
            end
        end
    end
end
end