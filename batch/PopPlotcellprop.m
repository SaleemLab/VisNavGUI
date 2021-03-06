function PopPlotcellprop(cellprop)
batch2p = true;
strlistvarname = {'2p data','electrophys data'};
[varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'dataset', 'SelectionMode', 'single', 'InitialValue', 1);
if ok && varnamesel == 1
    batch2p = true;
    strlistvarname = {'V1medial','V1lateral','PPC', 'AL'};
    [varnamesel,ok] = listdlg('ListString',strlistvarname, 'Name', 'area','InitialValue', 1);
    Cellpos2p = varnamesel;
    probestr{1} = cell2mat(strlistvarname(varnamesel));
    probestr{2} = 'none';
    nProbe = 1;
elseif ok
    batch2p = false;
    Cellpos2p = 0;
    if ~isfield(cellprop, 'Cellpos2p')
        cellprop.Cellpos2p = zeros(size(cellprop.Probe));
    end
    probestr{1} = 'CA1';
    probestr{2} = 'V1';
    nProbe = 2;
end

% cellprop.field2dXcorrthetamax = cellprop.field2dXcorrthetamax_toref;
% cellprop.field2dXRefcorrtheta = cellprop.field2dXcorrtheta_toref;


c{1} = 'c';
c{2} = 'k';
c{3} = 'm';
titlestr{1} = 'low';
titlestr{2} = 'med';
titlestr{3} = 'high';

field_orig = cellprop.field;

SSImin = -inf;%0.75;%
Xcorrmin = 0.75;%0.75;%
maxshift = 14;%+inf;%20;%30;%30;%30;%
Reliabilitymin = -inf;%0.75;%
thetaReliabilitymin = -inf;%0.5;%0.5;%
zth = 1.95;%1.95;%2.33;%
zth_lowhigh = 1.95;%1.95;%2.33;%-inf;%
Pth = 0.0;+inf;%
Pth_lowhigh = 0.0;%+inf;%0.05;%+inf;%
Pthgain = 0.05;
thetaPth = 0.05;
thetaPth_lowhigh = 0.05;

for iprobe = 1:nProbe
    if iprobe == 1
        chfocus = 1:32;
        xposfocus = 1:10;
        if batch2p
            maxminrate = +inf;%5;%
            minmaxrate = -inf;%3;%
            cellprop.bestchan = ones(size(cellprop.bestchan));
        else
            maxminrate = +inf;%8;%+inf;%10;%5;%
            minmaxrate = -inf;%1;%
            minspkwidth = 18;
        end
    else
        chfocus = 1:32;
        xposfocus = 1:10;
        if batch2p
            maxminrate = +inf;%5;%
            minmaxrate = -inf;%3;%cellprop.bestchan = ones(size(cellprop.bestchan));
        else
            maxminrate = inf;%+inf;%10;%5;%
            minmaxrate = -inf;%0.5;%
            minspkwidth = 18;
        end
    end
    figure('name',['Single cells : ' probestr{iprobe}]);
    iplot = 0;
    numBins = size(cellprop.field{2},2);
    
    cellprop.field = field_orig;
    for g = [2 1 3]
        shift{g} = (cellprop.fieldXcorrMax{g}-cellprop.fieldXcorrMax{2})';%numBins/(2*pi)*circ_dist(2*pi/numBins*cellprop.fieldPos{g},2*pi/numBins*cellprop.fieldPos{2})';%
        shift_early{g} = (cellprop.fieldXcorrMax_early{g}-cellprop.fieldXcorrMax_early{2})';
        shift_late{g} = (cellprop.fieldXcorrMax_late{g}-cellprop.fieldXcorrMax_late{2})';
        fieldZ{g} = zeros(1,size(cellprop.fieldSE{2},1));
        for icell = 1:size(cellprop.fieldSE{2},1)
            fieldZ{g}(icell) = (cellprop.field{g}(icell,min(floor(cellprop.fieldPos{g}(icell))+1,size(cellprop.field{g},2)))-mean(cellprop.field{g}(icell,:)))/cellprop.fieldSE{g}(icell,min(floor(cellprop.fieldPos{g}(icell))+1,size(cellprop.field{g},2)));
        end
    end
    
    cellprop.SSI{2} = abs(cellprop.fieldAmp{2})./mean(cellprop.field{2},2)';
    for g = [2 1 3]
        cellprop.phsmodulation{g} = cellprop.phsfieldsinAmp{g}./(mean(cellprop.phsfield{g},2)');
    end
    if iprobe == 1
        Finterneuron = (cellprop.min2maxSpkwf<minspkwidth & cellprop.Spatialmodulation <= 2);
    else
        Finterneuron = cellprop.min2maxSpkwf<minspkwidth;
    end
    
    numBinsX = size(cellprop.field{2},2);
    for g = [2 1 3]
%         if g == 2
            goodcells{g} = (cellprop.Goodcluster & ~isnan(cellprop.fieldPos{g}) & cellprop.Probe == iprobe & ismember(cellprop.bestchan,chfocus)... % & ismember(cellprop.XposPop,xposfocus)...
            & cellprop.SSI{2} >=SSImin & cellprop.fieldMin{2} <= maxminrate & cellprop.fieldMax{2} >= minmaxrate...
            & cellprop.reliabilityCorr{2}' >= Reliabilitymin...
            & abs(shift{g}')<maxshift & ismember(cellprop.Cellpos2p,Cellpos2p) & cellprop.fieldShfAmpZ{2} <= Pth & cellprop.fieldShfAmpZ{g} <= Pth_lowhigh);%& cellprop.fieldShfAmpZ{2} >= zth/2 & cellprop.fieldShfAmpZ{g} >= zth_lowhigh/2);
%         else
%             goodcells{g} = (cellprop.Goodcluster & ~isnan(cellprop.fieldPos{g}) & cellprop.Probe == iprobe & ismember(cellprop.bestchan,chfocus)... % & ismember(cellprop.XposPop,xposfocus)...
%             & cellprop.SSI{2} >=SSImin & cellprop.fieldMin{2} <= maxminrate & cellprop.fieldMax{2} >= minmaxrate...
%             & cellprop.reliabilityCorr{2}' >= Reliabilitymin...
%             & abs(shift{g}')<maxshift & ismember(cellprop.Cellpos2p,Cellpos2p) & cellprop.fieldShfAmpZ{2} <= Pth & cellprop.fieldShfAmpZ{g} <= Pth_lowhigh & cellprop.fieldShfXgainShiftZ{g} <= Pthgain);%& cellprop.fieldShfAmpZ{2} >= zth/2 & cellprop.fieldShfAmpZ{g} >= zth_lowhigh/2);
%         end
        if iprobe == 1
            goodcells{g} = goodcells{g} & ~Finterneuron;
        end
        reliablecells{g} = max(cellprop.fieldXcorr{g}(:,floor(numBinsX/4):floor(3*numBinsX/4)),[],2)'>Xcorrmin;
    end
    
    allcells = (cellprop.Goodcluster & ~cellprop.Finterneuron & cellprop.Probe == iprobe & ismember(cellprop.bestchan,chfocus));
        
    disp(['probe ' num2str(iprobe) ': ' num2str(num2str(100*sum(goodcells{2})/sum(allcells))) '%, ' num2str(sum(goodcells{2})) ' out of ' num2str(sum(allcells))])
    disp(['Pyramidal: ' num2str(num2str(100*sum(goodcells{2} & ~Finterneuron)/sum(goodcells{2}))) '% ; Interneuron: ' num2str(num2str(100*sum(goodcells{2} & Finterneuron)/sum(goodcells{2})))])
    disp(' ')
    maxpos = cellprop.fieldPos{2};
    
%     maxpos = [];
%     for icell = 1:size(cellprop.fieldunwrap{2},1)
%         xorig = 1:size(cellprop.fieldunwrap{2},2);
%         xinterp = 0:0.1:size(cellprop.fieldunwrap{2},2);
%         mapinterp = cellprop.fieldunwrap{2}(icell,:);%interp1(xorig,cellprop.fieldunwrap_set2{2}(icell,:),xinterp,'spline');%
%         [~,imax] = max(mapinterp);
%         maxpos(icell) = imax;%xinterp(imax);
%     end
%     cellprop.field = cellprop.fieldunwrap;
    
    [fPos,isort] = sort(maxpos,'ascend');
    
    for g = [2 1 3]
        iplot = iplot + 1;
        subplot(3,9,[(g-1)*3+1 (g-1)*3+2 (g-1)*3+1+9 (g-1)*3+2+9]);
%         normfields = (cellprop.field{g}(isort(goodcells{g}(isort)),:) - repmat(min(cellprop.field{g}(isort(goodcells{g}(isort)),:),[],2),[1 size(cellprop.field{g},2)]))./repmat(max(cellprop.field{g}(isort(goodcells{g}(isort)),:),[],2)-min(cellprop.field{g}(isort(goodcells{g}(isort)),:),[],2),[1 size(cellprop.field{g},2)]);
        if g == 2
            normfields = (cellprop.field{g}(isort(goodcells{g}(isort)),:) - repmat(mean(cellprop.field{g}(isort(goodcells{g}(isort)),:),2),[1 size(cellprop.field{g},2)]))./repmat(mean(cellprop.field{g}(isort(goodcells{g}(isort)),:),2),[1 size(cellprop.field{g},2)]);
        else
            normfieldsSigni = (cellprop.field{g}(isort(goodcells{g}(isort) & cellprop.fieldShfXgainShiftZ{g}(isort) <= Pthgain),:) - repmat(mean(cellprop.field{g}(isort(goodcells{g}(isort) & cellprop.fieldShfXgainShiftZ{g}(isort) <= Pthgain),:),2),[1 size(cellprop.field{g},2)]))./repmat(mean(cellprop.field{g}(isort(goodcells{g}(isort)  & cellprop.fieldShfXgainShiftZ{g}(isort) <= Pthgain),:),2),[1 size(cellprop.field{g},2)]);
            normfieldsNS = (cellprop.field{g}(isort(goodcells{g}(isort) & cellprop.fieldShfXgainShiftZ{g}(isort) > Pthgain),:) - repmat(mean(cellprop.field{g}(isort(goodcells{g}(isort) & cellprop.fieldShfXgainShiftZ{g}(isort) > Pthgain),:),2),[1 size(cellprop.field{g},2)]))./repmat(mean(cellprop.field{g}(isort(goodcells{g}(isort)  & cellprop.fieldShfXgainShiftZ{g}(isort) > Pthgain),:),2),[1 size(cellprop.field{g},2)]);
            normfields = [normfieldsNS;NaN(10,size(normfieldsSigni,2));normfieldsSigni];
        end
        
%         normfields = (cellprop.field{g}(isort(goodcells{g}(isort)),:))./repmat(max(cellprop.field{g}(isort(goodcells{g}(isort)),:),[],2),[1 size(cellprop.field{g},2)]);
        centeredfields{g} = zeros(size(cellprop.field{g}));
        for icell = 1:size(cellprop.field{g},1)
            if ~isnan(cellprop.fieldPos{g}(icell))
                centeredfields{g}(icell,:) = circshift(cellprop.field{g}(icell,:),-round(cellprop.fieldPos{g}(icell))+floor(numBinsX/2));
            end
        end
        

        imagesc(normfields);%imagesc(-normfields);
        if g == 2
            hold on;plot(maxpos(isort(goodcells{g}(isort))),1:sum(goodcells{g}),'w');
        else
            hold on;plot(maxpos(isort(goodcells{g}(isort) & cellprop.fieldShfXgainShiftZ{g}(isort) > Pthgain)),1:sum(goodcells{g} & cellprop.fieldShfXgainShiftZ{g} > Pthgain),'w');
            hold on;plot(maxpos(isort(goodcells{g}(isort) & cellprop.fieldShfXgainShiftZ{g}(isort) <= Pthgain)),(sum(goodcells{g} & cellprop.fieldShfXgainShiftZ{g} > Pthgain)+11):(sum(goodcells{g})+10),'w');
        end
%         hold on;plot(mod(maxpos(isort(goodcells{g}(isort)))+size(normfields,2)/2,size(normfields,2)),1:sum(goodcells{g}),'r');
%         imagesc(-centeredfields);
%         hold on;plot([numBinsX/2 numBinsX/2],[1 numel(goodcells{g})],'w');
        cmap = cmocean('dense');
        colormap(cmap);%colormap(bone);%colormap(parula);%(RedWhiteBlue);
        
%         hold on;plot(mod(fPos+numBinsX/2,numBinsX),1:numel(goodcells{g}),'r');
        set(gca,'Ydir','normal','Clim', [-0.5 0.5],'PlotBoxAspectRatio', [2 3 1]);
        title([probestr{iprobe} ' ' titlestr{g}]);
        xlabel('position');
        ylabel('cell#');
        
%         normfields = (cellprop.field{g}(goodcells{g}(isort),:))./repmat(max(cellprop.field{g}(goodcells{g}(isort),:),[],2),[1 size(cellprop.field{g},2)]);
        
        meanfield = nanmean(centeredfields{g});
        
        subplot(3,9,[(g-1)*3+3]);
        hcount{iprobe,g} = histcounts(maxpos(goodcells{g}),0:5:100);
        bar(5:5:100,hcount{iprobe,g}/sum(hcount{iprobe,g}),'Facecolor','k','Facealpha',0.5);
%         if iprobe == 2
%             hold on;plot(5:5:100,hcount{1,g}/sum(hcount{1,g}),'--','Color',[0.5 0.5 0.5]);
%         end
        set(gca,'Xlim',[0 105],'Ylim',[0 0.1])
        
        subplot(3,9,[(g-1)*3+3+9]);
        nbin = 10;
        vec{iprobe,g} = zeros(nbin,numBinsX);
        cellave = cell(1,nbin);
        for xx = 1:nbin
            cellave{xx} = find(ismember(floor(maxpos(isort(goodcells{g}(isort)))),mod(((xx-1)*(numBinsX/nbin)-(numBinsX/nbin)/2):((xx-1)*(numBinsX/nbin)+(numBinsX/nbin)/2),numBinsX)));
            vec{iprobe,g}(xx,:) = mean(normfields(ismember(floor(maxpos(isort(goodcells{g}(isort)))),mod(((xx-1)*(numBinsX/nbin)-(numBinsX/nbin)/2):((xx-1)*(numBinsX/nbin)+(numBinsX/nbin)/2),numBinsX)),:));
            hold on;
            if g == 2 && iprobe == 1
                plot(vec{iprobe,g}(xx,:)+xx,'k')
            elseif g == 2 && iprobe == 2
                plot(vec{iprobe,g}(xx,:)+xx,'k')
                plot(vec{1,g}(xx,:)+xx,'--','Color',[0.5 0.5 0.5])
            else
                plot(vec{iprobe,g}(xx,:)+xx,c{g})
                plot(vec{iprobe,2}(xx,:)+xx,c{2})
            end
        end
        set(gca,'Xlim', [1 numBinsX],'Ylim', [1 nbin+1]);%,'PlotBoxAspectRatio', [1 1.5 1]);
        
        subplot(3,9,[(g-1)*3+1 (g-1)*3+2 (g-1)*3+1+9 (g-1)*3+2+9]);
        for xx = 1:nbin
            hold on;plot(numBinsX*ones(1,numel(cellave{xx})),cellave{xx},'linewidth',2);
        end
        
%         disp(num2str(max(meanfield(1),meanfield(end))/meanfield(numBinsX/2)));
        
%         subplot(3,3,9);
%         hold on;plot(-numBinsX/2:49,mean(centeredfields,1),c{g});
%         [~,imax] = max(mean(centeredfields,1));
%         hold on;plot([imax-numBinsX/2+1 imax-numBinsX/2+1],[0 1],c{g});
%         set(gca,'Xlim',[-20 20],'Ylim',[0.3 1]);
    end    
    
%     [~,shift_low] = max(cellprop.fieldXcorr{1}(goodcells{g},:),[],2);
%     shift_low = shift_low - (floor(size(cellprop.fieldXcorr{1},2)/2) + 1);
    shift_low = cellprop.fieldXgainShift{1}-cellprop.fieldXgainShift{2};%shift{1}';%cellprop.fieldXcorrMax{1}(goodcells)-cellprop.fieldXcorrMax{2}(goodcells);
    shift_lowZ = cellprop.fieldShfXgainShiftZ{1};%cellprop.fieldXcorrMaxSE{1};
%     [~,shift_high] = max(cellprop.fieldXcorr{3}(goodcells,:),[],2);
%     shift_high = shift_high - (floor(size(cellprop.fieldXcorr{3},2)/2) + 1);
    shift_high = cellprop.fieldXgainShift{3}-cellprop.fieldXgainShift{2};%shift{3}';%cellprop.fieldXcorrMax{3}(goodcells)-cellprop.fieldXcorrMax{2}(goodcells);
    shift_highZ = cellprop.fieldShfXgainShiftZ{3};%cellprop.fieldXcorrMaxSE{3};
    
    disp('% neurons with significant response')
    disp(['low gain (out of signi medium): ' num2str(100*numel(shift_low(goodcells{1}))/numel(shift_low(goodcells{2})))])
    disp(['high gain (out of signi medium): ' num2str(100*numel(shift_high(goodcells{3}))/numel(shift_high(goodcells{2})))])
    disp(' ')
    
    disp('% neurons with significant shift')
    disp(['low gain (out of signi low): ' num2str(100*sum(shift_lowZ(goodcells{1} & ~isnan(shift_low))<=Pthgain & reliablecells{1}(goodcells{1} & ~isnan(shift_low)))/numel(shift_lowZ(goodcells{1} & ~isnan(shift_low))))])
    disp(['high gain (out of signi high): ' num2str(100*sum(shift_highZ(goodcells{3} & ~isnan(shift_high))<=Pthgain & reliablecells{3}(goodcells{3} & ~isnan(shift_high)))/numel(shift_high(goodcells{3} & ~isnan(shift_high))))])
    disp(['low gain (out of signi medium): ' num2str(100*sum(shift_lowZ(goodcells{1} & ~isnan(shift_low))<=Pthgain & reliablecells{1}(goodcells{1} & ~isnan(shift_low)))/numel(shift_lowZ(goodcells{2} & ~isnan(shift_low))))])
    disp(['high gain (out of signi medium): ' num2str(100*sum(shift_highZ(goodcells{3} & ~isnan(shift_high))<=Pthgain & reliablecells{3}(goodcells{3} & ~isnan(shift_high)))/numel(shift_high(goodcells{2} & ~isnan(shift_high))))])
    disp(' ')
    
    %to plot fraction of signi shifted neurons per position
% n = [];
% for xx = 0:9
% n(xx+1) = 100*sum(maxpos(goodcells{1} & ~isnan(shift_low))>=xx*10 & maxpos(goodcells{1} & ~isnan(shift_low))<=(xx+1)*10 & shift_lowZ(goodcells{1} & ~isnan(shift_low))<=Pthgain & reliablecells{1}(goodcells{1} & ~isnan(shift_low)))/sum(maxpos>=xx*10 & maxpos<=(xx+1)*10 & goodcells{1} & ~isnan(shift_low));
% end
% figure;plot(n)
% n = [];
% for xx = 0:9
% n(xx+1) = 100*sum(maxpos(goodcells{3} & ~isnan(shift_high))>=xx*10 & maxpos(goodcells{3} & ~isnan(shift_high))<=(xx+1)*10 & shift_lowZ(goodcells{3} & ~isnan(shift_high))<=Pthgain & reliablecells{3}(goodcells{3} & ~isnan(shift_high)))/sum(maxpos>=xx*10 & maxpos<=(xx+1)*10 & goodcells{3} & ~isnan(shift_high));
% end
% figure;plot(n)

    disp('% Exc neurons with significant shift')
    disp(['low gain (out of signi low): ' num2str(100*sum(shift_lowZ(goodcells{1} & ~isnan(shift_low) & ~Finterneuron)<=Pthgain & reliablecells{1}(goodcells{1} & ~isnan(shift_low) & ~Finterneuron))/numel(shift_lowZ(goodcells{1} & ~isnan(shift_low) & ~Finterneuron)))])
    disp(['high gain (out of signi high): ' num2str(100*sum(shift_highZ(goodcells{3} & ~isnan(shift_high) & ~Finterneuron)<=Pthgain & reliablecells{3}(goodcells{3} & ~isnan(shift_high) & ~Finterneuron))/numel(shift_high(goodcells{3} & ~isnan(shift_high) & ~Finterneuron)))])
    disp(' ');
    
    disp('% Inh neurons with significant shift')
    disp(['low gain (out of signi low): ' num2str(100*sum(shift_lowZ(goodcells{1} & ~isnan(shift_low) & Finterneuron)<=Pthgain & reliablecells{1}(goodcells{1} & ~isnan(shift_low) & Finterneuron))/numel(shift_lowZ(goodcells{1} & ~isnan(shift_low) & Finterneuron)))])
    disp(['high gain (out of signi high): ' num2str(100*sum(shift_highZ(goodcells{3} & ~isnan(shift_high) & Finterneuron)<=Pthgain & reliablecells{3}(goodcells{3} & ~isnan(shift_high) & Finterneuron))/numel(shift_high(goodcells{3} & ~isnan(shift_high) & Finterneuron)))])
    disp(' ');
%     for k = 1:floor(numel(goodcells)/20)+1
%         figure;
%         n = 0;
%         for icell = (1+(k-1)*20):k*20
%             if icell <= numel(goodcells)
%                 n = n+1;
%                 subplot(4,5,n);
%                 for g = 1:3
%                     hold on;plot(cellprop.field{g}(goodcells(icell),:),c{g})
%                     axis tight;
%                     text(1,0.9,[num2str(cellprop.fieldPos{2}(goodcells(icell))) '   ' num2str(shift_low(icell)) '   ' num2str(shift_high(icell))])
%                 end
%             end
%         end
%     end

%     shift_low = cellprop.fieldCOM{1}(goodcells) - cellprop.fieldCOM{2}(goodcells);%imax_low;%
%     shift_low(shift_low > numBinsX/2) = shift_low(shift_low > numBinsX/2) - numBinsX;
%     shift_low(shift_low < -numBinsX/2) = shift_low(shift_low < -numBinsX/2) + numBinsX;
%     shift_high = cellprop.fieldCOM{3}(goodcells) - cellprop.fieldCOM{2}(goodcells);%imax_high;%
%     shift_high(shift_high > numBinsX/2) = shift_high(shift_high > numBinsX/2) - numBinsX;
%     shift_high(shift_high < -numBinsX/2) = shift_high(shift_high < -numBinsX/2) + numBinsX;
    fieldpos = cellprop.fieldPos{2};
    subplot(3,6,13);
    scatter(fieldpos(goodcells{1} & shift_low>=0 & shift_lowZ>Pthgain & reliablecells{1}),shift_low(goodcells{1} & shift_low>=0 & shift_lowZ>Pthgain & reliablecells{1}),'MarkerEdgeColor','none','MarkerFaceColor',c{1},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5,'SizeData',10);
    hold on;scatter(fieldpos(goodcells{3} & shift_high<=0 & shift_highZ>Pthgain & reliablecells{3}),shift_high(goodcells{3} & shift_high<=0 & shift_highZ>Pthgain & reliablecells{3}),'MarkerEdgeColor','none','MarkerFaceColor',c{3},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.3,'SizeData',10);
    
    fieldpos = cellprop.fieldPos{2};
    hold on;scatter(fieldpos(goodcells{1} & shift_low<0 & shift_lowZ>Pthgain & reliablecells{1}),shift_low(goodcells{1} & shift_low<0 & shift_lowZ>Pthgain & reliablecells{1}),'MarkerEdgeColor','none','MarkerFaceColor',c{1},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5,'SizeData',10);
    hold on;scatter(fieldpos(goodcells{1} & shift_lowZ<=Pthgain & reliablecells{1}),shift_low(goodcells{1} & shift_lowZ<=Pthgain  & reliablecells{1}),'MarkerEdgeColor',c{1},'MarkerEdgeAlpha',0.7,'MarkerFaceColor',c{1},'MarkerFaceAlpha',0.5,'SizeData',10);
    fieldpossigni = round(fieldpos(goodcells{1} & shift_lowZ<=Pthgain & reliablecells{1})/10)+1;
    shiftlowsigni = shift_low(goodcells{1} & shift_lowZ<=Pthgain & reliablecells{1});
    shiftlowave = NaN(1,10);
    for xx = 1:10
        shiftlowave(xx) = median(shiftlowsigni(fieldpossigni==xx & ~isnan(shiftlowsigni)));
    end
    hold on;plot(5:10:95,shiftlowave);
    hold on;plot([0 numBinsX],[0 0],'k');
    set(gca,'Ylim',[-20 20]);
    xlabel('position');
    ylabel('response shift');
%     subplot(3,6,15);
    fieldpos = cellprop.fieldPos{2};
    scatter(fieldpos(goodcells{3} & shift_high>0 & shift_highZ>Pthgain & reliablecells{3}),shift_high(goodcells{3} & shift_high>0 & shift_highZ>Pthgain & reliablecells{3}),'MarkerEdgeColor','none','MarkerFaceColor',c{3},'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.5,'SizeData',10);
    hold on;scatter(fieldpos(goodcells{3} & shift_highZ<=Pthgain & reliablecells{3}),shift_high(goodcells{3} & shift_highZ<=Pthgain & reliablecells{3}),'MarkerEdgeColor',c{3},'MarkerEdgeAlpha',0.7,'MarkerFaceColor',c{3},'MarkerFaceAlpha',0.5,'SizeData',10);
    fieldpossigni = round(fieldpos(goodcells{1} & shift_highZ<=Pthgain & reliablecells{3})/10)+1;
    shifthighsigni = shift_high(goodcells{1} & shift_highZ<=Pthgain & reliablecells{3});
    shifthighave = NaN(1,10);
    for xx = 1:10
        shifthighave(xx) = median(shifthighsigni(fieldpossigni==xx & ~isnan(shifthighsigni)));
    end
    hold on;plot(5:10:95,shifthighave);
    hold on;plot([0 numBinsX],[0 0],'k');
    set(gca,'Ylim',[-20 20]);
    xlabel('position');
    ylabel('response shift');
    
    
    subplot(3,6,14);
    hlow = histogram(shift_low(goodcells{1}),-numBinsX/2+0.5:numBinsX/2+0.5,'Visible','off');
    Vlow = hlow.Values;
    hhigh = histogram(shift_high(goodcells{3}),-numBinsX/2+0.5:numBinsX/2+0.5,'Visible','off');
    Vhigh = hhigh.Values;
    hlow = histogram(shift_low(goodcells{1} & shift_lowZ<=Pthgain & reliablecells{1}),-numBinsX/2+0.5:numBinsX/2+0.5,'Visible','off');
    Vlowsigni = hlow.Values;
    hhigh = histogram(shift_high(goodcells{3} & shift_highZ<=Pthgain & reliablecells{3}),-numBinsX/2+0.5:numBinsX/2+0.5,'Visible','off');
    Vhighsigni = hhigh.Values;
    xbins = hhigh.BinEdges(1:end-1)+(hhigh.BinEdges(2)-hhigh.BinEdges(1))/2;
    
    
    subplot(3,6,14);
    b = barh(xbins,Vlow);
    b(1).FaceColor = c{1};
    b(1).FaceAlpha = 0.5;
    b(1).LineStyle = 'none';
    hold on;plot([0 numBinsX],[0 0],'k');
    set(gca,'Ylim',[-20 20],'Xlim',[0 max(max(Vlow),max(Vhigh))]);
    xlabel('# of cells');
    ylabel('response shift');
%     subplot(3,6,16);
    b = barh(xbins,Vhigh);
    b(1).FaceColor = c{3};
    b(1).FaceAlpha = 0.5;
    b(1).LineStyle = 'none';
    
    
    b = barh(xbins,Vlowsigni);
    b(1).EdgeColor = c{1};
    b(1).FaceColor = 'none';
    b(1).FaceAlpha = 0.5;
    hold on;plot([0 numBinsX],[0 0],'k');
    set(gca,'Ylim',[-20 20],'Xlim',[0 max(max(Vlowsigni),max(Vhighsigni))]);
    xlabel('# of cells');
    ylabel('response shift');
%     subplot(3,6,16);
    b = barh(xbins,Vhighsigni);
    b(1).EdgeColor = c{3};
    b(1).FaceColor = 'none';
    b(1).FaceAlpha = 0.5;
    hold on;plot([0 numBinsX],[0 0],'k');
    set(gca,'Ylim',[-20 20],'Xlim',[0 max(max(Vlowsigni),max(Vhighsigni))]);
    xlabel('# of cells');
    ylabel('response shift');
    
    hold on;plot([0 numBinsX],[0 0],'k');
    set(gca,'Ylim',[-20 20],'Xlim',[0 max([max(Vlow) max(Vhigh) max(Vlowsigni) max(Vhighsigni)])]);
    xlabel('# of cells');
    ylabel('response shift');
    
    subplot(3,3,8);
    scatter(shift_low(goodcells{1} & goodcells{3}),shift_high(goodcells{1} & goodcells{3}),'MarkerEdgeColor','k','MarkerFaceColor','k','SizeData',10,'Marker','.');
    hold on;scatter(shift_low(goodcells{1} & shift_lowZ<Pthgain & goodcells{3} & shift_highZ<Pthgain),shift_high(goodcells{1} & shift_lowZ<Pthgain & goodcells{3} & shift_highZ<Pthgain),'MarkerEdgeColor','k','MarkerFaceColor','none','SizeData',20);
    set(gca,'Ylim',[-20 20],'Xlim',[-20 20]);
    xlabel('low gain shift');
    ylabel('high gain shift');
    
    subplot(3,3,9);
    animal_list = unique(cellprop.animal);
    meanshift_low = zeros(1,numel(animal_list));
    meanshift_high = zeros(1,numel(animal_list));
    semshift_low = zeros(1,numel(animal_list));
    semshift_high = zeros(1,numel(animal_list));
    for ianimal = 1:numel(animal_list)
        meanshift_low(ianimal) = mean(shift_low(goodcells{1} & cellprop.animal == animal_list(ianimal)));
        meanshift_high(ianimal) = mean(shift_high(goodcells{3} & cellprop.animal == animal_list(ianimal)));
        semshift_low(ianimal) = sem(shift_low(goodcells{1} & cellprop.animal == animal_list(ianimal)));
        semshift_high(ianimal) = sem(shift_high(goodcells{3} & cellprop.animal == animal_list(ianimal)));
    end
    b = bar(1:numel(animal_list),meanshift_low);
    b(1).FaceColor = c{1};
    b(1).LineStyle = 'none';
    set(gca,'Ylim',[-20 20],'Xlim',[0 numel(animal_list)+1]);
    hold on;
    errorbar(1:numel(animal_list),meanshift_low,semshift_low,'Color',c{1},'LineStyle','none');
    b = bar(1:numel(animal_list),meanshift_high);
    b(1).FaceColor = c{3};
    b(1).LineStyle = 'none';
    errorbar(1:numel(animal_list),meanshift_high,semshift_high,'Color',c{3},'LineStyle','none');
    set(gca,'Ylim',[-20 20],'Xlim',[0 numel(animal_list)+1]);
    xlabel('animal #');
    ylabel('mean response shift');
    
%     %to plot mean peak firing across cells
%     subplot(3,3,9);
%     for g = 1:3
%         meanrate(g) = mean(cellprop.fieldMax{g}(goodcells{g}));
%         semrate(g) = sem(cellprop.fieldMax{g}(goodcells{g}));
%     end
%     b = bar(1:3,meanrate);
%     b(1).FaceColor = 'k';
%     b(1).LineStyle = 'none';
%     set(gca,'Ylim',[0 10],'Xlim',[0 4]);
%     hold on;
%     errorbar(1:3,meanrate,semrate,'Color','k','LineStyle','none');

% to plot preference for half1 or 2 in unit of SD
% max1 = max(cellprop.field{2}(:,1:numBinsX/2),[],2);
% max2 = max(cellprop.field{2}(:,numBinsX/2+1:numBinsX),[],2);
% for icell = 1:size(cellprop.field{2},1)
% SEfield(icell) = cellprop.fieldSE{2}(icell,min(size(cellprop.field{2},2),floor(cellprop.fieldPos{2}(icell))+1));
% end
% SEfield = SEfield(:);
% figure;
% histogram((max2(goodcells{2})-max1(goodcells{2}))./SEfield(goodcells{2}),-5.1:0.2:5.1,'EdgeColor','k','FaceColor','k');
% hold on;histogram((max2(goodcells{2} & (mod(maxpos,numBinsX/2)<5 | mod(maxpos,numBinsX/2)>45))-max1(goodcells{2} & (mod(maxpos,numBinsX/2)<5 | mod(maxpos,numBinsX/2)>45)))./SEfield(goodcells{2} & (mod(maxpos,numBinsX/2)<5 | mod(maxpos,numBinsX/2)>45)),-5.1:0.2:5.1,'EdgeColor','g','FaceColor','g');
% set(gca,'Ylim',[-5 120],'Xlim',[-5 5])

disp('% neurons with significant rate change')
disp(['low gain (out of signi low): ' num2str(100*sum(cellprop.rategainShfdiffZ{1}(goodcells{1} & ~isnan(shift_low))<=Pth)/sum(goodcells{1} & ~isnan(shift_low)))])
disp(['high gain (out of signi high): ' num2str(100*sum(cellprop.rategainShfdiffZ{3}(goodcells{3} & ~isnan(shift_high))<=Pth)/sum(goodcells{3} & ~isnan(shift_high)))])
disp(' ')
figure;
subplot(1,2,1)
hgaindiff = histogram((cellprop.rategain{1}(goodcells{1})-cellprop.rategain{2}(goodcells{1}))./cellprop.rategain{2}(goodcells{1}),-1:0.05:1,'Visible','off');
Vgaindiff = hgaindiff.Values;
xbins = hgaindiff.BinEdges(1:end-1)+(hgaindiff.BinEdges(2)-hgaindiff.BinEdges(1))/2;
hgaindiffsigni = histogram((cellprop.rategain{1}(goodcells{1} & cellprop.rategainShfdiffZ{1}<=Pth)-cellprop.rategain{2}(goodcells{1} & cellprop.rategainShfdiffZ{1}<=Pth))./cellprop.rategain{2}(goodcells{1} & cellprop.rategainShfdiffZ{1}<=Pth),-1:0.05:1,'Visible','off');
Vgaindiffsigni = hgaindiffsigni.Values;
hold on;
b = bar(xbins,Vgaindiff);
b(1).EdgeColor = c{1};
b(1).FaceColor = c{1};
b(1).FaceAlpha = 0.25;
b(1).LineStyle = 'none';
b = bar(xbins,Vgaindiffsigni);
b(1).EdgeColor = c{1};
b(1).FaceColor = c{1};
b(1).FaceAlpha = 1;
set(gca,'Xlim',[-1.05 1.05],'Ylim',[0 max(Vgaindiff)]);
ylabel('# of cells');
xlabel('Relative change in response between low and medium gain');
title('Low gain')

subplot(1,2,2)
hgaindiff = histogram((cellprop.rategain{3}(goodcells{3})-cellprop.rategain{2}(goodcells{3}))./cellprop.rategain{2}(goodcells{3}),-1:0.05:1,'Visible','off');
Vgaindiff = hgaindiff.Values;
xbins = hgaindiff.BinEdges(1:end-1)+(hgaindiff.BinEdges(2)-hgaindiff.BinEdges(1))/2;
hgaindiffsigni = histogram((cellprop.rategain{3}(goodcells{3} & cellprop.rategainShfdiffZ{3}<=Pth)-cellprop.rategain{2}(goodcells{3} & cellprop.rategainShfdiffZ{3}<=Pth))./cellprop.rategain{2}(goodcells{3} & cellprop.rategainShfdiffZ{3}<=Pth),-1:0.05:1,'Visible','off');
Vgaindiffsigni = hgaindiffsigni.Values;
hold on;
b = bar(xbins,Vgaindiff);
b(1).EdgeColor = c{3};
b(1).FaceColor = c{3};
b(1).FaceAlpha = 0.25;
b(1).LineStyle = 'none';
b = bar(xbins,Vgaindiffsigni);
b(1).EdgeColor = c{3};
b(1).FaceColor = c{3};
b(1).FaceAlpha = 1;
set(gca,'Xlim',[-1.05 1.05],'Ylim',[0 max(Vgaindiff)]);
ylabel('# of cells');
xlabel('Relative change in response between low and medium gain');
title('High gain')


figure('Name','Change in response between grating and plaid');
subplot(1,2,1)
hhalfdiff = histogram(cellprop.fieldunwraphalfdiff{2}(goodcells{2})./mean(cellprop.fieldunwrap{2}(goodcells{2},:),2)',-1.025:0.05:1.025,'Visible','off');
Vhalfdiff = hhalfdiff.Values;
xbins = hhalfdiff.BinEdges(1:end-1)+(hhalfdiff.BinEdges(2)-hhalfdiff.BinEdges(1))/2;
hhalfdiffsigni = histogram(cellprop.fieldunwraphalfdiff{2}(goodcells{2} & cellprop.fieldShfunwraphalfdiffZ{2}<=Pth & abs(cellprop.fieldunwraphalfdiff{2}./cellprop.fieldunwraphalfdiffSE{2})>=2)./mean(cellprop.fieldunwrap{2}(goodcells{2} & cellprop.fieldShfunwraphalfdiffZ{2}<=Pth & abs(cellprop.fieldunwraphalfdiff{2}./cellprop.fieldunwraphalfdiffSE{2})>=2,:),2)',-1.025:0.05:1.025,'Visible','off');
Vhalfdiffsigni = hhalfdiffsigni.Values;
disp('% neurons with significant diff between plaid and grating')
disp(['all neurons: ' num2str(100*sum(goodcells{2} & cellprop.fieldShfunwraphalfdiffZ{2}<=Pth  & abs(cellprop.fieldunwraphalfdiff{2}./cellprop.fieldunwraphalfdiffSE{2})>=2)/sum(goodcells{2}))]);
hold on;
b = bar(xbins,Vhalfdiff);
b(1).EdgeColor = 'k';
b(1).FaceColor = [0.5 0.5 0.5];
b(1).FaceAlpha = 0.5;
b(1).LineStyle = 'none';
b = bar(xbins,Vhalfdiffsigni);
b(1).EdgeColor = 'k';
b(1).FaceColor = 'k';
b(1).FaceAlpha = 1;
set(gca,'Xlim',[-1.05 1.05],'Ylim',[0 max(Vhalfdiff)]);
ylabel('# of cells');
xlabel('Relative change in response between plaid and grating');
title('all neurons')

subplot(1,2,2)
hhalfdiff = histogram(cellprop.fieldunwraphalfdiff{2}(goodcells{2} & (maxpos<10 | maxpos>90))./mean(cellprop.fieldunwrap{2}(goodcells{2} & (maxpos<10 | maxpos>90),:),2)',-1.025:0.05:1.025,'Visible','off');
Vhalfdiff = hhalfdiff.Values;
xbins = hhalfdiff.BinEdges(1:end-1)+(hhalfdiff.BinEdges(2)-hhalfdiff.BinEdges(1))/2;
hhalfdiffsigni = histogram(cellprop.fieldunwraphalfdiff{2}(goodcells{2} & (maxpos<10 | maxpos>90) & cellprop.fieldShfunwraphalfdiffZ{2}<=Pth  & abs(cellprop.fieldunwraphalfdiff{2}./cellprop.fieldunwraphalfdiffSE{2})>=2)./mean(cellprop.fieldunwrap{2}(goodcells{2} & (maxpos<10 | maxpos>90) & cellprop.fieldShfunwraphalfdiffZ{2}<=Pth  & abs(cellprop.fieldunwraphalfdiff{2}./cellprop.fieldunwraphalfdiffSE{2})>=2,:),2)',-1.025:0.05:1.025,'Visible','off');
Vhalfdiffsigni = hhalfdiffsigni.Values;
disp(['L1 neurons: ' num2str(100*sum(goodcells{2}  & (maxpos<10 | maxpos>90) & cellprop.fieldShfunwraphalfdiffZ{2}<=Pth  & abs(cellprop.fieldunwraphalfdiff{2}./cellprop.fieldunwraphalfdiffSE{2})>=2)/sum(goodcells{2} & (maxpos<10 | maxpos>90))) ' (n = ' num2str(sum(goodcells{2} & (maxpos<10 | maxpos>90))) ' neurons)']);
disp(' ')
hold on;
b = bar(xbins,Vhalfdiff);
b(1).EdgeColor = 'k';
b(1).FaceColor = [0.5 0.5 0.5];
b(1).FaceAlpha = 0.5;
b(1).LineStyle = 'none';
b = bar(xbins,Vhalfdiffsigni);
b(1).EdgeColor = 'k';
b(1).FaceColor = 'k';
b(1).FaceAlpha = 1;
set(gca,'Xlim',[-1.05 1.05],'Ylim',[0 max(Vhalfdiff)]);
ylabel('# of cells');
xlabel('Relative change in response between plaid and grating');
title('neurons around L1')

goodshiftedcells{1} = shift_lowZ<=Pthgain;%shift_low./shift_low_SE<=-zth;%true(size(goodcells{2}));%
goodshiftedcells{3} = shift_highZ<=Pthgain;%shift_high./shift_high_SE>=zth;%true(size(goodcells{2}));%
goodshiftedcells{2} = true(size(goodcells{2}));
nphs = size(cellprop.field2dXthetapos{2},2);
for g = [2 1 3]
%     goodthetacells{g} = goodcells{g};
    goodthetaRhocells{g} = cellprop.field2dShfrhoXYZ{2} <= thetaPth & ~isnan(cellprop.field2dslopeXY{g});
    goodthetacells{g} = (max(diff([cellprop.field2dXcorrthetamax{g} cellprop.field2dXcorrthetamax{g}(:,1)],1,2),[],2)' <= 5 & (cellprop.ZAmpRefthetaShfpos{2}<=thetaPth | cellprop.ZsinAmpRefthetaShfpos{2}<=thetaPth) & ~isnan(cellprop.field2dXRefcorrthetamaxsinAmp{2}));% & cellprop.ZsinAmpthetaShfpos{g}<=thetaPth_lowhigh);%& cellprop.thetareliabilityCorr{2}' >=thetaReliabilitymin);%abs(360/(2*pi)*circ_dist(2*pi/360*cellprop.field2dXRefcorrthetamaxOffset{g},2*pi/360*cellprop.field2dXRefcorrthetamaxOffset{2}))<=90);%& (cellprop.ZAmpthetaShfpos{g}>=thetaPth_lowhigh/2);% & (cellprop.thetareliabilityCorr{2}' >=thetaReliabilitymin);%& cellprop.field2dXRefcorrthetamaxOffsetSE{g}<=45);%abs(circ_dist(2*pi/360*cellprop.field2dXRefcorrthetamaxOffset{2},2*pi/360*cellprop.field2dXRefcorrthetamaxOffset{g}))./(2*pi/360*cellprop.field2dXRefcorrthetamaxOffsetSE{g})<=zth); %&...
%         %goodshiftedcells{g};% & abs(cellprop.field2dslopeXY{g}./cellprop.field2dslopeXYSE{g})>=zth_lowhigh;% & sum(cellprop.fieldXcorrtheta{1}>Xcorrmin,2)'==size(cellprop.fieldXcorrtheta{1},2)  & sum(cellprop.fieldXcorrtheta{3}>Xcorrmin,2)'==size(cellprop.fieldXcorrtheta{3},2));% & cellprop.thetareliabilityCorr{2}' >=thetaReliabilitymin);
    goodthetaphscells{g} = (cellprop.phsfieldShfsinAmpZ{2}<=thetaPth | cellprop.phsfieldShfsinAmpZ{2}<=thetaPth) & ~isnan(cellprop.phsfieldPos{2});%(cellprop.phsfieldShfAmpZ{2}<=thetaPth | cellprop.phsfieldShfsinAmpZ{2}<=thetaPth);% & (cellprop.phsfieldShfAmpZ{g}>=thetaPth_lowhigh/2);
    PvalTheta{g} = min(cellprop.ZAmpRefthetaShfpos{2},cellprop.ZsinAmpRefthetaShfpos{2});
    PvalTheta{g}(PvalTheta{g}==0) = 1/500;
    PvalThetaphs{g} = min(cellprop.phsfieldShfsinAmpZ{2},cellprop.phsfieldShfsinAmpZ{2});
    PvalThetaphs{g}(PvalThetaphs{g}==0) = 1/500;
end

pFishertheta = getFisherpvalue(PvalTheta{2}(goodcells{2}));
disp(['Theta precession: Fisher p-value = ' num2str(pFishertheta)])
disp(' ');
pFisherthetaphs = getFisherpvalue(PvalThetaphs{2}(goodcells{2}));
disp(['Theta modulation: Fisher p-value = ' num2str(pFisherthetaphs)])
disp(' ');

% goodthetacells = goodthetacells & goodthetaphscells;

% cellprop.field2dXRefcorrthetamaxOffset = cellprop.field2dShfphi0XY;
% cellprop.field2dXRefcorrthetamaxsinAmp = cellprop.field2drhoXY;

crossthetaprec = cellprop.field2dXRefcorrthetamaxOffset{2};%mean(cellprop.field2dXthetaposNorm{2}(:,3:7),2)-mean(cellprop.field2dXthetaposNorm{2}(:,12:16),2);
[~,thetasort] = sort(crossthetaprec,'ascend');
% thetasort = isort;
numBinsY = size(cellprop.field2dXPhstheta{g},2);

if iprobe == 1
    meanoffsetCA1 = 360/(2*pi)*mod(circ_mean(2*pi/360*cellprop.field2dXRefcorrthetamaxOffset{2}(goodcells{2} & goodthetacells{2} & ~isnan(cellprop.field2dXRefcorrthetamaxOffset{2}))'),2*pi);
else
    meanoffsetV1 = 360/(2*pi)*mod(circ_mean(2*pi/360*cellprop.field2dXRefcorrthetamaxOffset{2}(goodcells{2} & goodthetacells{2} & ~isnan(cellprop.field2dXRefcorrthetamaxOffset{2}))'),2*pi);
end


for g = [2 1 3]
    thetafield{g} = zeros(size(cellprop.field2dXPhstheta{g}));
    thetafield_centered{g} = zeros(size(cellprop.field2dXPhstheta{g}));
    field2dXcorrthetamax_centered{g} = zeros(size(cellprop.field2dXcorrthetamax{g}));
    for icell = 1:size(cellprop.field2dXPhstheta{g},1)
        if ~isnan(floor(cellprop.fieldCOM{2}(icell)))
            fieldtheta = squeeze(cellprop.field2dXPhstheta{g}(icell,:,:));
            fieldthetaref = squeeze(cellprop.field2dXPhstheta{g}(icell,:,:));
            if ~isnan(cellprop.fieldPos{g}(icell))
%                 thetafield{g}(icell,:,:) = circshift((fieldtheta-min(fieldthetaref(:)))/(max(fieldthetaref(:))-min(fieldthetaref(:))),numBinsX/2 - floor(cellprop.fieldPos{g}(icell)),2);
                thetafield{g}(icell,:,:) = circshift((fieldtheta-mean(cellprop.field{g}(icell,:)))/mean(cellprop.field{g}(icell,:)),numBinsX/2 - floor(cellprop.fieldPos{g}(icell)),2);
                fieldtheta_centered = squeeze(thetafield{g}(icell,:,:));
                %thetafield{g}(icell,:,:) = circshift((fieldtheta-min(fieldthetaref(:)))/(max(fieldthetaref(:))-min(fieldthetaref(:))),numBinsX/2 - floor(cellprop.fieldCOM{2}(icell)),2);
                phsoffset = cellprop.field2dXRefcorrthetamaxOffset{g}(icell);
                if ~isnan(phsoffset)
                    if iprobe == 1
                        thetafield_centered{g}(icell,:,:) = circshift((fieldtheta_centered-min(fieldtheta_centered(:)))/(max(fieldtheta_centered(:))-min(fieldtheta_centered(:))),round(meanoffsetCA1/360*numBinsY) - floor(phsoffset/360*numBinsY),1);
                        field2dXcorrthetamax_centered{g}(icell,:) = circshift(cellprop.field2dXcorrthetamax{g}(icell,:),round(meanoffsetCA1/360*numBinsY) - floor(phsoffset/360*numBinsY));
                    else
                        thetafield_centered{g}(icell,:,:) = circshift((fieldtheta_centered-min(fieldtheta_centered(:)))/(max(fieldtheta_centered(:))-min(fieldtheta_centered(:))),round(meanoffsetV1/360*numBinsY)-floor(phsoffset/360*numBinsY),1);
                        field2dXcorrthetamax_centered{g}(icell,:) = circshift(cellprop.field2dXcorrthetamax{g}(icell,:),round(meanoffsetV1/360*numBinsY)-floor(phsoffset/360*numBinsY));
                    end
                end
            end
        end
    end
    thetafieldcorr{g} = zeros(size(cellprop.field2dXRefcorrtheta{g}));
    thetafieldcorr_centered{g} = zeros(size(cellprop.field2dXRefcorrtheta{g}));
    for icell = 1:size(cellprop.field2dXcorrtheta{g},1)
        if ~isnan(cellprop.field2dXRefcorrthetamaxOffset{g}(icell))
            fieldtheta = squeeze(cellprop.field2dXcorrtheta{g}(icell,:,:));
            thetafieldcorr{g}(icell,:,:) = fieldtheta;
            phsoffset = cellprop.field2dXRefcorrthetamaxOffset{g}(icell);
            if ~isnan(phsoffset)
                if iprobe == 1
                    thetafieldcorr_centered{g}(icell,:,:) = circshift(thetafieldcorr{g}(icell,:,:),round(meanoffsetCA1/360*numBinsY) - floor(phsoffset/360*numBinsY),2);
                else
                    thetafieldcorr_centered{g}(icell,:,:) = circshift(thetafieldcorr{g}(icell,:,:),round(meanoffsetV1/360*numBinsY) - floor(phsoffset/360*numBinsY),2);
                end
            end
        end
    end
end
% for g = [2 1 3]
% [cellprop.Xmaxthetapos{g},~] = max(cellprop.field2dXcorrthetamax_toref{g},[],2);
% [cellprop.Xminthetapos{g},~] = min(cellprop.field2dXcorrthetamax_toref{g},[],2);
% end
figure('Name','Theta precession1')
for g = [2 1 3]
    thetamod = (cellprop.field2dXcorrthetamax{g}(thetasort(goodcells{g}(thetasort)  & goodthetacells{g}(thetasort)),:)-repmat(mean(cellprop.field2dXcorrthetamax{g}(thetasort(goodcells{g}(thetasort) & goodthetacells{g}(thetasort)),:),2),[1 nphs]))';
    thetamodref = (cellprop.field2dXcorrthetamax{2}(thetasort(goodcells{g}(thetasort)  & goodthetacells{g}(thetasort)),:)-repmat(mean(cellprop.field2dXcorrthetamax{2}(thetasort(goodcells{g}(thetasort) & goodthetacells{g}(thetasort)),:),2),[1 nphs]))';
    thetamod_centered = (field2dXcorrthetamax_centered{g}(thetasort(goodcells{g}(thetasort)  & goodthetacells{g}(thetasort)),:)-50)';
    thetamodref_centered = (field2dXcorrthetamax_centered{2}(thetasort(goodcells{g}(thetasort)  & goodthetacells{g}(thetasort)),:)-50)';
    subplot(4,9,(g-1)*3+1:(g-1)*3+2);
    imagesc(thetamod);
    set(gca,'Clim',[-2 2],'Ydir','normal');
    colormap(RedWhiteBlue);
    subplot(5,9,(g-1)*3+3);
    plot(thetamod,1:nphs,'color',[0.5 0.5 0.5]);
    hold on;plot(nanmean(thetamod,2),1:nphs,'k');
    hold on;plot(nanmean(thetamodref,2),1:nphs,'k--');
    set(gca,'Ylim',[1 nphs]);
    subplot(5,9,9+((g-1)*3+1));
    h = histogram(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g} & goodthetacells{g} & ~Finterneuron),0:40:360,'Visible','off');%-45:90:405,'Visible','off');%-22.5:45:382.5,'Visible','off');%
    thetaprecPhase1 = h.Values./sum(goodcells{g} & ~Finterneuron);
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
%     thetaprecPhase1(1) = thetaprecPhase1(1)+thetaprecPhase1(end);
%     thetaprecPhase1(end) = thetaprecPhase1(1);
    pRayleigh = circ_rtest(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g} & goodthetacells{g} & ~Finterneuron)/360*2*pi);
    if g==2
        disp('% Exc. neurons with significant theta phase shift')
        disp([num2str(100*sum(goodcells{g} & goodthetacells{g} & ~Finterneuron)/sum(goodcells{g} & ~Finterneuron)) ' (p = ' num2str(pRayleigh) ')']);
    end
%     hold on;bar(thetaprecPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
    h = histogram(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g} & Finterneuron),0:40:360,'Visible','off');%-45:90:405,'Visible','off');%
    thetaprecPhase2 = h.Values./sum(goodcells{g} & Finterneuron);
%     thetaprecPhase2(1) = thetaprecPhase2(1)+thetaprecPhase2(end);
%     thetaprecPhase2(end) = thetaprecPhase2(1);
    if sum(goodcells{g}  & goodthetacells{g} & Finterneuron) > 0
        pRayleigh = circ_rtest(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g} & Finterneuron)/360*2*pi);
    else
        pRayleigh = NaN;
    end
    if g==2
        disp('% Inh. neurons with significant theta phase shift')
        disp([num2str(100*sum(goodcells{g}  & goodthetacells{g} & Finterneuron)/sum(goodcells{g} & Finterneuron)) ' (p = ' num2str(pRayleigh) ')']);
    end
    if sum(goodcells{g}  & goodthetacells{g} & Finterneuron) > 0
        hold on;bar(xbins,[thetaprecPhase2;thetaprecPhase1]', 2, 'EdgeColor','none');
    else
        hold on;bar(xbins,thetaprecPhase1', 'Facecolor', 'k', 'EdgeColor','none');
    end
    
    if g==2
        disp('% neurons with significant linear-circular correlation')
        disp(num2str(100*sum(goodcells{g} & goodthetaRhocells{g})/sum(goodcells{g})));
    end
%     hold on;bar(thetaprecPhase,'FaceColor','k', 'EdgeColor','none');
%     [h,mu,ul,ll] = circ_mtest(2*pi/360*(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g})-cellprop.field2dXRefcorrthetamaxOffset{2}(goodcells{g}  & goodthetacells{g})), 0, 0.05);

%     h = histogram(cellprop.PhsthetaMod{g}(goodcells{g}  & goodthetaphscells{g})-1,0:0.05:2);
%     thetaPhase = h.Values./sum(goodcells{g});
%     bar(0:0.05:1.95,thetaPhase,'k');
%     set(gca,'Xlim',[-0.05 2]);
    
    subplot(5,9,9+((g-1)*3+2));
    if g == 2
        scatter(cellprop.min2maxSpkwf(goodcells{g} & goodthetacells{g}),cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g}));
        set(gca,'Ylim',[0 360]);
%         scatter(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g}),cellprop.field2dXRefcorrthetamaxsinAmp{g}(goodcells{g}  & goodthetacells{g}))
%         set(gca,'Xlim',[0 360]);
    end
    
    subplot(5,9,9+((g-1)*3+3));
    h = histogram(cellprop.field2dXRefcorrthetamaxsinAmp{g}(goodcells{g} & goodthetacells{g} & ~Finterneuron),0:1:10,'Visible','off');
    thetaPhase1 = h.Values./sum(goodcells{g} & ~Finterneuron);
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
%     hold on;bar(xbins,thetaPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
    h = histogram(cellprop.field2dXRefcorrthetamaxsinAmp{g}(goodcells{g}  & goodthetacells{g} & Finterneuron),0:1:10,'Visible','off');
    thetaPhase2 = h.Values./sum(goodcells{g} & Finterneuron);
    if sum(goodcells{g}  & goodthetacells{g} & Finterneuron) > 0
        hold on;bar(xbins,[thetaPhase2;thetaPhase1]',2,'EdgeColor','none');
    else
        hold on;bar(xbins,thetaPhase1,'Facecolor', 'k', 'EdgeColor','none');
    end
    set(gca,'Xlim',[-1 11]);
    
    subplot(5,9,18+((g-1)*3+1));
    scatter(cellprop.field2dXRefcorrthetamaxOffset{2}(goodcells{g}  & goodthetacells{g}),cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g}))
    [h,p] = ttest(cellprop.field2dXRefcorrthetamaxOffset{2}(goodcells{g}  & goodthetacells{g}),cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g}));
    set(gca,'Xlim',[0 360],'Ylim',[0 360]);
    hold on;plot([0 360],[0 360],'k');
    
    subplot(4,9,18+((g-1)*3+2));
    scatter(abs(cellprop.field2dXRefcorrthetamaxsinAmp{2}(goodcells{g}  & goodthetacells{g})),abs(cellprop.field2dXRefcorrthetamaxsinAmp{g}(goodcells{g}  & goodthetacells{g})))
    set(gca,'Xlim',[0 20],'Ylim',[0 20]);
    hold on;plot([0 20],[0 20],'k');
    
    
    subplot(5,9,18+((g-1)*3+3));
    thetamax_ahead{g} = (cellprop.field2dXcorrthetamax{g}-numBinsX/2);
    thetamax_ahead{g}(thetamax_ahead{2}>=0) = 0;
    thetamax_behind{g} = (cellprop.field2dXcorrthetamax{g}-numBinsX/2);
    thetamax_behind{g}(thetamax_behind{2}<=0) = 0;
    early = numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(field2dXcorrthetamax_centered{g}(:,1:9)-50),[],2);%min(field2dXcorrthetamax_centered{g}-50,[],2);%%getCircularAverage(-thetamax_ahead{g}',0,0.1,0.05);%
    late = numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(field2dXcorrthetamax_centered{g}(:,10:18)-50),[],2);%max(field2dXcorrthetamax_centered{g}-50,[],2);%numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(thetamax_behind{g}),[],2);%getCircularAverage(thetamax_behind{g}',0,0.1,0.05);%
    early_ref = numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(field2dXcorrthetamax_centered{2}(:,1:9)-50),[],2);%min(field2dXcorrthetamax_centered{2}-50,[],2);%numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(thetamax_ahead{2}),[],2);%getCircularAverage(-thetamax_ahead{2}',0,0.1,0.05);%
    late_ref = numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(field2dXcorrthetamax_centered{2}(:,10:18)-50),[],2);%max(field2dXcorrthetamax_centered{2}-50,[],2);%numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(thetamax_behind{2}),[],2);%getCircularAverage(thetamax_behind{2}',0,0.1,0.05);%
    
    iearly = getCircularAverage(-(cellprop.field2dXcorrthetamax{g}-repmat(mean(cellprop.field2dXcorrthetamax{g},2),[1 size(cellprop.field2dXcorrthetamax{g},2)]))',0,0.1,0.05);%numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(thetamax_ahead{g}),[],2);%getCircularAverage(-thetamax_ahead{g}',0,0.01,0.05);%
    ilate = getCircularAverage((cellprop.field2dXcorrthetamax{g}-repmat(mean(cellprop.field2dXcorrthetamax{g},2),[1 size(cellprop.field2dXcorrthetamax{g},2)]))',0,0.1,0.05);%numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(thetamax_behind{g}),[],2);%getCircularAverage(thetamax_behind{g}',0,0.01,0.05);%
    iearly_ref = getCircularAverage(-(cellprop.field2dXcorrthetamax{2}-repmat(mean(cellprop.field2dXcorrthetamax{2},2),[1 size(cellprop.field2dXcorrthetamax{2},2)]))',0,0.1,0.05);%numBinsX/(2*pi)*circ_mean(2*pi/numBinsX*(thetamax_ahead{2}),[],2);%getCircularAverage(-thetamax_ahead{2}',0,0.01,0.05);%
    ilate_ref = getCircularAverage((cellprop.field2dXcorrthetamax{2}-repmat(mean(cellprop.field2dXcorrthetamax{2},2),[1 size(cellprop.field2dXcorrthetamax{2},2)]))',0,0.1,0.05);
        

%     earlyshift = early(goodcells{g}  & goodthetacells{g})-early_ref(goodcells{g}  & goodthetacells{g});
%     lateshift = late(goodcells{g}  & goodthetacells{g})-late_ref(goodcells{g}  & goodthetacells{g});
%     plot([earlyshift lateshift],[1.5+0.01*randn(sum(goodcells{g}  & goodthetacells{g}),1) 1.7+0.01*randn(sum(goodcells{g}  & goodthetacells{g}),1)],'o')
%     hold on;barh([1.5 1.7],median([earlyshift lateshift],1))
%     hold on;errorbar(median([earlyshift lateshift],1),[1.5 1.7],std([earlyshift lateshift],[],1)/sqrt(sum((goodcells{g}  & goodthetacells{g}))),'horizontal')
%     set(gca,'Ylim',[1.4 1.8],'Xlim',[-20 20]);
    
    histogram(shift_early{g}(goodcells{g}  & goodthetacells{g}),-numBinsX/2+0.5:numBinsX/2+0.5, 'EdgeColor','none');
    hold on; histogram(shift_late{g}(goodcells{g}  & goodthetacells{g}),-numBinsX/2+0.5:numBinsX/2+0.5, 'EdgeColor','none');
    set(gca,'Xlim',[-20 20]);
%     scatter(shift_early{g}(goodcells{g}),shift_late{g}(goodcells{g}));
%     set(gca,'Xlim',[-20 20],'Ylim',[-20 20]);
    
%     histogram(early(goodcells{g}  & goodthetacells{g})-early_ref(goodcells{g}  & goodthetacells{g}),-19:2:21)
%     hold on;histogram(late(goodcells{g}  & goodthetacells{g})-late_ref(goodcells{g}  & goodthetacells{g}),-19:2:21)
%     scatter(early(goodcells{g}  & goodthetacells{g})-early_ref(goodcells{g}  & goodthetacells{g}),late(goodcells{g}  & goodthetacells{g})-late_ref(goodcells{g}  & goodthetacells{g}))
%     hold on;scatter(late_ref(goodcells{g}  & goodthetacells{g}),late(goodcells{g}  & goodthetacells{g}))
    
    %the following line is to test if the difference holds true when we
    %normalize by the visual distance to the center
%     scatter(late_ref(goodcells{g}  & goodthetacells{g})./abs(cellprop.fieldsize_behind{2}(goodcells{g}  & goodthetacells{g})) - early_ref(goodcells{g}  & goodthetacells{g})./abs(cellprop.fieldsize_ahead{2}(goodcells{g}  & goodthetacells{g})),...
%             late(goodcells{g}  & goodthetacells{g})./abs(cellprop.fieldsize_behind{g}(goodcells{g}  & goodthetacells{g})) - early(goodcells{g}  & goodthetacells{g})./abs(cellprop.fieldsize_ahead{g}(goodcells{g}  & goodthetacells{g})))
    
%     scatter((late_ref(goodcells{g}  & goodthetacells{g})-early_ref(goodcells{g}  & goodthetacells{g}))./cellprop.fieldsize{2}(goodcells{g}  & goodthetacells{g})',(late(goodcells{g}  & goodthetacells{g})-early(goodcells{g}  & goodthetacells{g}))./cellprop.fieldsize{g}(goodcells{g}  & goodthetacells{g})')

%     scatter(cellprop.PhsXminthetapos{2}(goodcells{g}  & goodthetacells{g}),cellprop.PhsXminthetapos{g}(goodcells{g}  & goodthetacells{g}));
%     scatter(cellprop.PhsXaheadthetapos{2}(goodcells{g}  & goodthetacells{g}),cellprop.PhsXaheadthetapos{g}(goodcells{g}  & goodthetacells{g}));
    
%     hold on;
%     scatter(late_ref(goodcells{g}  & goodthetacells{g}),early_ref(goodcells{g}  & goodthetacells{g}))
%     set(gca,'Xlim',[-20 20],'Ylim',[-20 20]);
%     hold on;plot([-20 20],[-20 20],'k');
%     hold on;plot([0 0],[-20 20],'k');
%     hold on;plot([-20 20],[0 0],'k');
    
    subplot(5,9,27+((g-1)*3+1:(g-1)*3+2));
    scatter(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetacells{g} & Finterneuron)/18*360+90-90,360),cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g} & Finterneuron),'b')
    hold on;scatter(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetacells{g} & ~Finterneuron)/18*360+90-90,360),cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g} & ~Finterneuron),'r')
    set(gca,'Xlim',[0 360],'Ylim',[0 360]);
    if g == 2
        [rho_circcorr, p_circcorr] = circ_corrcc(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetacells{g})/18*360+90-90,360)/360*2*pi, cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g})/360*2*pi);
        disp(['correlation between precession and modulation phase: rho = ' num2str(rho_circcorr) ' pval = ' num2str(p_circcorr)]);
    end
    
    subplot(5,9,36+((g-1)*3+1));
    h = histogram(cellprop.field2drhoXY{g}(goodcells{g}),-0.3:0.01:0.3,'Visible','off');%-45:90:405,'Visible','off');%-22.5:45:382.5,'Visible','off');%
    lincirccorr = h.Values./sum(goodcells{g});
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
    hold on;bar(xbins,lincirccorr,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    h = histogram(cellprop.field2drhoXY{g}(goodcells{g} & goodthetaRhocells{g}),-0.3:0.01:0.3,'Visible','off');%-45:90:405,'Visible','off');%-22.5:45:382.5,'Visible','off');%
    lincirccorr = h.Values./sum(goodcells{g});
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
    hold on;bar(xbins,lincirccorr,'FaceColor','k','EdgeColor','none');
    set(gca,'Xlim',[-0.3 0.3]);
    subplot(5,9,36+((g-1)*3+2));
    h = histogram(cellprop.field2drhoXY{g}(goodcells{g}),-0.3:0.01:0.3,'Visible','off');%-45:90:405,'Visible','off');%-22.5:45:382.5,'Visible','off');%
    lincirccorr = h.Values./sum(goodcells{g});
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
    hold on;bar(xbins,lincirccorr,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    h = histogram(cellprop.field2drhoXY{g}(goodcells{g} & goodthetacells{g}),-0.3:0.01:0.3,'Visible','off');%-45:90:405,'Visible','off');%-22.5:45:382.5,'Visible','off');%
    lincirccorr = h.Values./sum(goodcells{g});
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
    hold on;bar(xbins,lincirccorr,'FaceColor','k','EdgeColor','none');
    set(gca,'Xlim',[-0.3 0.3]);
end

figure('Name','Theta Modulation')
for g = [2 1 3]
    subplot(2,9,(g-1)*3+1:(g-1)*3+2);
    [~,phssort] = sort(cellprop.phsfieldPos{g},'ascend');
    normphsfield = cellprop.phsfield{g}(phssort(goodcells{g}(phssort)  & goodthetaphscells{g}(phssort)),:);
    normphsfield = (normphsfield - repmat(mean(normphsfield,2),[1 size(normphsfield,2)]))./repmat(mean(normphsfield,2),[1 size(normphsfield,2)]);
    imagesc(normphsfield);
    cmap = cmocean('dense');
    colormap(cmap);
    set(gca,'Clim',[-0.2 0.2],'Ydir','normal')
    
    subplot(2,9,9+((g-1)*3+1));
    h = histogram(mod(cellprop.phsfieldCOM{g}(goodcells{g} & goodthetaphscells{g} & ~Finterneuron)/18*360+90-90,360),0:40:360,'Visible','off');%-45:90:405,'Visible','off');
    thetaPhase1 = h.Values./sum(goodcells{g} & ~Finterneuron);
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
%     thetaPhase1(1) = thetaPhase1(1)+thetaPhase1(end);
%     thetaPhase1(end) = thetaPhase1(1);
%     hold on;bar(thetaPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
    h = histogram(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetaphscells{g} & Finterneuron)/18*360+90-90,360),0:40:360,'Visible','off');%-45:90:405,'Visible','off');
    thetaPhase2 = h.Values./sum(goodcells{g} & Finterneuron);
%     thetaPhase2(1) = thetaPhase2(1)+thetaPhase2(end);
%     thetaPhase2(end) = thetaPhase2(1);
    if sum(goodcells{g} & goodthetaphscells{g} & Finterneuron) > 0
        hold on;b = bar(xbins,[thetaPhase2;thetaPhase1]', 2, 'EdgeColor','none');
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
    else
        hold on;b = bar(xbins,thetaPhase1, 'Facecolor', 'k', 'EdgeColor','none');
    end
%     hold on;bar(thetaPhase,'FaceColor','k', 'EdgeColor','none');
    if g==2
        disp('% neurons with significant theta modulation of firing rate')
        disp(num2str(100*sum(goodthetaphscells{2}(goodcells{2}))/sum(goodcells{2})))
        disp('% Exc. neurons with significant theta modulation of firing rate')
        pRayleigh = circ_rtest(cellprop.phsfieldCOM{g}(goodcells{g} & goodthetaphscells{g} & ~Finterneuron)/18*2*pi);
        disp([num2str(100*sum(goodthetaphscells{2}(goodcells{2} & ~Finterneuron))/sum(goodcells{2} & ~Finterneuron)) ' (p = ' num2str(pRayleigh) ')'])
        disp('% Inh. neurons with significant theta modulation of firing rate')
        if sum(goodcells{g} & goodthetaphscells{g} & Finterneuron) > 0
            pRayleigh = circ_rtest(cellprop.phsfieldCOM{g}(goodcells{g} & goodthetaphscells{g} & Finterneuron)/18*2*pi);
        else
            pRayleigh = NaN;
        end
        disp([num2str(100*sum(goodthetaphscells{2}(goodcells{2} & Finterneuron))/sum(goodcells{2} & Finterneuron)) ' (p = ' num2str(pRayleigh) ')'])
        disp(' ')
    end

    subplot(2,9,9+((g-1)*3+2));
    if g ~= 2
        scatter(mod(cellprop.phsfieldCOM{2}(goodcells{g} & goodthetaphscells{g})/18*360+90-90,360),mod(cellprop.phsfieldCOM{g}(goodcells{g} & goodthetaphscells{g})/18*360+90-90,360))
        set(gca,'Xlim',[0 360],'Ylim',[0 360]);
        hold on;plot([0 360],[0 360]);
    else
        scatter(cellprop.min2maxSpkwf(goodcells{g} & goodthetaphscells{g}),mod(cellprop.phsfieldCOM{g}(goodcells{g} & goodthetaphscells{g})/18*360+90-90,360))
        set(gca,'Ylim',[0 360]);
    end
%     [h,mu,ul,ll] = circ_mtest(2*pi/360*(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetacells{g} & goodthetaphscells{g}-90,360)-mod(cellprop.phsfieldPos{2}(goodcells{g}  & goodthetacells{g} & goodthetaphscells{g})-90,360)), 0, 0.05);
    
    subplot(2,9,9+((g-1)*3+3));
    h = histogram(cellprop.phsmodulation{g}(goodcells{g} & goodthetaphscells{g} & ~Finterneuron),0:0.05:0.7,'Visible','off');
    thetaPhaseMod1 = h.Values./sum(goodcells{g} & ~Finterneuron);
    xbins = h.BinEdges(1:end-1)+(h.BinEdges(2)-h.BinEdges(1))/2;
%     hold on;bar(xbins,thetaPhaseMod,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
    h = histogram(cellprop.phsmodulation{g}(goodcells{g} & goodthetaphscells{g} & Finterneuron),0:0.05:0.7,'Visible','off');
    thetaPhaseMod2 = h.Values./sum(goodcells{g} & Finterneuron);
    if sum(goodcells{g} & goodthetaphscells{g} & Finterneuron) > 0
        hold on;b = bar(xbins,[thetaPhaseMod2;thetaPhaseMod1]', 2, 'EdgeColor','none');
        set(gca,'Xlim',[-0.05 0.75]);
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
    else
        hold on;b = bar(xbins,thetaPhaseMod1, 'Facecolor', 'k', 'EdgeColor','none');
        set(gca,'Xlim',[-0.05 0.75]);
    end
    %ciplot(gca,[],nanmean(normthetaPhsfield,1),nanstd(normthetaPhsfield,[],1)/sqrt(size(normthetaPhsfield,1)),0.5,'k');
    
%     thetamax = max(cellprop.field2dXcorrthetamax{g}(goodcells{g}  & goodthetacells{g},:),[],2);
%     thetamin = min(cellprop.field2dXcorrthetamax{g}(goodcells{g}  & goodthetacells{g},:),[],2);
%     polarplot(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g}),thetamax-thetamin,'.')
end

%to plot fast-spiking and non-fastspiking distri of pref theta phase
% figure;
% g=2;
% subplot(1,2,1);
% h = histogram(mod(cellprop.phsfieldCOM{g}(goodcells{g} & ~Finterneuron)/18*360,360),-22.5:45:382.5,'Visible','off');
%     thetaPhase = h.Values./sum(goodcells{g} & ~Finterneuron);
%     thetaPhase(1) = thetaPhase(1)+thetaPhase(end);
%     thetaPhase(end) = thetaPhase(1);
%     hold on;bar(thetaPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
%     h = histogram(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetaphscells{g} & ~Finterneuron)/18*360,360),-22.5:45:382.5,'Visible','off');
%     thetaPhase = h.Values./sum(goodcells{g} & ~Finterneuron);
%     thetaPhase(1) = thetaPhase(1)+thetaPhase(end);
%     thetaPhase(end) = thetaPhase(1);
%     hold on;bar(thetaPhase,'FaceColor','k', 'EdgeColor','none');
% subplot(1,2,2);
% h = histogram(mod(cellprop.phsfieldCOM{g}(goodcells{g} & Finterneuron)/18*360,360),-22.5:45:382.5,'Visible','off');
%     thetaPhase = h.Values./sum(goodcells{g} & Finterneuron);
%     thetaPhase(1) = thetaPhase(1)+thetaPhase(end);
%     thetaPhase(end) = thetaPhase(1);
%     hold on;bar(thetaPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
%     h = histogram(mod(cellprop.phsfieldCOM{g}(goodcells{g}  & goodthetaphscells{g} & Finterneuron)/18*360,360),-22.5:45:382.5,'Visible','off');
%     thetaPhase = h.Values./sum(goodcells{g} & Finterneuron);
%     thetaPhase(1) = thetaPhase(1)+thetaPhase(end);
%     thetaPhase(end) = thetaPhase(1);
%     hold on;bar(thetaPhase,'FaceColor','k', 'EdgeColor','none');

%to plot fast-spiking and non-fastspiking distri of theta modulation
% figure;
% g=2;
% subplot(1,2,1);
% h = histogram(cellprop.phsmodulation{g}(goodcells{g} & cellprop.min2maxSpkwf>=19),0:0.05:4,'Visible','off');
%     thetaPhaseMod = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf>=19);
%     hold on;bar(0.025:0.05:3.975,thetaPhaseMod,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
%     h = histogram(cellprop.phsmodulation{g}(goodcells{g} & goodthetaphscells{g} & cellprop.min2maxSpkwf>=19),0:0.05:4,'Visible','off');
%     thetaPhaseMod = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf>=19);
%     hold on;bar(0.025:0.05:3.975,thetaPhaseMod,'FaceColor','k', 'EdgeColor','none');
% subplot(1,2,2);
% h = histogram(cellprop.phsmodulation{g}(goodcells{g} & cellprop.min2maxSpkwf<19),0:0.05:4,'Visible','off');
%     thetaPhaseMod = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf<19);
%     hold on;bar(0.025:0.05:3.975,thetaPhaseMod,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
%     h = histogram(cellprop.phsmodulation{g}(goodcells{g} & goodthetaphscells{g} & cellprop.min2maxSpkwf<19),0:0.05:4,'Visible','off');
%     thetaPhaseMod = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf<19);
%     hold on;bar(0.025:0.05:3.975,thetaPhaseMod,'FaceColor','k', 'EdgeColor','none');


%to plot fast-spiking and non-fastspiking distri of theta phase offsets
% figure;
% g=2;
% subplot(1,2,1);
% h = histogram(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g} & cellprop.min2maxSpkwf>=19),0:40:360,'Visible','off');
%     thetaprecPhase = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf>=19);
%     hold on;bar(thetaprecPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
%     h = histogram(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g} & cellprop.min2maxSpkwf>=19),0:40:360,'Visible','off');
%     thetaprecPhase = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf>=19);
%     hold on;bar(thetaprecPhase,'FaceColor','k', 'EdgeColor','none');
% subplot(1,2,2);
% h = histogram(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g} & cellprop.min2maxSpkwf<19),0:40:360,'Visible','off');
%     thetaprecPhase = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf<19);
%     hold on;bar(thetaprecPhase,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6, 'EdgeColor','none');
%     h = histogram(cellprop.field2dXRefcorrthetamaxOffset{g}(goodcells{g}  & goodthetacells{g} & cellprop.min2maxSpkwf<19),0:40:360,'Visible','off');
%     thetaprecPhase = h.Values./sum(goodcells{g} & cellprop.min2maxSpkwf<19);
%     hold on;bar(thetaprecPhase,'FaceColor','k', 'EdgeColor','none');

figure('Name','Theta precession2')
for g = [2 1 3]
    subplot(5,3,g);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells{g}  & goodthetacells{g} & ((crossthetaprec >= 0 & crossthetaprec < 45) | (crossthetaprec <= 360 & crossthetaprec > 315)),:,:),1));
    thetamapref = squeeze(nanmean(thetafield{2}(goodcells{g}  & goodthetacells{g} & ((crossthetaprec >= 0 & crossthetaprec < 45) | (crossthetaprec <= 360 & crossthetaprec > 315)),:,:),1));
    if g==2
        disp(sum(goodcells{g}  & goodthetacells{g} & ((crossthetaprec >= 0 & crossthetaprec < 45) | (crossthetaprec <= 360 & crossthetaprec > 315)))/sum(goodcells{g}));
    end
    thetamap = circshift(thetamap,size(thetamap,1)/2);
    thetamap = repmat(thetamap,[2 1]);
%     [x0,y0] = meshgrid(1:size(thetamap,2),1:size(thetamap,1));
%     [xinterp,yinterp] = meshgrid(1:0.5:size(thetamap,2),1:0.5:size(thetamap,1));
%     numBins = numel(1:0.5:size(thetamap,2));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:2*nphs,'w');
    if g ~= 2
        hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    end
    set(gca,'Clim',[0 0.8],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+3);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 45 & crossthetaprec < 135,:,:),1));
    thetamapref = squeeze(nanmean(thetafield{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 45 & crossthetaprec < 135,:,:),1));
    if g==2
        disp(sum(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 45 & crossthetaprec < 135)/sum(goodcells{g}));
    end
    thetamap = circshift(thetamap,size(thetamap,1)/2);
    thetamap = repmat(thetamap,[2 1]);
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:2*nphs,'w');
    if g ~= 2
        hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    end
    set(gca,'Clim',[0 0.8],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+6);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 135 & crossthetaprec < 225,:,:),1));
    thetamapref = squeeze(nanmean(thetafield{2}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 135 & crossthetaprec < 225,:,:),1));
    if g==2
        disp(sum(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 135 & crossthetaprec < 225)/sum(goodcells{g}));
    end
    thetamap = circshift(thetamap,size(thetamap,1)/2);
    thetamap = repmat(thetamap,[2 1]);
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:2*nphs,'w');
    if g ~= 2
        hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    end
    set(gca,'Clim',[0 0.8],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+9);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 225 & crossthetaprec < 315,:,:),1));
    thetamapref = squeeze(nanmean(thetafield{2}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 225 & crossthetaprec < 315,:,:),1));
    if g==2
        disp(sum(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 225 & crossthetaprec < 315)/sum(goodcells{g}));
    end
    thetamap = circshift(thetamap,size(thetamap,1)/2);
    thetamap = repmat(thetamap,[2 1]);
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:2*nphs,'w');
    if g ~= 2
        hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    end
    set(gca,'Clim',[0 0.8],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+12);
    thetamap = squeeze(nanmean(thetafield{g}(goodcells{g}  & goodthetacells{g},:,:),1));
    thetamapref = squeeze(nanmean(thetafield{2}(goodcells{g}  & goodthetacells{g},:,:),1));
    if g==2
        disp('% neurons with significant theta phase shift')
        disp(100*sum(goodcells{g}  & goodthetacells{g})/sum(goodcells{g}));
    end
    thetamap = circshift(thetamap,size(thetamap,1)/2);
    thetamap = repmat(thetamap,[2 1]);
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:2*nphs,'w');
    if g ~= 2
        hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    end
    set(gca,'Clim',[0 0.8],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    cmap = cmocean('dense');
    colormap(cmap)
end

figure('Name','Theta precession3')
for g = [2 1 3]
    subplot(5,3,g);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells{g}  & goodthetacells{g} & ((crossthetaprec >= 0 & crossthetaprec < 45) | (crossthetaprec <= 360 & crossthetaprec > 315)),:,:),1));
    thetamapref = squeeze(nanmean(thetafieldcorr{2}(goodcells{g}  & goodthetacells{g} & ((crossthetaprec >= 0 & crossthetaprec < 45) | (crossthetaprec <= 360 & crossthetaprec > 315)),:,:),1));
%     [x0,y0] = meshgrid(1:size(thetamap,2),1:size(thetamap,1));
%     [xinterp,yinterp] = meshgrid(1:0.5:size(thetamap,2),1:0.5:size(thetamap,1));
%     numBins = numel(1:0.5:size(thetamap,2));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:nphs,'w');
    hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    set(gca,'Clim',[0 1],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+3);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 45 & crossthetaprec < 135,:,:),1));
    thetamapref = squeeze(nanmean(thetafieldcorr{2}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 45 & crossthetaprec < 135,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:nphs,'w');
    hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    set(gca,'Clim',[0 1],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+6);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 135 & crossthetaprec <= 225,:,:),1));
    thetamapref = squeeze(nanmean(thetafieldcorr{2}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 135 & crossthetaprec <= 225,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:nphs,'w');
    hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    set(gca,'Clim',[0 1],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+9);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 225 & crossthetaprec <= 315,:,:),1));
    thetamapref = squeeze(nanmean(thetafieldcorr{2}(goodcells{g}  & goodthetacells{g} & crossthetaprec >= 225 & crossthetaprec <= 315,:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:nphs,'w');
    hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    set(gca,'Clim',[0 1],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    subplot(5,3,g+12);
    thetamap = squeeze(nanmean(thetafieldcorr{g}(goodcells{g}  & goodthetacells{g},:,:),1));
    thetamapref = squeeze(nanmean(thetafieldcorr{2}(goodcells{g}  & goodthetacells{g},:,:),1));
    imagesc(thetamap);
    hold on;plot(getCircularAverage(thetamap',0,0.1,0.05),1:nphs,'w');
    hold on;plot(getCircularAverage(thetamapref',0,0.1,0.05),1:nphs,'w--');
    set(gca,'Clim',[0 1],'Xlim',[numBins/4 3*numBins/4],'Ydir','normal');
    colormap('jet')
end

end
end

%to look at single cell precession: run this in debug mode
% f1 = figure;f2 = figure;
% for icell = 1:17000
% if goodthetacells(icell) && goodcells(icell) && cellprop.animal(icell) == 7 && cellprop.iseries(icell) == 3
% map2d = squeeze(cellprop.field2dXPhstheta{2}(icell,:,:));
% figure(f1);
% colormap(bone)
% subplot(4,1,[2 3]);
% hold off;
% imagesc(-(map2d - min(map2d(:)))./(max(map2d(:)) - min(map2d(:))));
% hold on;plot(cellprop.field2dXthetapos{2}(icell,:),1:18,'w');
% set(gca,'Clim',[-1 0],'Ydir','normal');
% subplot(4,1,1);
% hold off;
% plot(mean(map2d,1));
% subplot(4,1,4);
% hold off;
% plot(mean(map2d,2));
% figure(f2);
% colormap(parula)
% hold off;
% imagesc(squeeze(cellprop.field2dXRefcorrtheta{2}(icell,:,:)));
% hold on;plot(cellprop.field2dXcorrthetamax{2}(icell,:),1:18,'w');
% hold on;plot(cellprop.field2dXcorrthetamax{2}(icell,:)+cellprop.field2dXcorrthetamaxSE{2}(icell,:),1:18,'y');
% hold on;plot(cellprop.field2dXcorrthetamax{2}(icell,:)-cellprop.field2dXcorrthetamaxSE{2}(icell,:),1:18,'y');
% set(gca,'Clim',[-1 1],'Ydir','normal');
% pause;
% end
% end