function PopPlotSpkCorr(popresSpkCorr,popresSpkCorr_blank,F_SEcells)
if nargin < 2
    popresSpkCorr_blank = [];
end
if nargin < 3
    F_SEcells = true;
end
figcorr = figure('Name','Spike Count Noise Correlation (r_sc)');
gbase = 2;
g = gbase;

c{1} = [0 0 0];
c{2} = [0.5 0.5 0.5];
c{3} = [1 0 0];
c{4} = [0.5 0 0];
Ylim = [-0.1 0.1];
for k = 1:4
    if ismember(k,[1 2])
        popspkcorr = popresSpkCorr;
        errstr = 'traj';
    else
        popspkcorr = popresSpkCorr_blank;
        errstr = 'blank';
    end
    if ~isempty(popspkcorr)
        goodcells = ~isnan(popspkcorr.CA1V1CrossCA1V1_All{g}) & popspkcorr.CA1V1CrossPvalCA1V1_All{g} < 0.05;% & popspkcorr.CA1V1CrossCA1V1_All{g}>0;
        if ismember(k,[1 3])
            Y = popspkcorr.CA1V1CrossCA1V1_All{g};
            shfstr = '{measured}';
        else
            Y = popspkcorr.CA1V1CrossCA1V1_shuffled_All{g};
            shfstr = '{shuffled}';
        end
        
        subplot(4,5,(k-1)*5 + 1);
%         scatter(popspkcorr.CA1V1fieldcorrCA1V1_All{g}(goodcells),Y(goodcells),'Marker','.','MarkerEdgeColor',[0.1 0.1 0.1],'MarkerEdgeAlpha',0.01);
        nbinscorr = 20;
        X = normalise1var(popspkcorr.CA1V1fieldcorrCA1V1_All{g}, nbinscorr, [], [-1 1]);
        if F_SEcells
            [corrmap,x,~,corrmapSE] = fast1Dmap(X(goodcells),Y(goodcells),1,1,nbinscorr/1,false);
        else
            if ismember(k,[1 3])
                corrmap = popspkcorr.CA1V1CrossAveCA1V1_Corr{g};
                corrmapSE = popspkcorr.CA1V1CrossAveCA1V1_Corr_SE{g};
                x = linspace(-1,1,numel(corrmap));
            else
                corrmap = popspkcorr.CA1V1CrossAveCA1V1_Corr_shuffled{g};
                corrmapSE = popspkcorr.CA1V1CrossAveCA1V1_Corr_shuffled_SE{g};
                x = linspace(-1,1,numel(corrmap));
            end
        end
        hold on;
        ciplot(gca,linspace(-1,1,numel(x)),corrmap,corrmapSE,0.5,'k');
        set(gca,'Xlim',[-1 1],'Ylim',Ylim);
        xlabel('r_{field}');
        ylabel(['r_{sc} ' errstr '_' shfstr]);
        
        subplot(4,5,(k-1)*5 + 2);
%         scatter(popspkcorr.CA1V1fielddistCA1V1_All{g}(goodcells)-50,Y(goodcells),'Marker','.','MarkerEdgeColor',[0.1 0.1 0.1],'MarkerEdgeAlpha',0.01);
        nbinscorr = 20;
        X = normalise1var(popspkcorr.CA1V1fielddistCA1V1_All{g}, nbinscorr, [], [0 100]);
        if F_SEcells
            [corrmap,x,~,corrmapSE] = fast1Dmap(X(goodcells),Y(goodcells),1,1,nbinscorr/1,true);
        else
            if ismember(k,[1 3])
                corrmap = popspkcorr.CA1V1CrossAveCA1V1_dPeak{g};
                corrmapSE = popspkcorr.CA1V1CrossAveCA1V1_dPeak_SE{g};
                x = linspace(-1,1,numel(corrmap));
            else
                corrmap = popspkcorr.CA1V1CrossAveCA1V1_dPeak_shuffled{g};
                corrmapSE = popspkcorr.CA1V1CrossAveCA1V1_dPeak_shuffled_SE{g};
                x = linspace(-1,1,numel(corrmap));
            end
        end
        hold on;
        ciplot(gca,linspace(-50,50,numel(x)),corrmap,corrmapSE,0.5,'k');
        set(gca,'Xlim',[-50 50],'Ylim',Ylim);
        xlabel('Dist_{field}');
        ylabel(['r_{sc} ' errstr '_' shfstr]);
        
        subplot(4,5,(k-1)*5 + 3);
%         scatter(popspkcorr.CA1V1fieldXcorrmaxPosCA1V1_All{g}(goodcells)-50,Y(goodcells),'Marker','.','MarkerEdgeColor',[0.1 0.1 0.1],'MarkerEdgeAlpha',0.01);
        nbinscorr = 20;
        X = normalise1var(popspkcorr.CA1V1fieldXcorrmaxPosCA1V1_All{g}, nbinscorr, [], [0 100]);
        if F_SEcells
            [corrmap,x,~,corrmapSE] = fast1Dmap(X(goodcells),Y(goodcells),1,1,nbinscorr/1,true);
        else
            if ismember(k,[1 3])
                corrmap = popspkcorr.CA1V1CrossAveCA1V1_dCorr{g};
                corrmapSE = popspkcorr.CA1V1CrossAveCA1V1_dCorr_SE{g};
                x = linspace(-1,1,numel(corrmap));
            else
                corrmap = popspkcorr.CA1V1CrossAveCA1V1_dCorr_shuffled{g};
                corrmapSE = popspkcorr.CA1V1CrossAveCA1V1_dCorr_shuffled_SE{g};
                x = linspace(-1,1,numel(corrmap));
            end
        end
        hold on;
        ciplot(gca,linspace(-50,50,numel(x)),corrmap,corrmapSE,0.5,'k');
        set(gca,'Xlim',[-50 50],'Ylim',Ylim);
        xlabel('Dist_{field(corr)}');
        ylabel(['r_{sc} ' errstr '_' shfstr]);
        
        try
        subplot(4,5,(k-1)*5 + 4);
%         scatter(popspkcorr.CA1V1fieldXcorrmaxPosCA1V1_All{g}(goodcells)-50,popspkcorr.CA1V1fieldXcorrmaxCA1V1_All{g}(goodcells),'Marker','.','MarkerEdgeColor',[0.1 0.1 0.1],'MarkerEdgeAlpha',0.01);
        nbinscorr = 20;
        X = normalise1var(popspkcorr.CA1V1fieldXcorrmaxPosCA1V1_All{g}, nbinscorr, [], [0 100]);
        if F_SEcells
            [corrmap,x,~,corrmapSE] = fast1Dmap(X(goodcells),popspkcorr.CA1V1fieldXcorrmaxCA1V1_All{g}(goodcells),1,1,nbinscorr/1,true);
        else
            if ismember(k,[1 3])
                corrmap = popspkcorr.CA1V1fieldXcorrAveCA1V1_dCorr{g};
                corrmapSE = popspkcorr.CA1V1fieldXcorrAveCA1V1_dCorr_SE{g};
                x = linspace(-1,1,numel(corrmap));
            else
                corrmap = popspkcorr.CA1V1fieldXcorrAveCA1V1_dCorr_shuffled{g};
                corrmapSE = popspkcorr.CA1V1fieldXcorrAveCA1V1_dCorr_shuffled_SE{g};
                x = linspace(-1,1,numel(corrmap));
            end
        end
        hold on;
        ciplot(gca,linspace(-50,50,numel(x)),corrmap,corrmapSE,0.5,'k');
        set(gca,'Xlim',[-50 50],'Ylim',[0.675 0.725]);
        xlabel('Dist_{field(corr)}');
        ylabel(['r_{field} ']);
        catch
        end
        subplot(2,5,floor((k-1)/2)*5 + 5);
        hold on;
        [f,x,flo,fup] = ecdf(Y(goodcells));
        ciplot(gca,x,f,f-flo,0.5,c{k});
        set(gca,'PlotBoxAspectRatio', [1 1 1]);
        set(gca,'Xlim',[-0.2 0.2],'Ylim',[0 1]);
        xlabel('r_{sc}');
        ylabel('cdf');
    end
end
end