function PopPlotJointCorrelation2D(popresCorr2D,popresCorr1D,popresCorr2D_blank,popresCorr1D_blank,FaveSessions)
if nargin < 3
    popresCorr2D_blank = [];
    popresCorr1D_blank = [];
end
if nargin < 5
    FaveSessions = false;
end
if FaveSessions
    CrossName = 'CA1V1Cross';
else
    CrossName = 'all_CA1V1Cross';
end
figCorr1 = figure('Name','Posterior Correlation (traj)');
titlestr{1} = 'low';
titlestr{2} = 'med';
titlestr{3} = 'high';
for g = [2 1 3]
    if ~isempty(popresCorr2D)
        subplot(3,3,g);
        imagesc(popresCorr2D.(CrossName){g});
        set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[-0.05 0.05],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
        xlabel('CA1 decoded position');
        ylabel('V1 decoded position');
        title(titlestr{g})
        hold on;plot([1 100],[1 100],'k');
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
        
        subplot(3,3,g + 3);
        imagesc(popresCorr2D.([CrossName '_shuffled']){g});
        set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[-0.05 0.05],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
        xlabel('CA1 decoded position');
        ylabel('V1 decoded position');
        title([titlestr{g} ' shf']);
        hold on;plot([1 100],[1 100],'k');
        hold on;plot([50 50],[1 100],'k');
        hold on;plot([1 100],[50 50],'k');
    end
    if ~isempty(popresCorr1D)
        subplot(3,3,g + 6);
        if FaveSessions
            ciplot(gca,[],popresCorr1D.(CrossName){g},popresCorr1D.([CrossName '_SE']){g},0.5,'k');
%             hold on;ciplot(gca,[],popresCorr1D.([CrossName '_shuffled']){g},popresCorr1D.([CrossName '_shuffled_SE']){g},0.5,[0.5 0.5 0.5]);
        else
            plot(popresCorr1D.(CrossName){g},'k');
            hold on;plot(popresCorr1D.([CrossName '_shuffled']){g},'k--');
        end
        set(gca,'Xlim',[1 100],'Ylim',[-0.1 0.1],'PlotBoxAspectRatio', [1 1 1]);
        xlabel('CA1-V1 decoded position');
        ylabel('noise correlation (Post.)');
        hold on;plot([50 50],[1 100],'k');
    end
end

if ~isempty(popresCorr2D_blank) || ~isempty(popresCorr1D_blank)
    figCorr2 = figure('Name','Posterior Correlation (blank)');
    for g = 2
        if ~isempty(popresCorr2D_blank)
            subplot(3,3,g);
            imagesc(popresCorr2D_blank.([CrossName]){g});
            set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[-0.05 0.05],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
            xlabel('CA1 decoded position');
            ylabel('V1 decoded position');
            title(titlestr{g})
            hold on;plot([1 100],[1 100],'k');
            hold on;plot([50 50],[1 100],'k');
            hold on;plot([1 100],[50 50],'k');
            
            subplot(3,3,g + 3);
            imagesc(popresCorr2D_blank.([CrossName '_shuffled']){g});
            set(gca,'Xlim',[1 100],'Ylim',[1 100],'Clim',[-0.05 0.05],'Ydir','normal','PlotBoxAspectRatio', [1 1 1]);
            xlabel('CA1 decoded position');
            ylabel('V1 decoded position');
            title([titlestr{g} ' shf']);
            hold on;plot([1 100],[1 100],'k');
            hold on;plot([50 50],[1 100],'k');
            hold on;plot([1 100],[50 50],'k');
        end
        if ~isempty(popresCorr1D_blank)
            subplot(3,3,g + 6);
            if FaveSessions
                ciplot(gca,[],popresCorr1D_blank.(CrossName){g},popresCorr1D_blank.([CrossName '_SE']){g},0.5,'k');
%                 hold on;ciplot(gca,[],popresCorr1D_blank.([CrossName '_shuffled']){g},popresCorr1D_blank.([CrossName '_shuffled_SE']){g},0.5,[0.5 0.5 0.5]);
            else
                plot(popresCorr1D_blank.(CrossName){g},'k');
                hold on;plot(popresCorr1D_blank.([CrossName '_shuffled']){g},'k--');
            end
            set(gca,'Xlim',[1 100],'Ylim',[-0.1 0.1],'PlotBoxAspectRatio', [1 1 1]);
            xlabel('CA1-V1 decoded position');
            ylabel('noise correlation (Post.)');
            hold on;plot([50 50],[1 100],'k');
        end
    end
end
end