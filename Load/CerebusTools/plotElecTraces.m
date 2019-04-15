function plotElecTraces(fignum,tagname,A,X,iXraw,Xraw,xAxisInMsec,plotColors,k,...
   class_idx,threshold,yAxisMin,yAxisMax,elecs)
% Plot electrode traces function
    if nargin < 14
       elecs = A';
       elecs = elecs(~isnan(elecs))';
    end
    
   for elec = elecs
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
      if ~isnan(elec) %&& ~isempty(X{elec})
         h = findobj('Tag',['elec',num2str(elec),tagname{fignum}]);
         set(0,'CurrentFigure',get(h,'Parent'))
         set(gcf,'CurrentAxes',h(1));

         hmenu = findobj('Tag','expnev1.menu.subplotfcns(1)');

         if ~isempty(X{elec})
            % Plot Subthreshold traces, or area if more than 100
            if nnz(~iXraw{elec}) > 100
                  fill([xAxisInMsec xAxisInMsec(end:-1:1)],...
                     [mean(X{elec}(       :,~iXraw{elec})   ,2) +    ...
                       std(X{elec}(       :,~iXraw{elec}),[],2)   ;  ...
                      mean(X{elec}(end:-1:1,~iXraw{elec})   ,2) -    ...
                       std(X{elec}(end:-1:1,~iXraw{elec}),[],2)    ],...
                     plotColors{5}+[0.8 0.8 0.8],        'UIContextMenu',hmenu);
                  plot(xAxisInMsec,mean(X{elec}(:,~iXraw{elec}),2),...
                     'Color',plotColors{5},'LineWidth',2,'UIContextMenu',hmenu);
            elseif nnz(~iXraw{elec}) > 0
                  plot(xAxisInMsec,     X{elec}(:,~iXraw{elec})   ,...
                     'Color',plotColors{5},              'UIContextMenu',hmenu);
            end

            % Plot all supra-threshold traces
            for a = 1:k(elec)
               if ~isempty(intersect(class_idx{elec,a}  ,find(iXraw{elec}) ))
                  plot(xAxisInMsec,X{elec}(:,intersect(class_idx{elec,a}  ,...
                     find(iXraw{elec}) )),'Color',plotColors{a},'UIContextMenu',hmenu);
               end
            end
         else % All traces are sub-threshold--use Xraw instead of X
            % Plot Subthreshold traces, or area if more than 100
            if size(Xraw{elec},2) > 100
               fill([xAxisInMsec xAxisInMsec(end:-1:1)]',...
                  [mean(Xraw{elec}            ,2)+std(Xraw{elec}            ,[],2); ...
                   mean(Xraw{elec}(end:-1:1,:),2)-std(Xraw{elec}(end:-1:1,:),[],2)],...
                  plotColors{5}+[0.8 0.8 0.8],        'UIContextMenu',hmenu);
%                      alpha(0.5);
               plot(xAxisInMsec,mean(Xraw{elec},2),...
                  'Color',plotColors{5},'LineWidth',2,'UIContextMenu',hmenu);
            elseif ~isempty(Xraw{elec})
               plot(xAxisInMsec,     Xraw{elec}' ,...
                  'Color',plotColors{5},'LineWidth',2,'UIContextMenu',hmenu);
            else
               text(0.5,0.5,'No Data','Units','Normalized',...
                  'HorizontalAlignment','center','UIContextMenu',hmenu)
            end
         end

         %AZ20090721
         if ~isempty(Xraw{elec})
            if threshold(elec) % Plot Threshold line, if defined
               plot([xAxisInMsec(1) xAxisInMsec(end)], abs(threshold(elec))*ones(1,2),...
                  '-k','UIContextMenu',hmenu,'Tag',sprintf('elec%dplotThresh',elec));
               plot([xAxisInMsec(1) xAxisInMsec(end)],-abs(threshold(elec))*ones(1,2),...
                  '-k','UIContextMenu',hmenu);
            end

            axis([xAxisInMsec(1) xAxisInMsec(end) ...
                  [min(yAxisMin(elec),-abs(1.5*threshold(elec))) ...
                   max(yAxisMax(elec), abs(1.5*threshold(elec)))]   ]);
         end
      end
%         end
    end
end
% END Plot electrode traces function