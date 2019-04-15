function f1h = plotsortedunits(f1h,X,sidx_tograph,prev,xAxisInMsec,col,sidx,d)
% PLOTSORTEDUNITS
% Change log:
% plotunclassified.m and plotclassunits.m created by Dario Ringach
%   2009-02 AZ Merged plotunclassified2.m and plotclassunits2.m

yl = zeros(4,2);
for i = 1:4 % the sorted waveforms
   set(f1h.fig.f1,'CurrentAxes',f1h.axes.unit(i)); cla; hold on;
   
   if ~isempty(sidx_tograph{i})
      if howManyClasses(sidx_tograph(1:4)) > 1
         %AZ20090106: plot (lightly) all other traces
         plot(xAxisInMsec,X(:,vertcat(sidx_tograph{setdiff(1:4,i)})),...
            'Color', [0.8 0.8 0.8]);
      end
      f1h.plots(i).unit = line(xAxisInMsec,X(:,sidx_tograph{i}),'Color',col{i});

      yl = [yl; ylim]; %#ok<AGROW>
      My = max(max(X(:,sidx_tograph{i})));
      my = min(min(X(:,sidx_tograph{i})));
      ry = My-my;
%AZ20081210: change xaxis to msec
      axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
      
      %AZ20090127: plot mean waveform  %AZ20090317: for ALL waveforms
      f1h.plots(i).unitavg = plot(xAxisInMsec,mean(X(:,sidx{i}),2)','w-','LineWidth',1.5);

      %AZ20090212: plot previously sorted mean waveform
      for j = 1:size(prev,2)
         if exist('prev','var') && sum(size(prev)) > 1 && j <= size(prev,2) && ...
               isfield(prev(j),'unit') && isfield(prev(j).unit,'gmm') && ...
               size(prev(j).unit.gmm.icell,2) >= i && prev(j).unit.gmm.icell(i)
            switch size(prev(j).unit.prototype,2)
               case 48
                  if     length(xAxisInMsec) == 48
                     f1h.plots(i).prevunitavg = plot(xAxisInMsec,prev(j).unit.prototype      ,'k--','LineWidth',1.5);
                  elseif length(xAxisInMsec) == 28
                     f1h.plots(i).prevunitavg = plot(xAxisInMsec,prev(j).unit.prototype(7:34),'k--','LineWidth',1.5);
                  end
               case 28
                  if     length(xAxisInMsec) == 48
                     f1h.plots(i).prevunitavg = plot(xAxisInMsec(7:34),prev(j).unit.prototype,'k--','LineWidth',1.5);
                  elseif length(xAxisInMsec) == 28
                     f1h.plots(i).prevunitavg = plot(xAxisInMsec      ,prev(j).unit.prototype,'k--','LineWidth',1.5);
                  end
               otherwise
                  warning('uhoh') %#ok<WNTAG>
            end
         end
      end
      
      set(gca, 'Visible', 'on' ,'UserData',[0 i],'ButtonDownFcn',@highlightPointTrace);%, 'Box', 'off', 'Color', get(gcf, 'Color'));

      %AZ20090107: insert d' above chart
%       d(i) = dprime(V,sidx{i},vertcat(sidx{setdiff(1:4,i)}));
      set(f1h.txt.unit(i).status,'String',['d'' = ',num2str(d(i),'%0.4g')]);
   else
      set(gca, 'Visible', 'off');
      set(f1h.txt.unit(i).status,'String','');
   end
end


% % Global scaling...
% if get(f1h.gs,'Value') && nnz(yl)
%    yl = [min(yl(:,1)) max(yl(:,2))];
%    for i = 1:4
% %       axes(f1h.axes.unit(i));
%       ylim(f1h.axes.unit(i),yl);
%    end
%    clear yl; %don't want to carry over scaling from previous views
% end