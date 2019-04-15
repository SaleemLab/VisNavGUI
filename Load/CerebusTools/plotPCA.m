function f1h = plotPCA(f1h,V,sidx_tograph,gmm,col,sidx,sidx_tograph_only)
% PLOTPCA
% Plots (un)sorted units in cluster plots
%
% Change log:
% plotunclassified.m and plotclassunits.m created by Dario Ringach
%   2008-08 LB added plotting of unclassified units in cluster plots
%   2008-12 AZ color-coded cluster axes. made prettier.
%   2009-02 AZ Merged plotunclassified2.m and plotclassunits2.m

clearHighlights;

f1h.V_PCA = [];

sidx_nonempty = [];
for i = 1:5
   if ~isempty(sidx_tograph{i})
      sidx_nonempty = [sidx_nonempty i]; %#ok<AGROW>
   end
end

% PLOT
if isfield(gmm,'mu')
   k = size(gmm.mu,1);
else
   k = get(f1h.pop.nclusters,'Value');
end
dim = 1:3;
for j = dim
   n = dim(dim~=j);
   set(f1h.fig.f1,'CurrentAxes',f1h.axes.pc(4-j)); cla; hold on; % cluster views (bottom axis: p2 vs p3)
   for i = sidx_nonempty
       % AZ20090710: PLOT sub-threshold separately
       line(V(sidx_tograph_only,n(1)),V(sidx_tograph_only,n(2)), ...
           'Marker', '.', 'markersize', 3, 'LineStyle','none', 'Color', 'k', ...
           'Parent',f1h.axes.pc(4-j));
   %    f1h.plots(i).pc{3} = plot(V(sidx_tograph{i},2),V(sidx_tograph{i},3), '.', ...
   %          'markersize', 3, 'Color', col{i}, ...
   %          'buttondownfcn',@highlightPointTrace,'Parent',f1h.axes.pc(3));
      f1h.V_PCA{i} = [V(sidx_tograph{i},:) (1:length(sidx_tograph{i}))'];
%       if f1h.plotInteractive
%          for a = 1:length(sidx_tograph{i})
%             f1h.plots(i).pc{4-j}(a) = line(V(sidx_tograph{i}(a),n(1)),V(sidx_tograph{i}(a),n(2)), ...
%                'Marker', '.', 'markersize', 3, 'Color', col{i}, ...
%                'buttondownfcn',@highlightPointTrace,'Parent',f1h.axes.pc(4-j));
%          end
%       else
         f1h.plots(i).pc{4-j} = line(V(sidx_tograph{i},n(1)),V(sidx_tograph{i},n(2)), ...
            'Marker', '.', 'markersize', 3, 'LineStyle','none', 'Color', col{i}, ...
            'ButtonDownFcn',@highlightPointTrace,'Parent',f1h.axes.pc(4-j));
%       end
   %    set(f1h.plots(i).pc{1}, 'ZData', V(sidx_tograph{i},1));
   end % end for nonempty
   
   % PLOT GMM ellipses: they exist only if there are classified points
   if sidx_nonempty(1) <= k && sidx_nonempty(1) > 0
      f1h = plotGMMellipses(f1h,gmm,k,col,4-j);
   end

   % if ~isempty(sidx_tograph{5})
   xr = range(V(:,n(1))); yr = range(V(:,n(2)));
   axis([min(V(:,n(1)))-0.1*xr max(V(:,n(1)))+0.15*xr min(V(:,n(2)))-0.15*yr max(V(:,n(2)))+0.1*yr])
   % else
   line(0,0,'Marker','.','Color','k','markersize',16);
   line(0,0,'Marker','.','Color','r','markersize',6);
   % end
end % end for each pca plot

% AZ20090318: display proportion of spikes shown
title(sprintf('%d/%d (%0.1f%%) spikes shown',...
          [ size(vertcat(sidx_tograph{:},sidx_tograph_only),1),size(V,1), ...
        100*size(vertcat(sidx_tograph{:},sidx_tograph_only),1)/size(V,1) ]));
       

drawnow;

end % end function plotPCA