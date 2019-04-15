function [f1h,gmm,d,sidx,sidx_tograph,C,S,U,V,Xrealigned] = clickEllipse( ...
          f1h,gmm,d,sidx,sidx_tograph,  S,U,V,Xrealigned,                 ...
          iXraw,idx_tograph,xAxisInMsec,SamplingRateInKHZ,col,threshold,prev,action)
% CLICKELLIPSE.M: resize ellipse on demand
currentObject = gco;

infos      = get(currentObject,'Tag');
unit_class = str2double(infos(7));
plotNum    = str2double(regexp(get(gca,'Tag'),'[0-9]','match'));
% If selected small (1SD) ellipse, pretend we selected the large (3SD) one
if str2double(infos(end-1))
   currentObject = f1h.plots(unit_class).ellipse(plotNum,2);
end

if strcmp(action,'edit')
   mx = min(get(currentObject,'XData'));
   Mx = max(get(currentObject,'XData'));
   my = min(get(currentObject,'YData'));
   My = max(get(currentObject,'YData'));

   % delete 1SD ellipse, make 3SD ellipse dashed
   delete(f1h.plots(unit_class).ellipse(plotNum,1));
   set(f1h.plots(unit_class).ellipse(plotNum,2),'LineStyle','--')

   fcn = makeConstrainToRectFcn('imellipse',...
                              get(gca,'XLim'),get(gca,'YLim'));
   ellipse = imellipse(gca,[mx my Mx-mx My-my]);
   setPositionConstraintFcn(ellipse,fcn);
   setColor(ellipse,col{unit_class});
   position = wait(ellipse);


   Alim(:,setdiff(1:3,4-plotNum),unit_class) = ...
      [ min(position(:,1)) min(position(:,2)); ...
        max(position(:,1)) max(position(:,2))    ];
   %% AZ20090529
   if     ceil(1-(plotNum)/3)+1 == 1, caLim = 'XData';
   elseif ceil(1-(plotNum)/3)+1 == 2, caLim = 'YData'; end
   unchangedLimit = get(f1h.plots(unit_class).ellipse(round(1-(plotNum)/3)+1,2),caLim);
   Alim(:,4-plotNum,unit_class) = [min(unchangedLimit); max(unchangedLimit)];

   gmm.mu(unit_class,setdiff(1:3,4-plotNum)) = mean(position);

   delete(ellipse);
   % % BUG-HACK: pointer for some reason remains a "fleur" <+>
   % set(gcf,'Pointer','arrow')
   plot(position(:,1),position(:,2),'Color',col{unit_class});

   A = eye(3); % because ellipses aren't rotated
      % axisLengths(1) = length of ellipse axis 1
      axisLengths = (abs(Alim(2,:,unit_class)-Alim(1,:,unit_class))./2).^2;
   axisLengths = axisLengths ./ 3^2; % SCALE BY 1/9
   % gmm.Sigma = Covariance Matrix

   %% AZ20090529
   gmm.mu(unit_class,4-plotNum ) = mean(unchangedLimit);
   gmm.Sigma(:,:,unit_class    ) = A*diag(axisLengths)/A;
   k = size(gmm.mu,1);
   
   for selectedAxes = 1:3
      set(f1h.fig.f1,'CurrentAxes',f1h.axes.pc(selectedAxes)); hold on;
      f1h = plotGMMellipses(f1h,gmm,k,col,selectedAxes);
   end
   
elseif strcmp(action,'remove')
   k = size(gmm.mu,1);
   gmm.mu    = gmm.mu(       setdiff(1:k,unit_class),:);
   gmm.Sigma = gmm.Sigma(:,:,setdiff(1:k,unit_class)  );
   k = k - 1;
   set(f1h.pop.nclusters,'Value',max([k 1]));
end

if ~isempty(gmm.mu)
   [f1h,gmm,d,sidx,sidx_tograph,C,S,U,V,Xrealigned] = cluster_gmm3sd( ...
    f1h,gmm,d,sidx,sidx_tograph,S,U,V,Xrealigned,iXraw,idx_tograph,   ...
      xAxisInMsec,SamplingRateInKHZ,col,threshold,prev);
else
   C = [];
   clearEllipses;
   gmm = struct; k = []; sidx = cell(1,5); sidx_tograph = cell(1,5);
   
   f1h = sortnev2_controls_update(false,f1h);
   [Xrealigned,idx_tograph,idx_tograph_only,sidx_tograph,sidx,f1h] = sortnev2_gui_update(...
    Xrealigned,V,iXraw,sidx_tograph,gmm,prev,xAxisInMsec,col,sidx,d,threshold,f1h,...
    'Ellipses Cleared.');
end

end