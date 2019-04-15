function [f1h,gmm,k] = createPCellipse(f1h,gmm,col)
% CREATEPCELLIPSE.M

selectedAxes = zeros(1,2);

if isfield(gmm,'mu')
   iclass = size(gmm.mu,1) + 1;
else
   iclass = 1;
end

for          axisNum =  1:2
   % select plot
   [f1h,selectedAxes] = choosePCAxis(f1h,selectedAxes,axisNum);
   if isempty(nonzeros(selectedAxes))
      gmm = []; k = [];
      return;
   end
   axes(f1h.axes.pc(selectedAxes(axisNum))); hold on; drawnow;
   
   set(f1h.txt.status,'String',['Adjust ellipses (plot #',...
      num2str(selectedAxes(axisNum)),').']);
   set(f1h.txt.nchan ,'String','Double-click ellipse when done.');
   
   ctr = get(f1h.axes.pc(selectedAxes(axisNum)),'CurrentPoint');
   ctr = ctr(1,1:2);
      
   if        axisNum == 1
%       ctr = myginput(1,'crosshair'); % max 1 points, custom pointer
      k = size(ctr,1);
      
      % MORE INITIALIZATION
      Alim = zeros(2,3,k);
      starting_ellipses = {zeros(k,4); zeros(k,4);};

      gmm.mu(iclass,:)      = zeros(k,3);
      gmm.Sigma(:,:,iclass) = zeros(3,3,k);
      

      % First calculate max widths possible (hardcoded to 1/4 view)
      starting_ellipses{axisNum}(:,3)   = ...
            min( [ones(k,1)*diff(get(gca,'XLim'))/4 ...
         ctr(:,1)-ones(k,1)*min( get(gca,'XLim'))     ],[],2);
      starting_ellipses{axisNum}(:,4)   = ...
            min( [ones(k,1)*diff(get(gca,'YLim'))/4 ...
         ctr(:,2)-ones(k,1)*min( get(gca,'YLim'))     ],[],2);
      
      % Now append selected centers
      starting_ellipses{axisNum}(:,1:2) = ...
         ctr-(starting_ellipses{axisNum}(:,3:4)/2);
   elseif    axisNum == 2
      lastDimensionNum    = find(sum(sum(Alim==0,3),1)==k*2);
      currentPCAxes       = setdiff(1:3,4-selectedAxes(2));
      axisConstrainByPrev = find(currentPCAxes~=lastDimensionNum);
      % 1 = x-, 2 = y-axis
      
      % Constrain size by previous ellipse placement.  If axisConstrainByPrev =
      % 1:   previous ellipse x-axis => starting_ellipses{2}(:,3)  or
      % 2:   previous ellipse y-axis => starting_ellipses{2}(:,4)
      starting_ellipses{axisNum}(:, 2+axisConstrainByPrev) = squeeze(           ...
         2*(mean(Alim(:,currentPCAxes(axisConstrainByPrev),:)) -                ...
                 Alim(1,currentPCAxes(axisConstrainByPrev),:)));

      % Constrain size by current axes.                If axisConstrainByPrev =
      % 1:   current y-axis => starting_ellipses{2}(:,4)  or
      % 2:   current x-axis => starting_ellipses{2}(:,3)
      if     axisConstrainByPrev == 1, caLim = 'YLim';
      elseif axisConstrainByPrev == 2, caLim = 'XLim'; end
      starting_ellipses{axisNum}(:, 5-axisConstrainByPrev) =                    ...
         ones(k,1)*min( [ diff(get(gca,caLim))/4                                ...
                          mean(get(gca,caLim)) - min( get(gca,caLim)) ],[],2);

      % Center at mean of constraints, based on the 2 set directly above
      % If axisConstrainByPrev =
      % 1:   [prev(x) y-axis] - s_e{2}(:,3:4) => starting_ellipses{2}(:,1:2)  or
      % 2:   [prev(y) x-axis] - s_e{2}(:,4:3) => starting_ellipses{2}(:,2:1)
      starting_ellipses{axisNum}(:,[  axisConstrainByPrev 3-axisConstrainByPrev]) =...
         [ squeeze(mean(Alim(:,currentPCAxes(axisConstrainByPrev),:)))          ...
                                       ones(k,1)*mean(get(gca,caLim)) ] -       ...
     (starting_ellipses{axisNum}(:,[2+axisConstrainByPrev 5-axisConstrainByPrev])/2);
  
      %% TODO: merge lines below with the line directly above
      %% TODO: FIX lines below
      % AZ20090623: Overwrite position of axis unconstrained by prev
      % ellipse with clicked point
      if     ctr(:, 3-axisConstrainByPrev) < ...
             min( get(gca,caLim))          + ...
             starting_ellipses{axisNum}(:, 5-axisConstrainByPrev)/2
         % use min if clicked below it
         starting_ellipses{axisNum}(:, 3-axisConstrainByPrev) =              ...
                                         min( get(gca,caLim)) +              ...
         starting_ellipses{axisNum}(:, 5-axisConstrainByPrev)/2;
         
      elseif ctr(:, 3-axisConstrainByPrev) > ...
             max( get(gca,caLim))          + ...
             starting_ellipses{axisNum}(:, 5-axisConstrainByPrev)/2
         % use max if clicked above it
         starting_ellipses{axisNum}(:, 3-axisConstrainByPrev) =              ...
                                         max( get(gca,caLim)) +              ...
         starting_ellipses{axisNum}(:, 5-axisConstrainByPrev)/2;
         
      else % place mean where clicked
         starting_ellipses{axisNum}(:, 3-axisConstrainByPrev) =              ...
                                ctr(:, 3-axisConstrainByPrev) -              ...
         starting_ellipses{axisNum}(:, 5-axisConstrainByPrev)/2;
      end
   end
   
%    for iclass = 1:k
      if     axisNum == 1
         fcn = makeConstrainToRectFcn('imellipse',...
                           get(gca,'XLim'),get(gca,'YLim'));
      elseif axisNum == 2 % CONSTRAIN, BASED ON 1st PASS
         if     axisConstrainByPrev == 1,
            fcn = makeConstrainToRectFcn('imellipse',          ...
               Alim(:,currentPCAxes(axisConstrainByPrev),1),   ... %1 was iclass
               get(f1h.axes.pc(selectedAxes(axisNum)),'YLim')    );
         elseif axisConstrainByPrev == 2,
            %% FIXED?: Can sometimes extend PC3 in plot 2 after doing plot 3.
            fcn = makeConstrainToRectFcn('imellipse',          ...
               get(f1h.axes.pc(selectedAxes(axisNum)),'XLim'), ...
               Alim(:,currentPCAxes(axisConstrainByPrev),1)      ); %1 was iclass
         end
      end

      ellipse = imellipse(f1h.axes.pc(selectedAxes(axisNum)),   ...
                                 starting_ellipses{axisNum}(1,:)); %1 was iclass
      setPositionConstraintFcn(ellipse,fcn);
      setColor(ellipse,col{iclass});
      position = wait(ellipse);

      % Alim(:,1) = [min(ellipse axis 1); max(ellipse axis 1)]
%       if     axisNum == 1
         Alim(:,setdiff(1:3,4-selectedAxes(axisNum)),1) = ... %1 was iclass
                 [ min(position(:,1)) min(position(:,2)); ...
                   max(position(:,1)) max(position(:,2))    ];
%          [ min(position(logical(position(:,1)==min(position(:,1))),1))  ...
%            min(position(logical(position(:,2)==min(position(:,2))),2)); ...
%            max(position(logical(position(:,1)==max(position(:,1))),1))  ...
%            max(position(logical(position(:,2)==max(position(:,2))),2))     ];
%       elseif axisNum == 2 % CONSTRAIN, BASED ON 1st PASS
%          lastAxis = find(currentPCAxes==lastDimensionNum);
%          Alim(:,lastDimensionNum                  ,iclass) =            ...
%          [ min(position(logical(position(:,lastAxis) ==            ...
%                             min(position(:,lastAxis))),lastAxis)); ...
%            max(position(logical(position(:,lastAxis) ==            ...
%                             max(position(:,lastAxis))),lastAxis))     ];
%       end
      gmm.mu(iclass,setdiff(1:3,4-selectedAxes(axisNum))) = mean(position);

      % ASK FOR PC3 CENTROID & VERT STD
      delete(ellipse);
      % BUG-HACK: pointer for some reason remains a "fleur" <+>
      set(gcf,'Pointer','arrow')
      plot(position(:,1),position(:,2),'Color',col{iclass});
%    end % END for class
end    % END for selected axes

f1h = choosePCAxis(f1h,selectedAxes,3);
set(f1h.txt.nchan ,'String','');
set(f1h.txt.status,'String',''); drawnow;

%% Build GMM struct (means & covariance matrices)
% for iclass = 1:k
   % square root of the eigenvalues of cov mat correspond to ellipse axis lengths,
   % and the corresponding eigenvectors define the ellipse axis directions
   % Construct A with ellipse axes
   A = eye(3); % because ellipses aren't rotated
   % d(1) = length of ellipse axis 1
   d = (abs(Alim(2,:,:)-Alim(1,:,:))./2).^2;
   d = d ./ 3^2; % SCALE BY 1/9
   % gmm.Sigma = Covariance Matrix
   gmm.Sigma(:,:,iclass) = A*diag(d)*inv(A);
% end

k = iclass;

for selectedAxes = 1:3
   axes(f1h.axes.pc(selectedAxes)); hold on;
   f1h = plotGMMellipses(f1h,gmm,k,col,selectedAxes);
end