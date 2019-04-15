function f1h = plotGMMellipses(f1h,gmm,k,col,selectedAxes,linewidth)
% PLOTGMMELLIPSES.M: Plot ellipses from GMM struct

if ~exist('linewidth','var') || isempty(linewidth)
   linewidth = 1;
end

% PC's the x-(1) and y-(2) axes represent for the selectedAxes
n = setdiff(1:3,4-selectedAxes);

%generate a vector of angles from 0 to 2*pi
theta = 0:.02:2*pi;

for unit_class = 1:k
%         plot(gmm.obj.mu(unit_class,n(1)),gmm.obj.mu(unit_class,n(2)),...
%             col{unit_class},'Marker','x','MarkerSize',16,'LineWidth',5);

   %diagonalize the covariance matrix
%         [u,lam] = eig(inv(gmm.Sigma(:,:,unit_class)));
   [u,lam] = eig(gmm.Sigma(n,n,unit_class));

   %calculate the x component of the ellipsoid for all angles
   r(:,1) = gmm.mu(unit_class,n(1)) + (  sqrt(lam(1,1)))*u(1,1)*cos(theta) + ...
                                      (  sqrt(lam(2,2)))*u(1,2)*sin(theta);
   r(:,3) = gmm.mu(unit_class,n(1)) + (3*sqrt(lam(1,1)))*u(1,1)*cos(theta) + ...
                                      (3*sqrt(lam(2,2)))*u(1,2)*sin(theta);
   %calculate the y component of the ellipsoid for all angles
   r(:,2) = gmm.mu(unit_class,n(2)) + (  sqrt(lam(1,1)))*u(2,1)*cos(theta) + ...
                                      (  sqrt(lam(2,2)))*u(2,2)*sin(theta);
   r(:,4) = gmm.mu(unit_class,n(2)) + (3*sqrt(lam(1,1)))*u(2,1)*cos(theta) + ...
                                      (3*sqrt(lam(2,2)))*u(2,2)*sin(theta);
   %plot(x,y)
   if exist('f1h','var') && ~isempty(f1h)
      f1h.plots(unit_class).ellipse(selectedAxes,1) = line(r(:,1),r(:,2),...
         'Color',col{unit_class},'UIContextMenu',f1h.menu.pcaellipse(1), ...
         'Tag',['plots(',num2str(unit_class),').ellipse(',num2str(selectedAxes),',1)'],...
         'LineWidth',linewidth);
                                %'ButtonDownFcn',@clickEllipse);
      f1h.plots(unit_class).ellipse(selectedAxes,2) = line(r(:,3),r(:,4),...
         'Color',col{unit_class},'UIContextMenu',f1h.menu.pcaellipse(1), ...
         'Tag',['plots(',num2str(unit_class),').ellipse(',num2str(selectedAxes),',2)'],...
         'LineWidth',linewidth);
                                %'ButtonDownFcn',@clickEllipse);
   else
      line(r(:,1),r(:,2),'Color',col{unit_class}, ...
         'Tag',['plots(',num2str(unit_class),').ellipse(',num2str(selectedAxes),',1)'],...
         'LineWidth',linewidth);
                                %'ButtonDownFcn',@clickEllipse);
      line(r(:,3),r(:,4),'Color',col{unit_class}, ...
         'Tag',['plots(',num2str(unit_class),').ellipse(',num2str(selectedAxes),',2)'],...
         'LineWidth',linewidth);
   end
end