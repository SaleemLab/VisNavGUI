function [f1h,selectedAxes] = choosePCAxis(f1h,selectedAxes,axisNum)
% CHOOSEPCAXIS.M

figSize = get(f1h.fig.f1,'Position'); figSize = figSize(3:4);
currentPoint = zeros(1,2);
ctr = 0;

% AZ20090526: Click PC Axis to choose first one.
if     axisNum == 1
   set(f1h.txt.status,'String','Click to place 1st ellipse.');
elseif axisNum == 2
   set(f1h.txt.status,'String','Click to place 2nd ellipse.');
elseif axisNum == 3
   for i = 1:3
       set(f1h.axes.pc(                    i),'Box','off','LineWidth',0.5,...
          'FontWeight','normal');
   end
   return;
else
   warning('sortnev2:createPCellipse:choosePCAxis:TooManyAxes',...
           'Pick up to 2 axes only')
end
pause(0.3); % HACK: keeps pointer as circle 2nd time
set(f1h.fig.f1,'Pointer','circle'); drawnow;

% Get PC Plot Boundaries
pcbound1x = get(f1h.axes.pc(1),'Position');
pcbound1y = pcbound1x(2) + pcbound1x(4);

pcbound2y = pcbound1x(2);
pcbound1x = pcbound1x(1);

pcbound3y = get(f1h.axes.pc(2),'Position');
pcbound3y = pcbound3y(2);

% Keep waiting until user clicks a PCA plot or presses ESC
while (ctr == 0 && ~(currentPoint(1) > pcbound1x && currentPoint(2) < pcbound1y) ) || ...
      (ctr == 1 &&   double(get(f1h.fig.f1,'CurrentCharacter')) ~= 27    )
   ctr = waitforbuttonpress;

   currentPoint = get(f1h.fig.f1,'CurrentPoint')./figSize;
end
% Quit if user pressed ESC
if ctr == 1 && double(get(f1h.fig.f1,'CurrentCharacter')) == 27
   set(f1h.txt.status,'String','Aborted. Try again.');
   set(f1h.fig.f1,'Pointer','arrow'); drawnow;
   return;
% Assign PCA plot # (values dependent on those defined in sortnev2)
elseif ctr == 0 && (currentPoint(1) > pcbound1x || currentPoint(2) < pcbound1y)
   if     currentPoint(2      ) < pcbound3y
          selectedAxes(axisNum) = 3;
   elseif currentPoint(2      ) < pcbound2y
          selectedAxes(axisNum) = 2;
   else
          selectedAxes(axisNum) = 1;
   end
   set(f1h.fig.f1,'Pointer','arrow'); drawnow;
end

% Repeat execution if selected the same plot twice
if selectedAxes(1) == selectedAxes(2)      
   set(f1h.txt.status,'String',['You already selected plot #',...
   num2str(selectedAxes(axisNum)),'. Try again.']); drawnow;
   pause(1);

   [f1h,selectedAxes] = choosePCAxis(f1h,selectedAxes,axisNum);
else
       set(f1h.axes.pc(selectedAxes(axisNum)),'Box','on' ,'LineWidth',1.5,...
          'FontWeight','bold'  );
   for i = setdiff(1:3,selectedAxes(axisNum))
       set(f1h.axes.pc(                    i),'Box','off','LineWidth',0.5,...
          'FontWeight','normal');
   end
end

end