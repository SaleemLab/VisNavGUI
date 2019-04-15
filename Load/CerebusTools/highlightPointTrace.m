function highlightPointTrace(hObject,eventdata)
% HIGHLIGHTPOINTTRACE.M
% When you click on a trace, highlight its PCA point, and vice versa
%% TODO: un-globalize
global f1h

% SamplingRateInKHZ = 30;
% xAxisInMsec = (0:47)/SamplingRateInKHZ;
col    = {[1 0 0] [0 0 1] [0 1 0] [1 0 1] [0.5 0.5 0.5]};

% sel_idx = get(hObject,'UserData');
% unit_class = sel_idx(2); sel_idx = sel_idx(1);

type = get(get(hObject,'Parent'),'Tag');
if strcmp(type,'panel.unitcontrols')
   type = get(hObject,'Tag');
   if strcmp(type,'btngrp.bringtotop')
      type = get(eventdata.NewValue,'Tag');
   end
end
% type_paren_idx = regexp(type,'[()]');
% unit_class = str2double(type(type_paren_idx(1)+1:type_paren_idx(2)-1));
% type = type(1:type_paren_idx(1)-1);



switch type
   case 'unclassified'
      % Clicked waveform trace
      % HACK: plot traces, add in hidden points: (10,unit_class) and
      % (11,sel_idx)
      xi = get(hObject,'YData');
      sel_idx = xi(end-1);
      unit_class = xi(end);
   case {'pc(1).ax'; 'pc(2).ax'; 'pc(3).ax'}
      % Clicked PCA point
      % identify sel_idx = index of trace, unit_class
      PCAselAx = str2double(type(4));
      for i = 1:size(f1h.plots,2)
         if ~isempty(f1h.V_PCA{i}) && hObject == f1h.plots(i).pc{PCAselAx}
            unit_class = i;
            break;
         end
      end
      
      PCAselXY = get(get(hObject,'Parent'),'CurrentPoint');
      dim = 1:3;
      n = dim(dim~=(4-PCAselAx));
      
      if ~exist('unit_class','var')
         unit_class = [];
         if ~isempty(f1h.V_PCA{5}) % from plotunclassified2
            unit_class = 5;
         else%if ~isempty(sidx_tograph) % from plotclassunits2
            for i = 1:4
               if ~isempty(f1h.V_PCA{i})
                  unit_class = [unit_class i];
               end
            end
         end
      end
%       unit_class = find(hObject == [PCAh{:,PCAselAx}]);
      minval  =  ones(5,1);
      sel_idx = zeros(5,1);
      for i = unit_class
         [minval(i),sel_idx(i)] = min(sqrt( ( f1h.V_PCA{i}(:,n(1)) - PCAselXY(1,1) ).^2 + ...
                                            ( f1h.V_PCA{i}(:,n(2)) - PCAselXY(1,2) ).^2   ));
      end
      
      [minval,unit_class] = min(minval); clear minval; %#ok<ASGLU>
      sel_idx             = sel_idx(unit_class);
   case {'radio.bringtotop(1)'; 'radio.bringtotop(2)'; 'radio.bringtotop(3)'; ...
         'radio.bringtotop(4)'; 'radio.bringtotop(5)'; ''}
      % Clicked a BringToTop radio button
%       type = get(hObject,'Tag');
      if strcmp(type(1:16),'radio.bringtotop')
         unit_class = str2double(regexp(type,'[0-9]','match'));
      else
         warning('couldn''t identify which graph you clicked') %#ok<WNTAG>
      end
   otherwise
      warning('couldn''t identify which graph you clicked') %#ok<WNTAG>
end


% Bring selected class to top [TRACES]
         top_children = f1h.plots(unit_class).unclassified;
unclassified_children = get(f1h.axes.unclassified,'Children');
% Find duplicates
[duplicates,dupi]     = setdiff(top_children,unclassified_children); %#ok<ASGLU>
% Clear duplicates
top_children(dupi)    = [];
       other_children = setdiff(unclassified_children,top_children);
% set(f1h.axes.unclassified,'Children',[top_children; other_children(end:-1:1)]);

% MEAN TRACES to top
if isfield(f1h.means,'old') && unit_class <= size(f1h.means.old.line,2) % if    loaded units for this class
        mean_children = ismember(other_children,...
           [ f1h.means.new.line(unit_class)  f1h.means.new.bord(unit_class)   ...
             f1h.means.old.line{unit_class}' f1h.means.old.bord{unit_class}' ] );
else                                        % if NO loaded units for this class
       mean_children = ismember(other_children,...
           [ f1h.means.new.line(unit_class)  f1h.means.new.bord(unit_class)  ] );
end
        % must flip so borders stay below lines
        mean_children = flipud(other_children(mean_children));
       other_children = setdiff(other_children,mean_children);
set(f1h.axes.unclassified,'Children',[  mean_children; top_children; ...
                                       other_children(end:-1:1)      ]);


% Bring selected class to top [PCA POINTS]
for i = 1:3
            top_children = f1h.plots(unit_class).pc{i};
   unclassified_children = get(f1h.axes.pc(i),'Children');
   % Find duplicates
   [duplicates,dupi]     = setdiff(top_children,unclassified_children); %#ok<ASGLU>
   % Clear duplicates
   top_children(dupi) = [];
          other_children = setdiff(unclassified_children,top_children);
   set(f1h.axes.pc(i),'Children',[top_children; other_children(end:-1:1)]);
end
clear duplicates;

clearHighlights;

% Change Radio Button Selection
set(f1h.radio.bringtotop(unit_class),'Value',1);

% Don't highlight anything if no sel_idx specified: user clicked in a ch* axis
% if ~f1h.plotInteractive || ~exist('sel_idx','var'), return; end
if ~exist('sel_idx','var'), return; end

% Read data out of gui
xAxisInMsec = get(f1h.plots(unit_class).unclassified(sel_idx),'XData')';
xi          = get(f1h.plots(unit_class).unclassified(sel_idx),'YData')';
xAxisInMsec = xAxisInMsec(1:end-2);
xi          =          xi(1:end-2);

% V(1)        = get(f1h.plots(unit_class).pc{1}                 ,'XData')';
% V(2)        = get(f1h.plots(unit_class).pc{1}                 ,'YData')';
% V(3)        = get(f1h.plots(unit_class).pc{2}                 ,'YData')';
V =                      get(f1h.plots(unit_class).pc{1}      ,'XData')';
V = [V(     sel_idx   ); get(f1h.plots(unit_class).pc{1}      ,'YData')'];
V = [V([1   sel_idx+1]); get(f1h.plots(unit_class).pc{2}      ,'YData')'];
V =  V([1:2 sel_idx+2]);
% V = f1h.V_PCA{unit_class}(sel_idx,1:3);

% V           = [V{1}(sel_idx) V{2}(sel_idx) V{3}(sel_idx)];

% Plot highlights
set(gcf,'CurrentAxes',f1h.axes.unclassified); hold on;
f1h.highlights.unc(1) = plot(xAxisInMsec,xi,'w','LineWidth',4,...
   'Parent',f1h.axes.unclassified);
f1h.highlights.unc(2) = plot(xAxisInMsec,xi,'Color',col{unit_class},...
   'Parent',f1h.axes.unclassified);
   
dim = 1:3;
for j = fliplr(dim)
   n = dim(dim~=j);

   set(gcf,'CurrentAxes',f1h.axes.pc(4-j)); hold on;
   f1h.highlights.pc{4-j}(2) = plot(V(n(1)),V(n(2)),'.','Color',col{unit_class},...
      'markersize',16);
   f1h.highlights.pc{4-j}(1) = plot(V(n(1)),V(n(2)), ['w' 'o'],...
      'markersize',10,'LineWidth',1,'MarkerEdgeColor','k');
end

end %end function