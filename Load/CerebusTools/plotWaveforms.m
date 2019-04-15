function f1h = plotWaveforms(f1h,X,sidx_tograph,prev,xAxisInMsec,col,sidx,...
   sidx_tograph_only,threshold)
% PLOTWAVEFORMS
% Plots unsorted units in cluster plots
%
% Change log:
% plotunclassified.m and plotclassunits.m created by Dario Ringach
%   2008-08 LB corrected some matlab complaints
%   2008-12 AZ modified unclassified xaxis to display in msec, assuming a
%       certain sampling rate. made prettier
%   2009-01 AZ classunit plots now plot all other (non-selected) waveforms
%       in light gray, behind selected waveforms.  d' calculated for each
%       nonempty unit, value is displayed above the given classunit plot
%   2008-12 AZ made unclassified plot show all classified waveforms if no 
%       unclassified ones are left. modified unclassified xaxis to display
%       in msec, assuming a certain sampling rate. made prettier.
%   2009-02 AZ Merged plotunclassified2.m and plotclassunits2.m

set(f1h.fig.f1,'CurrentAxes',f1h.axes.unclassified); cla; hold on;

sidx_nonempty = [];
for i = 1:5
   if ~isempty(sidx_tograph{i})
      sidx_nonempty = [sidx_nonempty i]; %#ok<AGROW>
   end
end

%% PLOT %%
My = max(max(X(:,vertcat(sidx_tograph{sidx_nonempty}))));
my = min(min(X(:,vertcat(sidx_tograph{sidx_nonempty}))));
ry = My-my;

% AZ20090710: PLOT sub-threshold separately
if ~isempty(sidx_tograph_only)
    line(xAxisInMsec,X(:,sidx_tograph_only),'color','k','Parent',f1h.axes.unclassified);
end

for i = fliplr(sidx_nonempty) % the sorted waveforms
   % HACK: plot traces, add in hidden points: (10,i) and (11,a),
   %         where i = class, and a = index of trace
   f1h.plots(i).unclassified = line([xAxisInMsec 10 11],...
      [X(:,sidx_tograph{i}); 1:size(X(:,sidx_tograph{i}),2); i*ones(1,size(sidx_tograph{i},1))],...
        'color',col{i},'buttondownfcn',@highlightPointTrace,'Parent',f1h.axes.unclassified);
end

% AZ20081210: change xaxis to msecs
axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
xlabel('msec');

% % Show number of traces displayed
% if isfield(    f1h,'unclassified_annotation') && ...
%       ~isempty(f1h.unclassified_annotation)   && ...
%       ishandle(f1h.unclassified_annotation)
%     delete(    f1h.unclassified_annotation);
% end
% f1h.unclassified_annotation = ...
%    annotation('textbox',[0.365 0.76 0.1 0.025],'Units','normalized',...
%    'LineStyle','none','String',['N = ',num2str(size(vertcat(sidx_tograph{sidx_nonempty}),1))]);

if ~isempty(sidx_tograph{5})
   title(sprintf('%1g/%1g Unclassified (%2.1f%%)',...
      [size(sidx{5},1) size(X,2) 100*size(sidx{5},1)/size(X,2)]))
else %AZ20081211: Plot classified trains together in 'unclassified' box
   title('All Classified')
   %AZ20090107: disable 'rest' buttons
%    for j = 1:4
%       set(f1h.pshbtn.unit(j).unc,'Enable','Off');
%    end
end

%% Mean Waveforms
for i = fliplr(sidx_nonempty)
   %AZ20090127: plot mean waveform  %AZ20090317: for ALL waveforms
   f1h.means.new.bord(i) = line(xAxisInMsec,mean(X(:,sidx_tograph{i}),2)',...
      'Color','k'                                ,'LineWidth',7  ,'Visible','off');
   f1h.means.new.line(i) = line(xAxisInMsec,mean(X(:,sidx_tograph{i}),2)',...
      'Color',col{i}+0.6*double(~logical(col{i})),'LineWidth',4.5,'Visible','off');
end

%% AZ20090212/0603: plot previously sorted mean waveforms
if size(prev,2) > 0
   UnitNum = zeros(1,size(prev,2));
   for j = 1:size(prev,2)
      UnitNum(j) = find(prev(j).unit.gmm.icell==prev(j).unit.icell);
   end
else
   UnitNum = [];
end

% Initialize
for i = unique(UnitNum)
   f1h.means.old.line{i} = zeros(sum(UnitNum==i),1);
   f1h.means.old.bord{i} = zeros(sum(UnitNum==i),1);
end

for j = 1:size(prev,2)
   if       exist('prev','var')                    && ...
          sum(size(prev)) > 1                      && ...
         j <= size(prev,2)                         && ...
           isfield(prev(j),'unit')                 && ...
           isfield(prev(j).unit,'gmm')             && ...
              size(prev(j).unit.gmm.icell,2) >= i  && ...
                   prev(j).unit.gmm.icell(i)
      
      switch num2str([length(prev(j).unit.prototype) length(xAxisInMsec)])
         case {'48  48'; '28  28'}
            f1h.means.old.line{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))) = ...
                          line(xAxisInMsec      ,prev(j).unit.prototype      );
%                line( [xAxisInMsec       10 11] , [prev(j).unit.prototype       0 0] );
         case {'48  28'}
            f1h.means.old.line{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))) = ...
                          line(xAxisInMsec      ,prev(j).unit.prototype(7:34));
%                line( [xAxisInMsec       10 11] , [prev(j).unit.prototype(7:34) 0 0] );
         case {'28  48'}
            f1h.means.old.line{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))) = ...
                          line(xAxisInMsec(7:34),prev(j).unit.prototype      );
%                line( [xAxisInMsec(7:34) 10 11] , [prev(j).unit.prototype       0 0] );
         otherwise
            warning('uhoh')
      end
      set(   f1h.means.old.line{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))),...
         'LineStyle','-' ,'LineWidth',5  ,'Color','w'                ,'Visible','off');
      f1h.means.old.bord{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))) = line(...
         get(f1h.means.old.line{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))),'XData'),...
         get(f1h.means.old.line{UnitNum(j)}(sum(UnitNum(1:j)==UnitNum(j))),'YData'),...
         'LineStyle','--','LineWidth',2.5,'Color',col{UnitNum(j)}*0.75,'Visible','off');
   end
end

%%
% Comment this out if using noise estimation for Dynamic Multiphasic Filter
if threshold % Plot Threshold line, if defined
   plot([xAxisInMsec(1) xAxisInMsec(end)], abs(threshold(1))*ones(1,2),'-k');
   plot([xAxisInMsec(1) xAxisInMsec(end)],-abs(threshold(1))*ones(1,2),'-k');
   
% COMMENTED OUT: not a terribly informative implementation
%    % Plot dynamic multiphasic filter secondary thresholds
% 	  [xr,yr] = min(X,[],1);
%    a = [max([floor(mean(yr)-std(yr))-8 1]) min([ceil(mean(yr)+std(yr))+8 48])];
%    b = [min(xr + 2*threshold) max(xr + 2*threshold)];
%    plot(xAxisInMsec(a(1))*ones(1,2),ylim,'-y')
%    plot(xAxisInMsec(a(2))*ones(1,2),ylim,'-y')
%    plot(xAxisInMsec(a),b(1)*ones(1,2),'-y')
%    plot(xAxisInMsec(a),b(2)*ones(1,2),'-y')
%    
%    [xr,yr] = max(X,[],1);
%    a = [max([floor(mean(yr)-std(yr))-8 1]) min([ceil(mean(yr)+std(yr))+8 48])];
%    b = [min(xr - 2*threshold) max(xr - 2*threshold)];
%    plot(xAxisInMsec(a(1))*ones(1,2),ylim,'-m')
%    plot(xAxisInMsec(a(2))*ones(1,2),ylim,'-m')
%    plot(xAxisInMsec(a),b(1)*ones(1,2),'-m')
%    plot(xAxisInMsec(a),b(2)*ones(1,2),'-m')
end

end % end function plotWaveforms