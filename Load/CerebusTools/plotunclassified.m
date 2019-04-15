function plots = plotunclassified(plots)
% PLOTUNCLASSIFIED
% Plots unsorted units in cluster plots and as waveforms
%
% Change log:
% Created by Dario Ringach
% 2008-08 LB added plotting of unclassified units in cluster plots
% 2008-12 AZ color-coded cluster axes. made unclassified plot show all
%   classified waveforms if no unclassified ones are left. modified
%   unclassified xaxis to display in msec, assuming a certain sampling rate.
%   made prettier.

global X V sidx uidx lowidx highidx;

% cluster views (top axis: p1 vs p2)
h = guidata(sortnev);

% unclassified
uidx = 1:size(V,1);
%AZ20081211: %Make uidx have contain all spike train indices, then take out
             %sorted ones
uidx = setdiff(uidx,union(union(sidx{1},sidx{2}),union(sidx{3},sidx{4})));

if ~isempty(uidx)
   axes(h.p1p2);
   plot(V(uidx,1),V(uidx,2), '.', 'markersize', 3, 'Color', [0.5 0.5 0.5]);
   set(gca,'XColor',[0.4 0 0],'YColor',[0 0.4 0])

   axes(h.p1p3);
   plot(V(uidx,1),V(uidx,3), '.', 'markersize', 3, 'Color', [0.5 0.5 0.5]);
   set(gca,'XColor',[0.4 0 0],'YColor',[0 0 0.4])

   axes(h.p2p3);
   plot(V(uidx,2),V(uidx,3), '.', 'markersize', 3, 'Color', [0.5 0.5 0.5]);
   set(gca,'XColor',[0 0.4 0],'YColor',[0 0 0.4])
end

axes(h.unclassified);
cla;
uidx = 1:size(X,2);
uidx = setdiff(uidx,union(union(sidx{1},sidx{2}),union(sidx{3},sidx{4})));

if(~isempty(uidx))
    title('Unclassified')
%AZ20081210: scale xaxis to msec (assuming sampling at 30kHz)
    SamplingRateInKHZ = 30;
    xAxisInMsec = (1:48)/SamplingRateInKHZ;
    plots(5).unclassified = plot(xAxisInMsec,X(:,uidx),'color',[0.5 0.5 0.5]);
    set(plots(5).unclassified,'UserData',0);
%AZ20081210%    plot(X(:,uidx),'color',[0.5 0.5 0.5]); set(h.unclassified,'visible','off');
    My = max(max(X(:,uidx)));
    my = min(min(X(:,uidx)));
    ry = My-my;
    %AZ20081210: change xaxis to msec
    axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
    set(gca, 'Visible', 'on', 'Box', 'off', 'Color', get(gcf, 'Color'));
%AZ20081210%   axis([1 48 my-0.1*ry My+0.1*ry]);
    xlabel('msec');
    hold on;
    yl = get(gca,'ylim');
    
    set(plots(5).unclassified,'buttondownfcn','highlightTracePoint');
else %AZ20081211: Plot classified trains together in 'unclassified' box
    axes(h.unclassified);  hold on;
    col = {'r' 'b' 'g' 'm'};
    for i=1:length(sidx) % the sorted waveforms
        if(~isempty(sidx{i}))
            %AZ20081210: scale xaxis to msec (assuming sampling at 30kHz)
            SamplingRateInKHZ = 30;
            xAxisInMsec = (1:48)/SamplingRateInKHZ;
            plots(i).unclassified = plot(xAxisInMsec,X(:,sidx{i}),col{i},'Parent',h.unclassified);
            set(plots(i).unclassified,'UserData',i);
            
            % TODO: Uncomment once sortnev converted to non-guide
    %          set(plots(i).unit   ,'buttondownfcn',@(plots(i).unclassified(1))BringToTop(h));
    %          set(plots(i).unitavg,'buttondownfcn',@(plots(i).unclassified(1))BringToTop(h));
               set(plots(i).unclassified,'buttondownfcn','highlightTracePoint');
        else
            %AZ20090107: disable 'rest' buttons
            for j = 1:4
                eval(sprintf('set(h.unit%drest,''Enable'',''Off'')',j));
            end
        end
    end
    title('All Classified')
    hold off;
    %AZ20081211    
end