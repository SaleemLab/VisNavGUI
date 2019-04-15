function f1h = plotunclassified2(f1h,X,V,sidx,uidx)
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

% global X V sidx uidx lowidx highidx;
% global f1h

% cluster views (top axis: p1 vs p2)
% h = guidata(sortnev);

%AZ20081210: scale xaxis to msec (assuming sampling at 30kHz)
SamplingRateInKHZ = 30;
xAxisInMsec = (1:48)/SamplingRateInKHZ;
col = {'r' 'b' 'g' 'm'};

% Clear existing highlights
if isfield(f1h,'unclassified_high') && ~isempty(f1h.unclassified_high) && ishandle(f1h.unclassified_high)
    delete(f1h.unclassified_high);
    delete(f1h.unclassified_high2);

    refresh;
end
if isfield(f1h,'p1p2_high2') && ~isempty(f1h.p1p2_high2) && ishandle(f1h.p1p2_high2)
    delete(f1h.p1p2_high2);
    delete(f1h.p1p2_high1);
    delete(f1h.p1p3_high2);
    delete(f1h.p1p3_high1);
    delete(f1h.p2p3_high2);
    delete(f1h.p2p3_high1);
    
    refresh;
end

% unclassified
uidx = 1:size(V,1);
%AZ20081211: %Make uidx have contain all spike train indices, then take out
             %sorted ones
uidx = setdiff(uidx,union(union(sidx{1},sidx{2}),union(sidx{3},sidx{4})));

if ~isempty(uidx)
   axes(f1h.p1p2); cla; hold on;
   for a = 1:length(uidx)
       f1h.plots(5).p1p2(a) = plot(V(uidx(a),1),V(uidx(a),2), '.', ...
           'markersize', 3, 'Color', [0.5 0.5 0.5], 'UserData', [a 5], ...
           'buttondownfcn',@highlightPointTrace,'Parent',f1h.p1p2);
   end
   set(gca,'XColor',[0.4 0 0],'YColor',[0 0.4 0])
   xr = range(V(:,1)); yr = range(V(:,2));
   axis([min(V(:,1))-0.1*xr max(V(:,1))+0.15*xr min(V(:,2))-0.15*yr max(V(:,2))+0.1*yr])
   set(gca, 'Visible', 'on', 'Box', 'off', 'Color', get(gcf, 'Color')); %AZ20081211

   axes(f1h.p1p3); cla; hold on;
   for a = 1:length(uidx)
       f1h.plots(5).p1p3(a) = plot(V(uidx(a),1),V(uidx(a),3), '.', ...
           'markersize', 3, 'Color', [0.5 0.5 0.5], 'UserData', [a 5], ...
           'buttondownfcn',@highlightPointTrace,'Parent',f1h.p1p3);
   end
   set(gca,'XColor',[0.4 0 0],'YColor',[0 0 0.4])
   xr = range(V(:,1)); yr = range(V(:,3));
   axis([min(V(:,1))-0.1*xr max(V(:,1))+0.15*xr min(V(:,3))-0.15*yr max(V(:,3))+0.1*yr])
   set(gca, 'Visible', 'on', 'Box', 'off', 'Color', get(gcf, 'Color')); %AZ20081211
   
   axes(f1h.p2p3); cla; hold on;
   for a = 1:length(uidx)
       f1h.plots(5).p2p3(a) = plot(V(uidx(a),2),V(uidx(a),3), '.', ...
           'markersize', 3, 'Color', [0.5 0.5 0.5], 'UserData', [a 5], ...
           'buttondownfcn',@highlightPointTrace,'Parent',f1h.p2p3);
   end
   set(gca,'XColor',[0 0.4 0],'YColor',[0 0 0.4])
   xr = range(V(:,2)); yr = range(V(:,3));
   axis([min(V(:,2))-0.1*xr max(V(:,2))+0.15*xr min(V(:,3))-0.15*yr max(V(:,3))+0.1*yr])
   set(gca, 'Visible', 'on', 'Box', 'off', 'Color', get(gcf, 'Color')); %AZ20081211
   
end

axes(f1h.unclassified);
cla; hold on;
uidx = 1:size(X,2);
uidx = setdiff(uidx,union(union(sidx{1},sidx{2}),union(sidx{3},sidx{4})));

if(~isempty(uidx))
    title('Unclassified')
    for a = 1:length(uidx)
        f1h.plots(5).unclassified(a) = plot(xAxisInMsec,X(:,uidx(a)),...
           'color',[0.5 0.5 0.5], 'UserData', [a 5], ...
           'buttondownfcn',@highlightPointTrace,'Parent',f1h.unclassified);
    end
    My = max(max(X(:,uidx)));
    my = min(min(X(:,uidx)));
    ry = My-my;
    %AZ20081210: change xaxis to msec
    axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
    set(gca, 'Visible', 'on', 'Box', 'off', 'Color', get(gcf, 'Color'));
    xlabel('msec');
    hold on;
    yl = get(gca,'ylim');
    
else %AZ20081211: Plot classified trains together in 'unclassified' box
    for i=1:length(sidx) % the sorted waveforms
        if(~isempty(sidx{i}))
            for a = 1:length(sidx{i})
                f1h.plots(i).unclassified(a) = plot(xAxisInMsec,X(:,sidx{i}(a)),...
                   'color',col{i}, 'UserData', [a i], ...
                   'buttondownfcn',@highlightPointTrace,'Parent',f1h.unclassified);
            end
        else
            %AZ20090107: disable 'rest' buttons
            for j = 1:4
                eval(sprintf('set(f1h.unit%drest,''Enable'',''Off'')',j));
            end
        end
    end
    title('All Classified')
    hold off;
    %AZ20081211    
end