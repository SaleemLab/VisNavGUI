% SAVE_DEBUG.M: Shows PCA points' re-classification during saving
% AZ 2009-01-29

% idx_subthreshold = setdiff(sidx{5},find(~iXraw));

i = 1;%spoly{1,1}{1}; %AZ20090204
a = floor(i/3)+1;
b = floor(i/2)+2;
figure(4); cla;
% set(gcf,'Renderer','zbuffer');

% centroid = gmm.obj.mu(:,1:2);

subplot(3,2,[1 3 5]); cla; hold on
% Plot classified
for i = k:-1:1
    line(V(sidx{i},a),V(sidx{i},b),'Color',col{i},...
        'LineStyle','none','Marker','.','MarkerSize',10);
end

% AZ20090710: PLOT sub-threshold separately
line(V(~iXraw,a),V(~iXraw,b),'Color','k',...
    'Marker','.','LineStyle','none','MarkerSize',10);
% Plot unclassified
line(V(setdiff(sidx{5},find(~iXraw)),a),V(setdiff(sidx{5},find(~iXraw)),b),...
    'Color',col{5},'Marker','.','LineStyle','none','MarkerSize',10);
axis tight;
title({'Classification applied to all traces in .nev file,';...
   [' electrode # ',...
    num2str(elec),', N = ',num2str(size(Xrealigned,2))]})

a = 0;
subplot(3,2,[2 4]); cla; hold on
if ~isempty(vertcat(sidx{1:4}))
   for i = k:-1:1
   %     if isempty(spoly{1,i}), i = i -1; break; end
   %     switch size(Xrealigned,1)
   %        case 48
             line(xAxisInMsec,Xrealigned(:,sidx{i}),'Color',col{i});
   %        case 28
   %           plot(xAxisInMsec,Xrealigned(7:34,sidx{i}),'Color',col{i});
   %     end
      % AZ20090317: mean waveform
       line(xAxisInMsec,mean(Xrealigned(:,sidx{i}),2)','Color','w',...
          'LineStyle','-','LineWidth',1.5);
       My = max(max(Xrealigned(:,sidx{i})));
       my = min(min(Xrealigned(:,sidx{i})));
       ry = My-my;
   %AZ20081210: change xaxis to msec
       axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
       a = a + size(sidx{i},1);
   end
end
title(['All Classified, N = ',num2str(size(vertcat(sidx{1:4}),1))])

subplot(3,2,6); cla; hold on
title(['All Unclassified, N = ',num2str(size(vertcat(setdiff(sidx{5},find(~iXraw))),1))])
if ~isempty(sidx{5})
%     switch size(Xrealigned,1)
%        case 48
          line(xAxisInMsec,Xrealigned(:,setdiff(sidx{5},find(~iXraw))),'Color',col{5});
%        case 28
%           plot(xAxisInMsec,Xrealigned(7:34,sidx{5}),'Color',[0.5 0.5 0.5]);
%     end
    My = max(max(Xrealigned(:,setdiff(sidx{5},find(~iXraw)))));
    my = min(min(Xrealigned(:,setdiff(sidx{5},find(~iXraw)))));
    ry = My-my;
    %AZ20081210: change xaxis to msec
    axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
end