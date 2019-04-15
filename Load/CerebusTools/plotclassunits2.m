function [f1h,d] = plotclassunits2(f1h,X,V,sidx,uidx,col,gmm,C)
% PLOTCLASSUNITS
% Plots sorted units in cluster views and as waveforms
% change log:
%   created by Dario Ringach
%   2008-08 LB corrected some matlab complaints
%   2008-12 AZ modified unclassified xaxis to display in msec, assuming a
%       certain sampling rate. made prettier
%   2009-01 AZ classunit plots now plot all other (non-selected) waveforms
%       in light gray, behind selected waveforms.  d' calculated for each
%       nonempty unit, value is displayed above the given classunit plot

% global V sidx X Cm spoly gmm d prev;

% cluster views (top axis: p1 vs p2)
axes(f1h.p1p2); cla;
for i=1:length(sidx)
    hold on;
    if(~isempty(sidx{i}))
        f1h.plots(i).p1p2 = plot(V(sidx{i},1),V(sidx{i},2), [col{i} '.'],'markersize',1); 
        set(f1h.plots(i).p1p2,'UserData',i);
    end
%AZ20081211%    set(f1h.p1p2,'visible','off')
    xr = range(V(:,1)); yr = range(V(:,2));
    axis([min(V(:,1))-0.1*xr max(V(:,1))+0.15*xr min(V(:,2))-0.15*yr max(V(:,2))+0.1*yr])
    hold on, plot(0,0,'r.','markersize',7);
end

% cluster views (middle axis: p1 vs p3)
axes(f1h.p1p3); cla;
for i=1:length(sidx)
    hold on;
    if(~isempty(sidx))
        f1h.plots(i).p1p3 = plot(V(sidx{i},1),V(sidx{i},3), [col{i} '.'],'markersize',1);
        set(f1h.plots(i).p1p3,'UserData',i);
    end
%AZ20081211%    set(f1h.p1p2,'visible','off')
%     set(f1h.p1p3,'visible','off') %AZ20081211
    xr = range(V(:,1)); yr = range(V(:,3));
    axis([min(V(:,1))-0.1*xr max(V(:,1))+0.15*xr min(V(:,3))-0.15*yr max(V(:,3))+0.1*yr])
    hold on, plot(0,0,'r.','markersize',7);
end


% cluster views (bottom axis: p2 vs p3)
axes(f1h.p2p3); cla;
for i=1:length(sidx)
    hold on;
    if(~isempty(sidx))
        f1h.plots(i).p2p3 = plot(V(sidx{i},2),V(sidx{i},3), [col{i} '.'],'markersize',1);
        set(f1h.plots(i).p2p3,'UserData',i);
    end
%AZ20081211%    set(f1h.p1p2,'visible','off')
%     set(f1h.p2p3,'visible','off') %AZ20081211
    xr = range(V(:,3)); yr = range(V(:,3));
    axis([min(V(:,3))-0.1*xr max(V(:,2))+0.15*xr min(V(:,3))-0.15*yr max(V(:,3))+0.1*yr])
    hold on, plot(0,0,'r.','markersize',7);
end

yl = []; j = 1;
for i=1:length(sidx) % the sorted waveforms
    axes(eval(sprintf('f1h.unit%d',i)));
    cla; set(gca, 'Visible', 'off');
    
    %AZ20090106: nidx = indices of all other waveforms other than sidx{i}
    nidx = 1:length(X);
    nidx(sidx{i}) = [];
    nidx = nidx';
    if(~isempty(sidx{i}))
%AZ20081210: scale xaxis to msec (assuming sampling at 30kHz)
        SamplingRateInKHZ = 30;
        xAxisInMsec = (1:48)/SamplingRateInKHZ;
        %AZ20090106: plot (lightly) all other traces
        hold on;
        plot(xAxisInMsec,X(:,nidx   ),'Color', [0.8 0.8 0.8]);
        f1h.plots(i).unit = plot(xAxisInMsec,X(:,sidx{i}),col{i});
%AZ20081210%        plot(X(:,sidx{i}),col{i});

        %AZ20090212: plot previously sorted mean waveform
        if exist('prev','var') && sum(size(prev)) > 1 && j <= size(prev,2) && ...
                isfield(prev(j),'unit') && prev(j).unit.gmm.icell(i)
            f1h.plots(i).prevunitavg = plot(xAxisInMsec,prev(j).unit.prototype,'k-','LineWidth',1.5);
            j = j + 1;
        end

        %AZ20090127: plot mean waveform
        f1h.plots(i).unitavg = plot(xAxisInMsec,mean(X(:,sidx{i}),2)','w-','LineWidth',1.5);
        
        yl = [yl; get(gca,'ylim')]; %#ok<AGROW>
        My = max(max(X(:,sidx{i})));
        my = min(min(X(:,sidx{i})));
        ry = My-my;
%AZ20081210: change xaxis to msec
        axis([xAxisInMsec(1) xAxisInMsec(end) my-0.1*ry My+0.1*ry]);
        set(gca, 'Visible', 'on', 'Box', 'off', 'Color', get(gcf, 'Color'));
%AZ20081210%        axis([1 48 my-0.1*ry My+0.1*ry]);
%AZ20081210%        set(gca, 'Visible', 'on', 'Box', 'off', 'YColor', get(gcf, 'Color'), 'Color', get(gcf, 'Color'), 'XLim', [1 48]);

        %AZ20090107: calculate d', insert above chart
        d(i) = dprime(V,sidx{i},nidx);
        eval(sprintf('set(f1h.unit%dstatus,''String'',''d'''' = %0.2g'')',[i d(i)]));
    end
end

% Global scaling...

if(get(f1h.gs,'Value') && ~isempty(yl))
    yl = [min(yl(:,1)) max(yl(:,2))];
    for(i=1:length(sidx))
        axes(eval(sprintf('f1h.unit%d',i)));
        set(gca,'ylim',yl);
    end
end

%AZ20090204: Plot ellipses of gmm fits
k = get(f1h.nclusters,'Value')+1;
dim = 3:-1:1;
for j = dim
    n = fliplr(dim(find(dim~=j)));
    fit(j).gmm = gmdistribution(gmm.obj.mu(:,n),gmm.obj.Sigma(n,n,:),gmm.obj.PComponents);
    
    eval(sprintf('axes(f1h.p%dp%d)',n(1),n(2))); hold on;
%     scatter(V(:,n(1)),V(:,n(2)),10,'.')
    for i = 1:k
%         plot(gmm.obj.mu(i,n(1)),gmm.obj.mu(i,n(2)),col{i},'Marker','x',...
%             'MarkerSize',16,'LineWidth',5);
    
        %diagonalize the covariance matrix
        [u,lam]=eig(inv(fit(j).gmm.Sigma(:,:,i)));
        %generate a vector of angles from 0 to 2*pi
        theta=(0:.01:2*pi);
        %calculate the x component of the ellipsoid for all angles
        r(:,1)=fit(j).gmm.mu(i,1) + (1/sqrt(lam(1,1)))*u(1,1)*cos(theta) + ...
            (1/sqrt(lam(2,2)))*u(1,2)*sin(theta);
        r(:,3)=fit(j).gmm.mu(i,1) + (sqrt(5)/sqrt(lam(1,1)))*u(1,1)*cos(theta) + ...
            (sqrt(5)/sqrt(lam(2,2)))*u(1,2)*sin(theta);
        %calculate the y component of the ellipsoid for all angles
        r(:,2)=fit(j).gmm.mu(i,2) + (1/sqrt(lam(1,1)))*u(2,1)*cos(theta) + ...
            (1/sqrt(lam(2,2)))*u(2,2)*sin(theta);
        r(:,4)=fit(j).gmm.mu(i,2) + (sqrt(5)/sqrt(lam(1,1)))*u(2,1)*cos(theta) + ...
            (sqrt(5)/sqrt(lam(2,2)))*u(2,2)*sin(theta);
        %plot(x,y)
        f1h.plots(i).ellipse  = plot(r(:,1),r(:,2),col{i});
        f1h.plots(i).ellipse2 = plot(r(:,3),r(:,4),col{i});
    end
%     conth(j) = ezcontour(@(x,y)pdf(fit(j).gmm,[x y]),xlim,ylim);
    
%     Cn(j,:) = cluster(fit(j).gmm,V(:,n));
end
%AZ20090204: END Plot ellipses of gmm fits

drawnow;
% plots = plotunclassified(plots);

