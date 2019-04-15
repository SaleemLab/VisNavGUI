function [Xrealigned,idx_tograph,idx_tograph_only,sidx_tograph,sidx,f1h] = ...
   sortnev2_gui_update(...
   Xrealigned,V,iXraw,sidx_tograph,gmm,prev,xAxisInMsec,col,sidx,d,threshold,f1h,MESSAGE)
% SORTNEV2_GUI_UPDATE.M: Updates PCA, Waveforms, and Controls in sortnev2 window

%       % RECALCULATE PCA
%       [U,S,V] = svds(Xrealigned,3);
%% GRAPH ONLY FIRST 3000 POINTS: Xrealigned
%if not coming from clickEllipse
if ~strcmp(MESSAGE,'Ellipses Cleared.') || ~strcmp(MESSAGE,'Clustering done!')
   idx_tograph = randperm(size(Xrealigned,2));
   idx_tograph = idx_tograph(1:min(3000,size(idx_tograph,2)));
   % AZ20090707: separate points
   idx_tograph_only = intersect(idx_tograph,find(~iXraw));
   idx_tograph = setdiff(idx_tograph,idx_tograph_only);
end

% unclassified
sidx{5} = find(iXraw);
%AZ20081211: %Make sidx{5} have contain all spike train indices, then take out
%sorted ones
sidx{5} = setdiff(sidx{5},union(union(sidx{1},sidx{2}),union(sidx{3},sidx{4})));

if ~isempty(cat(1,sidx_tograph{:})) %cluster_gmm3sd
   sidx_tograph{5} = sidx{5}(ismember(sidx{5},idx_tograph));
else                         %loadelec
   for i = 1:5
      sidx_tograph{i} = intersect(sidx{i},idx_tograph)';
   end
end
sidx_tograph_only = intersect(find(~iXraw),idx_tograph_only)';

f1h = plotPCA(        f1h,           V,sidx_tograph,gmm,                 col,sidx,sidx_tograph_only);
f1h = plotWaveforms(  f1h,Xrealigned,  sidx_tograph,    prev,xAxisInMsec,col,sidx,sidx_tograph_only,threshold);
%    f1h = plotsortedunits(f1h,Xrealigned,  sidx_tograph,    prev,xAxisInMsec,col,sidx,d);

f1h = sortnev2_controls_update(true,f1h,d,col,sidx_tograph);

%AZ20081217
set(f1h.txt.status,'String',MESSAGE);
set(f1h.fig.f1,'Pointer','arrow');
drawnow;

end