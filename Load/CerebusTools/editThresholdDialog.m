function editThresholdDialog()
% editThresholdDialog.m

   % Get updated vars from base workspace
   varsToLoad = {'fh','X','Xraw','iXraw','threshold','d','nClusters','A','k',...
      'V','class_idx','calculated','xAxisInMsec','plotColors','scaling_params'};
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end

   spname = get(get(gco,'Parent'),'Tag');
   elec   = str2double(spname(5:end-4));

	tagname{1} = 'plot';
   tagname{2} = 'pca';

   answer = inputdlg(sprintf('New threshold (elec %2.0f):',elec),...
      sprintf('Change Threshold (elec %2.0f)',elec),1,{num2str(threshold(elec))});
   
   if ~isempty(answer{1})
      threshold(elec) = str2double(answer{1});
      [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},fh);
   end
   
   [class_idx,d,fh,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
      class_idx,A,iXraw,fh,calculated,V,elec);
   
   if ~isempty(X{elec})
      scaling_params.yAxisMax(elec) = 1.1*max(max(X{elec}(:,1:end)));
      scaling_params.yAxisMin(elec) = 1.1*min(min(X{elec}(:,1:end)));
   end
                  
   set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Refreshing view, please wait...'); drawnow;
   calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated,elec);
%    calculated = elecSubplotCreate(         1,A,          calculated,elec);
   if get(fh.expnev1.show.box.Num,'Value')
      plotElecNum(   1,tagname,A,d,elec);
   end
   if get(fh.expnev1.show.box.Traces,'Value')
      plotElecTraces(1,tagname,A,X,iXraw,Xraw,xAxisInMsec,plotColors,k,class_idx,...
         threshold,scaling_params.yAxisMin,scaling_params.yAxisMax,elec);
   end
   if get(fh.expnev1.show.box.D,'Value')
      plotElecD(     1,tagname,A,d);
   end
   set(findobj(gcf,'Tag','expnev1.txt.status'),'String',...
      sprintf('Elec %g threshold set to %g',elec,threshold(elec)));
   set(fh.expnev1.fig,'Pointer','arrow');
   drawnow;
   
   varsToSave = {'fh',...
      'X','d','class_idx','iXraw','V','threshold','scaling_params',...
      'calculated','k',};
   for thevar = varsToSave
      eval(sprintf('assignin(''base'',''%s'',%s)',thevar{:},thevar{:}));
   end

end
% END editThresholdDialog function