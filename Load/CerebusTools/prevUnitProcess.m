function [Xrealigned,V,xAxisInMsec,sidx,d,threshold,fh,iXraw,gmm,prev,S,U,C,k] = ...
   prevUnitProcess( Xrealigned,gmm,prev,xAxisInMsec,Xraw,iXraw,nevopen_outcome,...
          SamplingRateInKHZ,already_sorted_to_load,animaltype,sorted_series,...
          animal_spikesdir,prevSavedUnits,iseries,iexp,threshold,fh,nSamplesPerPacket)

global DIRS

if ~isempty(already_sorted_to_load)
   d = struct;
   C    = cell(size(prev,2),1);
   V    = cell(size(prev,2),1);
   S    = []; %AZ20091006
   U    = []; %AZ20091006
   sidx = cell(size(prev,2),5);
   
   % Read in parameters from all other saved units (INTO MEMORY ONLY)
   for i = 1:size(prevSavedUnits,1)
      unit = load([DIRS.spikes filesep        animaltype filesep ...
                        num2str(sorted_series(i)) filesep ...
                   animal_spikesdir{prevSavedUnits(i)}]);
      prev(i).unit = unit.unit;
      if isfield(prev(i).unit,'gmm')
         %AZ20090502: Read in old GMM data, but interpret with 3SD cutoff method
         if isfield(prev(i).unit.gmm,'obj')
            if ~isfield(prev(i).unit.gmm,'dprime')
               % TODO: FIX this
               [prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,sorted_series] = ...
                  eliminateUnit(prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,...
                  sorted_series,i);
               return;
            end
            prev(i).unit.gmm.mu    = prev(i).unit.gmm.obj.mu;
            prev(i).unit.gmm.Sigma = prev(i).unit.gmm.obj.Sigma;
%             delete(prev(i).unit.gmm.obj);
         end
         if prev(i).unit.iexp == iexp
            % If this unit is for the same exp as the one we are analyzing,
            % populate the replace_dprime field with the unit's dprime, so
            % we know what we should strive to beat while sorting
            prev(i).replace_dprime = prev(i).unit.gmm.dprime;
         end
         prev(i).loaded = prev(i).unit.gmm.icell;
      else % AZ 2009-10-06
         % TODO: FIX this
         [prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,sorted_series] = ...
            eliminateUnit(prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,...
            sorted_series,i);
         error('Old style unit ARCHIVED.  Restart sortnev.')
      end
   end

   % If any replace_dprime field exists, make it the max across all units found
   for i = 1:size(prevSavedUnits,1)
      if ~isempty(prev) && isfield(prev(i), 'replace_dprime')
         prev(i).replace_dprime = max(vertcat(prev.replace_dprime));
      end
   end
   
%    % AZ 2009-10-06: If no gmm field, skip
%    a = false(size(prevSavedUnits,1),1);
%    for i = 1:size(prevSavedUnits,1)
%       a(i) = isfield(prev(i).unit,'gmm');
%    end
%    prevSavedUnits = prevSavedUnits(a);
%    prev           = prev(a);
%    if isempty(prevSavedUnits)
%       return;
%    end
   
   UnitsToKeep = true(size(prev));
   
   % Iterate through each sorted unit, calculate d'
   for unitIDX = 1:size(prev,2)

      probThreshold = zeros(1,5); P = zeros(size(Xrealigned,2),5);
      % RESIZE
      sizediff = size(Xrealigned,1) - size(prev(unitIDX).unit.gmm.U,1);
      if sizediff ~= -20  %% TODO: FIX: don't discard units length 28
         switch sizediff
            case 0
               % do nothing
               V{unitIDX} = ( prev(unitIDX).unit.gmm.S\...
                             (prev(unitIDX).unit.gmm.U'*Xrealigned        ) )';
            case 20
               % take only samples (7:34) (t = -0.1:1/30:0.8)
               V{unitIDX} = ( prev(unitIDX).unit.gmm.S\...
                             (prev(unitIDX).unit.gmm.U'*Xrealigned(7:end-14,:)) )';
            otherwise
               warning('sortnev2:prevUnitProcess:BadSizeUprev','Weird size of loaded PCA''s');
         end
         
         k = size(prev(unitIDX).unit.gmm.icell,2);
         
          % AZ20090502: Clustering method replaced with new 3SD cutoff method
    %       C{unitIDX} = cluster(prev(unitIDX).unit.gmm.obj,V{unitIDX});

         for i = 1:k
            probThreshold(i) = mvnpdf(prev(unitIDX).unit.gmm.mu(i,:),      ...
                                      prev(unitIDX).unit.gmm.mu(i,:),      ...
                                      prev(unitIDX).unit.gmm.Sigma(:,:,i));
            P(iXraw,i)       = mvnpdf(   V{unitIDX}(iXraw,:),              ...
                                      prev(unitIDX).unit.gmm.mu(i,:),      ...
                                      prev(unitIDX).unit.gmm.Sigma(:,:,i));
         end
          % Threshold at 3 SD's (of largest distribution): ~99.7%
          probThreshold = (1-erf(3/sqrt(2)))*max(probThreshold);
          % C = index of maximum probability
          [a,C{unitIDX}] = max(P,[],2);
          %indices where all probs for all classes are below (prob)
          %threshold, and above (spike) threshold
          sidx{unitIDX,5} = find(sum(P(:,1:k)<probThreshold,2)==k);
          sidx{unitIDX,5} = intersect(sidx{unitIDX,5},find(iXraw));
          
          clear a P;

          for i = 1:4
             sidx{unitIDX,i} = intersect(find(C{unitIDX}==i),find(iXraw));
             sidx{unitIDX,i} = setdiff(sidx{unitIDX,i},sidx{unitIDX,5});
          end
          for i = 1:4
             if ~isempty(sidx{unitIDX,i})
                d.candidate{unitIDX}(i) = dprime(V{unitIDX},sidx{unitIDX,i},...
                   vertcat(sidx{unitIDX,setdiff(1:4,i)}));
                if isnan(d.candidate{unitIDX}(i)), d.candidate{unitIDX}(i) = []; end
             end
          end

          try   d.max(unitIDX) = max(d.candidate{unitIDX});
          catch d.max(unitIDX) = 0;
          end
      else
         % TODO: fix this
         k = 2;
         [prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,sorted_series] = ...
            eliminateUnit(prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,...
            sorted_series,unitIDX);
         UnitsToKeep(unitIDX) = ~isempty(prev(unitIDX).unit);
      end
   end

   % unload previously loaded units that have been archived
   prev                   = prev(                  UnitsToKeep  );
   already_sorted_to_load = already_sorted_to_load(UnitsToKeep,:);
   prevSavedUnits         = prevSavedUnits(        UnitsToKeep  );
   sorted_series          = sorted_series(         UnitsToKeep  );
   
   % Find saved unit with highest d'
   try
      [a,i] = max(d.max);
   catch
      error(['This is a particularly nasty bug I haven''t found a way to fix.  ',...
         'If you restart sortnev, it should work.']);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECTION GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Generate selection text
   fdialogh.cellSortedUnits    = cell(1,size(prev,2)+1);

   unitIDX = 1;
   while unitIDX <= size(prev,2)
      fdialogh.cellSortedUnits{unitIDX} = sprintf(                          ...
         's%02u e%02u elec:%04u cell:%03u  --  d'' = %2.4f %2.4f %2.4f %2.4f',...
         [already_sorted_to_load(unitIDX,1:4),d.candidate{unitIDX}]              );
      if iseries == already_sorted_to_load(unitIDX,2)
         fdialogh.cellSortedUnits{unitIDX} = ['**' ...
            fdialogh.cellSortedUnits{unitIDX} '**'];
      end
      unitIDX = unitIDX + 1;
   end
   
   % find last line with series < current series
   doNotLoadIDX = find(already_sorted_to_load(:,2) < iexp,1,'last')+1;
   if isempty(doNotLoadIDX)
      doNotLoadIDX = 1;
   end
%    fdialogh.cellSortedUnits{end+1} = [];
   fdialogh.cellSortedUnits([1:doNotLoadIDX-1,doNotLoadIDX+1:end]) = fdialogh.cellSortedUnits(1:end-1);
   fdialogh.cellSortedUnits{doNotLoadIDX} = 'SELECT ME TO NOT LOAD ANY SAVED PARAMETERS';

   % Draw dialog figure
   set(0,'defaultFigureColor',get(0,'defaultUicontrolBackgroundColor'));

   fdialog = figure('MenuBar','none','Units','Pixels','Position',[200 200 520 225],...
      'NumberTitle','off','Name',sprintf(['sortnev2: Load saved parameters (',...
      animaltype,'_s%02u_e%02u)?'],[iseries iexp]),...
      'CloseRequestFcn',@closedialog,'KeyPressFcn',@closeOnEnterKey);

   fdialogh.listSortedUnits  = uicontrol('Style','listbox','Value',i+1,...
      'Units','Pixels','Tag','listSortedUnits','Position',[20 65 480 120],...
      'String',fdialogh.cellSortedUnits,'Enable','on','FontName','FixedWidth',...
      'KeyPressFcn',@closeOnEnterKey);
   fdialogh.txtInst = uicontrol('Parent',fdialog,'Units','Pixels',...
     'HorizontalAlignment','left','KeyPressFcn',@closeOnEnterKey,...
     'String',{'There are units from other files saved for this animal and electrode.';...
     'Would you like to load their sorting parameters?';...
     'Default selection: The saved unit yielding the highest d'' for this expt.'},...
     'Style','text','FontWeight','bold','Tag','txtInst','Position',[20 8 480 45]);

   fdialogh.pbLoad = uicontrol('Parent',fdialog,'Units','Pixels',...
     'HorizontalAlignment','left',...
     'String','Load','Style','pushbutton','Tag','pbLoad',...
     'Position',[420 35 80 20],'Callback',@closedialog,'KeyPressFcn',@closeOnEnterKey);
   fdialogh.pbNoLoad = uicontrol('Parent',fdialog,'Units','Pixels',...
     'HorizontalAlignment','left',...
     'String','Do Not Load','Style','pushbutton','Tag','pbNoLoad',...
     'Position',[420 10 80 20],'Callback',@closedialog,'KeyPressFcn',@closeOnEnterKey);
  
   fdialogh.txtWarn = uicontrol('Parent',fdialog,'Units','Pixels',...
     'HorizontalAlignment','left','KeyPressFcn',@closeOnEnterKey,...
     'String',{'This window automatically closes after 5 seconds.';...
     'Press any key (except return or escape) or click inside the window to keep it open.'},...
     'Style','text','FontWeight','bold','ForegroundColor',[1 0 0],...
     'Tag','txtWarn','Position',[20 190 480 30]);

   % updateEditPath;   % Update text box with RootDir\SubDir

   uiwait(fdialog,5);  % Wait for user input
   
   if ishandle(fdialog)
      if (~isempty(get(fdialog,'CurrentCharacter')) || ~isempty(gco(fdialog)))
         set(fdialogh.txtWarn,'Visible','off'); drawnow;
         uiwait(fdialog);  % Wait for user input
      else
         closedialog(findobj('Tag','pbLoad'));
      end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SELECTION GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if unitIDX ~= 0 && unitIDX ~= doNotLoadIDX
      if unitIDX > doNotLoadIDX
         unitIDX = unitIDX - 1;
      end
      
      % LOAD SELECTED PARAMETERS
      sidx      = sidx(unitIDX,:);
      d         = d.candidate{unitIDX};
      V         = V{unitIDX};
      C         = C{unitIDX};
      S         = prev(unitIDX).unit.gmm.S;
      U         = prev(unitIDX).unit.gmm.U;
      gmm.mu    = prev(unitIDX).unit.gmm.mu;
      gmm.Sigma = prev(unitIDX).unit.gmm.Sigma;
      
      % Set k, based on unit to load
      k = size(prev(i).unit.gmm.icell,2);
      set(fh.expnev1.cluster.popup.n,'Value',k);
      
      % Change threshold if the loaded unit has a different threshold
      if isfield(prev(unitIDX).unit.gmm, 'threshold') && ...
                 prev(unitIDX).unit.gmm.threshold ~= threshold
         set(fh.expnev1.edit.thresh,'String',num2str(prev(unitIDX).unit.gmm.threshold))
         [Xrealigned,iXraw] = changeThreshold(threshold,Xraw,fh,nevopen_outcome);
      end

      % RESIZE Xrealigned
      % take only samples (7:34) (t = -0.1:1/30:0.8)
      Xrealigned = Xrealigned(max(1,size(Xrealigned,1)-41):...
                              min(  size(Xrealigned,1),nSamplesPerPacket-14),:);
      xAxisInMsec = ((7:nSamplesPerPacket-14)-10)/SamplingRateInKHZ; % center about thresholding point (10)

      % Text feedback
      disp(['Parameters loaded from ', animaltype filesep ...
                        num2str(sorted_series(unitIDX)) filesep ...
                   animal_spikesdir{prevSavedUnits(unitIDX)}]);
      if isfield(prev(unitIDX).unit.gmm,'dprime')
         disp(['-->d'' = ', num2str(max(prev(unitIDX).unit.gmm.dprime))]);
      end
   else % Don't load any parameters
      % Must clear temporarily loaded parameters from above;
      % MATLAB gets angry if output variables are not initialized
      C    = [];  d    = [];  gmm  = [];  sidx = cell(1,5);
      
      % RESIZE
      % take only samples (7:34) (t = -0.1:1/30:0.8)
      Xrealigned = Xrealigned(max(1,size(Xrealigned,1)-41):...
                              min(  size(Xrealigned,1),nSamplesPerPacket-14),:);
      xAxisInMsec = ((7:nSamplesPerPacket-14)-10)/SamplingRateInKHZ; % center about thresholding point (10)

      [U,S,V] = svds(Xrealigned,3);
   end
else % else if there are no saved units to load
   % MATLAB gets angry if output variables are not initialized
   C    = [];  d    = [];  gmm  = [];  sidx = cell(1,5);
   % RESIZE
   % take only samples (7:34) (t = -0.1:1/30:0.8)
   Xrealigned = Xrealigned(max(1,size(Xrealigned,1)-41):...
                           min(  size(Xrealigned,1),nSamplesPerPacket-14),:);
   xAxisInMsec = ((7:nSamplesPerPacket-14)-10)/SamplingRateInKHZ; % center about thresholding point (10)

   [U,S,V] = svds(Xrealigned,3);
end


function closeOnEnterKey(hObject,eventdata)
    if strcmp(eventdata.Key,'return')
       % Pretend we clicked 'Load'        in fdialog
       closedialog(findobj('Tag','pbLoad')  ,eventdata);
    elseif strcmp(eventdata.Key,'escape')
       % Pretend we clicked 'Do Not Load' in fdialog
       closedialog(findobj('Tag','pbNoLoad'),eventdata);
    elseif strcmp(eventdata.Key,'downarrow')
       if get(fdialogh.listSortedUnits,'Value') < size(prev,2)+1
          set(fdialogh.listSortedUnits,'Value',get(fdialogh.listSortedUnits,'Value')+1)
       end
    elseif strcmp(eventdata.Key,'uparrow')
       if get(fdialogh.listSortedUnits,'Value') > 1
          set(fdialogh.listSortedUnits,'Value',get(fdialogh.listSortedUnits,'Value')-1)
       end
    end
end

function closedialog(hObject,eventdata)
   set(fdialog,'Pointer','watch'); drawnow;
   if hObject == fdialogh.pbLoad % Clicked 'Load', read selected parameters
      unitIDX = get(fdialogh.listSortedUnits,'Value');
   else  % Didn't click 'OK', CLEAR selected path
      unitIDX = 0;
   end
   delete(fdialog);
end

end