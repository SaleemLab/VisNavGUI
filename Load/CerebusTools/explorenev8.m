function explorenev8(varargin)
% PLEASE PLEASE WRITE A HELP TEXT
% No longer nested functions
%
%
% 2010-06-14 MC replaced ~ in function arguments with ignoreme to avoid crash

%% Begin Main Figure Main Function
% toggle: 'Menu'/'Toolbar' = 'none'/'figure'
if exist('fh','var') && isfield(fh,'fig') && isfield(fh.fig,'expnev1') ...
      && ishandle(fh.expnev1.fig)
   set(0,'CurrentFigure',fh.expnev1.fig); clf;
% Condition above assumes fh is a global, below doesn't
elseif ~isempty(findobj('Tag','expnev1.fig'))
   fh.expnev1.fig = findobj('Tag','expnev1.fig');
   set(0,'CurrentFigure',fh.expnev1.fig); clf;
else
   fh.expnev1.fig = figure('Position',[520 380 672 420],'Name','explorenev8',...
      'Toolbar','none','MenuBar','none','NumberTitle','off','Tag','expnev1.fig');
end
% Maximize
% set(gcf,'units','normalized','outerposition',[0 0.05 1 0.95]);
set(gcf,'units','normalized','outerposition',[0.05 0.45 0.5 0.5]);
% Create a variable to store the handle of the new GUI figure
% f2 becomes non-empty when the new figure is created.
f2 = []; f2h = [];

% Initialize variables
nevopen_outcome   = [];

%% Initialize GUI objects
% Create the 'load' pushbutton
fh.expnev1.button.load = uicontrol(fh.expnev1.fig,'Style', 'pushbutton','String','Load',...
    'Units','normalized','Tag','loadbutton','Position',[0.82 0.93 0.06 0.05],...
    'callback',@load_button_press);
% Create the 'status' text area
fh.expnev1.txt.status = uicontrol('Style','text','Tag','expnev1.txt.status','Units','normalized',...
    'Position',[0.8 0 0.2 0.075],'BackgroundColor',get(fh.expnev1.fig,'Color'));
% Create the 'Scaling' panel, containing scaling buttons
fh.expnev1.scaling.panel = uibuttongroup('Title','Scaling','Tag','expnev1.scaling.panel',...
    'Units','normalized','Position',[0.8 0.075 0.2 0.15],...
    'BackgroundColor',get(fh.expnev1.fig,'Color'),'SelectionChangeFcn',@scaling_button_press);
% Create the 'Global' scaling button
fh.expnev1.scaling.button.global = uicontrol('Parent',fh.expnev1.scaling.panel,'String','Global',...
    'Style','radiobutton','Tag','expnev1.scaling.button.global','Value',0,...
    'Units','pixels','BackgroundColor',get(fh.expnev1.fig,'Color'),'Position',[11 24 87 23]);
% Create the 'Fit'    scaling button
fh.expnev1.scaling.button.fit    = uicontrol('Parent',fh.expnev1.scaling.panel,'String','Fit',...
    'Style','radiobutton','Tag','expnev1.scaling.button.fit'   ,'Value',1,...
    'Units','pixels','BackgroundColor',get(fh.expnev1.fig,'Color'),'Position',[11 4 87 23]);
% Create the 'Show' panel, containing show boxes
fh.expnev1.show.panel = uipanel('Title','Show','Tag','expnev1.show.panel',...%'Visible','off',...
    'Units','normalized','Position',[0.8 0.225 0.2 0.15],...
    'BackgroundColor',get(fh.expnev1.fig,'Color'));
fh.expnev1.show.box.D = uicontrol('Parent',fh.expnev1.show.panel,'String','d''',...
    'Style','checkbox','Tag','expnev1.show.box.D','Value',1,'Enable','off',...
    'Units','pixels','BackgroundColor',get(fh.expnev1.fig,'Color'),'Position',[11 24 40 23],...
    'callback',@showorhide);
fh.expnev1.show.box.Num = uicontrol('Parent',fh.expnev1.show.panel,'String','Elec #',...
    'Style','checkbox','Tag','expnev1.show.box.Num','Value',1,...
    'Units','pixels','BackgroundColor',get(fh.expnev1.fig,'Color'),'Position',[11 4 60 23],...
    'callback',@showorhide);
fh.expnev1.show.box.Axes = uicontrol('Parent',fh.expnev1.show.panel,'String','Axes',...
    'Style','checkbox','Tag','expnev1.show.box.Axes','Value',1,...
    'Units','pixels','BackgroundColor',get(fh.expnev1.fig,'Color'),'Position',[71 24 60 23],...
    'callback',@showorhide);
fh.expnev1.show.box.Traces = uicontrol('Parent',fh.expnev1.show.panel,'String','Traces',...
    'Style','checkbox','Tag','expnev1.show.box.Traces','Value',1,...
    'Units','pixels','BackgroundColor',get(fh.expnev1.fig,'Color'),'Position',[71 4 60 23],...
    'callback',@showorhide);
% Create the 'Clustering' panel
fh.expnev1.cluster.panel = uipanel('Title','Clustering','Tag','expnev1.cluster.panel',...
    'Units','normalized','Position',[0.8 0.375 0.2 0.15],...%'Visible','off',...
    'BackgroundColor',get(fh.expnev1.fig,'Color'));
% fh.textClusters = uicontrol('Parent',fh.expnev1.cluster.panel,'Style','text',...
%     'Tag','textClusters',...
%     'Units','normalized','String','Clusters','Position',[0.43 0.35 0.5 0.5],...
%     'HorizontalAlignment','left','BackgroundColor',get(fh.expnev1.fig,'Color'));
fh.expnev1.cluster.button = uicontrol('Parent',fh.expnev1.cluster.panel,'Style','pushbutton',...
    'Tag','expnev1.cluster.button',...
    'Units','normalized','String','Cluster!','Position',[0.525 0.525 0.4 0.45],...
    'HorizontalAlignment','left','callback',@ClusterButton_Callback);
fh.expnev1.cluster.popup.n = uicontrol('Parent',fh.expnev1.cluster.panel,'Style','popupmenu',...
    'Tag','expnev1.cluster.popup.n',...%'BackgroundColor',get(fh.expnev1.fig,'Color'),...
    'String',{'1';'2';'3';'4';'best'},'Units','normalized',...
    'Position',[0.08 0.525 0.375 0.45],'callback',@ClusterPopupN_Callback);
fh.expnev1.cluster.popup.nwaves = uicontrol('Parent',fh.expnev1.cluster.panel,'Style',...
    'popupmenu','Tag','expnev1.cluster.popup.nwaves',...%'BackgroundColor',get(fh.expnev1.fig,'Color'),...
    'String',{'100';'500';'1000';'2500'},'Units','normalized',...
    'Position',[0.08 0.10 0.375 0.45],'callback',@ClusterPopupN_Callback);
fh.expnev1.cluster.text.nwaves = uicontrol('Parent',fh.expnev1.cluster.panel,'Style',...
    'text','Tag','expnev1.cluster.text.nwaves','Units','normalized',...
    'Position',[0.525 0 0.4 0.45],'BackgroundColor',get(fh.expnev1.fig,'Color'),...
    'String','# Waves');

fh.expnev1.PCABox = uicontrol('String','PCA (New Window)',...
    'Style','checkbox','Tag','expnev1.PCABox','Value',0,'Enable','on',...
    'Units','normalized','BackgroundColor',get(fh.expnev1.fig,'Color'),...
    'Position',[0.81 0.525 0.2 0.05],...
    'callback',@pcafigure);
% fh.SortnevLinksBox = uicontrol('String','Sortnev links',...
%     'Style','checkbox','Tag','expnev1.PCABox','Value',1,'Enable','on',...
%     'Units','normalized','BackgroundColor',get(fh.expnev1.fig,'Color'),...
%     'Position',[0.81 0.575 0.2 0.05],...
%     'callback',@pcafigure);
 
% AZ20090220: Thresholding tool
% Create the 'Threshold' panel
% fh.ThreshPanel = uipanel('Title','Threshold','Tag','ThreshPanel',...%'Visible','off',...
%     'Units','normalized','Position',[0.8 0.625 0.2 0.15],...
%     'BackgroundColor',get(fh.expnev1.fig,'Color'));
fh.expnev1.txt.thresh = uicontrol('Parent',fh.expnev1.fig,'Units','normalized',...
   'HorizontalAlignment','left','BackgroundColor',get(fh.expnev1.fig,'Color'),...
   'ListboxTop',0,'String','Threshold: +/-','Style','text','Tag','expnev1.txt.thresh',...
   'Position',[0.810 0.700 0.140 0.026]);
fh.expnev1.edit.thresh = uicontrol('Parent',fh.expnev1.fig,'Units','normalized',...
   'BackgroundColor',[1 1 1],'FontSize',9,...
   'Callback',@editThreshold_Callback,...%'KeyPressFcn',@shortcutKeys,...
   'ListboxTop',0,'String','32','Style','edit','Tag','expnev1.edit.thresh',...
   'Position',[0.920 0.696 0.040 0.032]);
fh.expnev1.pshbtn.save998 = uicontrol('Parent',fh.expnev1.fig,'Style','pushbutton',...
   'String','Save 998','Units','normalized','Tag','expnev1.pshbtn.save998',...
   'Callback',@save998_Callback,...
   'Position',[0.890 0.930 0.085 0.050]);

%% TODO: figure out a way to avoid this quoted code absurdity
% FOR subplots created in elecSubplotCreate
% Define a context menu; it is not attached to anything
fh.expnev1.menu.subplotfcns(1) = uicontextmenu('Tag','expnev1.menu.subplotfcns(1)');
% Define the context menu items and install their callbacks
fh.expnev1.menu.subplotfcns(2) = uimenu(fh.expnev1.menu.subplotfcns(1),  ...
   'Label','Open in sortnev2','Tag','expnev1.menu.subplotfcns(2)', ...
   'Callback', ...
   ['spname = get(get(gco,''Parent''),''Tag'');',...
    'elec = str2double(spname(5:end-4));',...
    'threshold(elec) = max(get(findobj(''Tag'',sprintf(''elec%dplotThresh'',elec)),''YData''));',...
    'fname = regexp(get(gcf,''Name''),''[\s,]'',''split'');',...
    'fname = fname(~strcmp(fname,''''));',...
    'sortnev2(fname{2},str2double(fname{end-1}),str2double(fname{end}),elec,threshold(elec));']);
fh.expnev1.menu.subplotfcns(3) = uimenu(fh.expnev1.menu.subplotfcns(1),...
   'Label','Change Threshold','Tag','expnev1.menu.subplotfcns(3)',...
   'Callback','editThresholdDialog();');

%    ['spname=get(get(gco,''Parent''),''Tag'');',...
%     'elec = str2double(spname(5:6));',...
%     ' answer = inputdlg(sprintf(''New threshold: for elec %2.0f'',elec),sprintf(''Change Threshold (elec %f)'',elec),1,{''32''});',...
%     ' if ~isempty(answer{1});',...
%     '    threshold(elec) = str2double(answer{1});',...
%     '    [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},fh);',...
%     ' end']);
% fh.expnev1.menu.subplotfcns(4) = uimenu(fh.expnev1.menu.subplotfcns(1),...
%    ...%'Callback', @clickEllipse_Callback,...
%    'Label','fprintf(''Electrode #%f'',elec)','Tag','expnev1.menu.subplotfcns(4)');


% fh.sidebar = uipanel('Title','Controls','Tag','expnev1.cluster.panel',...%'Visible','off',...
%     'Units','normalized','Position',[0.8 0.375 0.2 0.095],...
%     'BackgroundColor',get(fh.expnev1.fig,'Color'));


%% Initialize variables
SetDefaultDirs;

nClusters = get(fh.expnev1.cluster.popup.n,'String');
nClusters = nClusters{get(fh.expnev1.cluster.popup.n,'Value')};          % nClusters = # of classes/neurons
nClusters = str2double(nClusters);
try k = nClusters*ones(96,1); catch, k = zeros(96,1); end

nwaves = get(fh.expnev1.cluster.popup.nwaves,'String');
nwaves = str2double(nwaves{get(fh.expnev1.cluster.popup.nwaves,'Value')});% # of waves
% SamplingRateInKHZ = 30;
% xAxisInMsec = (1:48)/SamplingRateInKHZ;

threshold = str2double(get(fh.expnev1.edit.thresh,'String'))*ones(96,1);

% class_idx_temp    = cell(1,4);
% d_temp    = []; k_temp    = [];
% deltat    = [];

% Keep track of what was calculated during a given clustering session...
calculated(1).D      = 0;
calculated(1).Num    = 0;
calculated(1).Traces = 0;
calculated(1).Axes   = 0;   % make sure this is always 0
calculated(2).Num    = 0;
calculated(2).PCA    = 0;
calculated(2).Axes   = 0;   % make sure this is always 0

% Plot parameters: X,Y = Entire array, x,y = inidividual 
Xoffset = 0.03;
Yoffset = 0.10;
xoffset = 0.077;
yoffset = 0.097;
xsize   = 0.076;
ysize   = 0.096;
plotColors = {[0 0 1] [0 1 0] [1 1 0] [0 1 1] [0 0 0]};
% dMin = []; dMax = [];
tagname{1} = 'plot';
tagname{2} = 'pca';

% for i = 1:10
%     for j = 1:10
%         fh.plot{i,j} = uicontrol('Style','text','Tag',['plot',num2str((i-1)*10+j)],'Units','normalized',...
%         'HorizontalAlignment','Center',...
%         'FontWeight','Bold','FontUnits','normalized','FontSize',1.25,...
%         'BackgroundColor',[rand(1,3)],...
%         'String',num2str(round(rand*100)),'ForegroundColor',[1 1-rand 1-rand],...
%         'Position',[Xoffset+(j-1)*xoffset 1-Yoffset-(i-1)*yoffset xsize ysize]);
%         drawnow;
%     end
% end


%% Setup GUI default values

set(fh.expnev1.cluster.popup.n,'Value',2);
set(fh.expnev1.cluster.popup.nwaves,'Value',3);
% Run ClusterPopupN Callback once, to make sure d' box is appropriately
% enabled/disabled
[fh,nClusters,k,nwaves] = ClusterPopupN_Callback(...
   fh.expnev1.cluster.popup.n,[],fh,nevopen_outcome);
fname = [];
% AZ20090502 from: LB 080808 to allow automatic file specification
if length(varargin) == 3
   animaltype = varargin{1};
   iseries    = varargin{2};
   iexp       = varargin{3};
   fdir       = [DIRS.Cerebus filesep animaltype filesep];
   fname      = sprintf('u%03d_%03d.nev',iseries,iexp);
   set(fh.expnev1.fig,'Name',['explorenev8: ',animaltype,',',num2str(iseries),',',...
      num2str(iexp)]);
   
   [fh,nevopen_outcome,animaltype,iseries,iexp,nevsorted,mwaves,...
   SamplingRateInKHZ,xAxisInMsec,A,nelecs,X,d,class_idx,Xraw,iXraw,V,...
   threshold,scaling_params,calculated,k,animal_spikesdir,already_sorted_998] = ...
   load_button_press([],[],fh,nevopen_outcome,animaltype,iseries,iexp,threshold,...
   nwaves,calculated,nClusters,fname,fdir,tagname,plotColors);
end

% Add variables just created to base workspace
varnames = whos;
for ivar = 1:numel(varnames)
   eval(sprintf('assignin(''base'',''%s'',%s)',...
      varnames(ivar).name,varnames(ivar).name));
end

end
%% % End Main function (the following functions are NO LONGER NESTED)

% %% CREATE (UI Controls) FUNCTIONS

%% Load .NEV File Function
function [fh,nevopen_outcome,animaltype,iseries,iexp,nevsorted,mwaves,...
   SamplingRateInKHZ,xAxisInMsec,A,nelecs,X,d,class_idx,Xraw,iXraw,V,...
   threshold,scaling_params,calculated,k,animal_spikesdir,already_sorted_998] = ...
   load_button_press(hObject,eventdata,fh,nevopen_outcome,animaltype,iseries,iexp,threshold,...
   nwaves,calculated,nClusters,fname,fdir,tagname,plotColors)

   global DIRS
   
fromBaseWS = false;
if ~exist('fh','var')
   % Get updated vars from base workspace
   fromBaseWS = true;
   varsToLoad = {'fh','nevopen_outcome','calculated','nClusters','tagname',...
      'plotColors','A','nwaves'};
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end
end

   % IF RELOADING, CLEAR EXISTING DATA
   if exist('nevopen_outcome','var') && ~isempty(nevopen_outcome)
      calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated);
      fname = []; fdir = [];
   end

   [fh,nevopen_outcome,animaltype,fname,nevsorted,mwaves,SamplingRateInKHZ,...
      nelecs,nSamplesPerPacket] = opennevfile(fh,fname,fdir);
   
   xAxisInMsec = (1:nSamplesPerPacket)/SamplingRateInKHZ;

   a = regexp(fname,'[u_\.]');
   iseries = str2double(fname(a(1)+1:a(2)-1));
   iexp    = str2double(fname(a(2)+1:a(3)-1));
   %     i = strfind(fdir,filesep);
   %     animaltype = fdir(i(end-1)+1:i(end)-1);

   [A,C] = UtahGetLayout(animaltype,iseries); A
   elecs = A(~isnan(A))';
   nelecs = nnz(elecs);
   X         =  cell(nelecs,1);
   d         =  cell(nelecs,1);
   k         = zeros(nelecs,1);
   class_idx =  cell(nelecs,5);  %CLASS_IDX: set second dimension to be max nClusters possible
   elec      = [];
   Xraw      =  cell(nelecs,1);
   iXraw     =  cell(nelecs,1);
%    Xrealigned = cell(nelecs,1);
   V         =  cell(nelecs,1);
   S         =  cell(nelecs,1);
   U         =  cell(nelecs,1);
   C         =  cell(nelecs,1);
   threshold = str2double(get(fh.expnev1.edit.thresh,'String'))*ones(nelecs,1);

%   [fh,already_sorted,already_sorted_998] = ...
%   sortedUnitCheck(fh,animaltype,iseries)
    sortedUnitCheck;
   threshold = load998s(animaltype,iseries,animal_spikesdir,threshold,...
      already_sorted_998);
   
   [fh,X,Xraw,iXraw,scaling_params] = process_elecs(...
      fh,A,X,Xraw,iXraw,threshold,nwaves,nSamplesPerPacket);
%     set(fh.expnev1.show.panel,'Visible','on');
%% AZ 2010-05-05 Load units!
   gmm = [];
   prev         = struct([]);
   %AZ20090204: LOAD SAVED PARAMETERS
   for elec = elecs
   already_sorted_to_load = already_sorted(already_sorted(:,1)==iseries & ...
                                           already_sorted(:,2)==iexp    & ...
                                           already_sorted(:,3)==(1000+elec),:);
   sorted_series{elec}  = already_sorted_to_load(:,1);
%    currUnitNum    = already_sorted_to_load(:,4);
   prevSavedUnits{elec} = already_sorted_to_load(:,5);
      [X{elec},V{elec},xAxisInMsec,class_idx(elec,:),d{elec},threshold(elec),fh,...
           iXraw{elec}] = ...
       prevUnitProcess(X{elec},gmm,prev,xAxisInMsec,Xraw{elec},iXraw{elec},nevopen_outcome,...
                       SamplingRateInKHZ,already_sorted_to_load,animaltype,...
                       sorted_series,animal_spikesdir,prevSavedUnits,iseries,...
                       iexp,threshold(elec),fh,nSamplesPerPacket);
   end
    [class_idx,d,fh,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
     class_idx,A,iXraw,fh,calculated,V);

    set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Plotting Electrode data...');
    drawnow;
%     plotAllElecs;
    calculated = elecSubplotCreate(         1,tagname,A,  calculated);
   [fh,calculated] = showorhide(fh.expnev1.show.box.Num   ,[],...
      fh,[],nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,...
      Xraw,xAxisInMsec,plotColors,k,class_idx,threshold);
   [fh,calculated] = showorhide(fh.expnev1.show.box.Traces,[],...
      fh,[],nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,...
      Xraw,xAxisInMsec,plotColors,k,class_idx,threshold);
   [fh,calculated] = showorhide(fh.expnev1.show.box.D     ,[],...
      fh,[],nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,...
      Xraw,xAxisInMsec,plotColors,k,class_idx,threshold);
   [fh,calculated] = showorhide(fh.expnev1.show.box.Axes  ,[],...
      fh,[],nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,...
      Xraw,xAxisInMsec,plotColors,k,class_idx,threshold);
    drawnow;
    set(fh.expnev1.fig,'Name',sprintf('explorenev8: %s,%2d,%2d',animaltype,iseries,iexp));
    set(findobj(gcf,'Tag','expnev1.txt.status'),'String','');
   set(gcf,'Pointer','arrow');
   
   if fromBaseWS
      varsToSave = {'fh','nevopen_outcome','animaltype','iseries','iexp',...
         'nevsorted','mwaves','SamplingRateInKHZ','xAxisInMsec','A','nelecs',...
         'X','d','class_idx','Xraw','iXraw','V','threshold','scaling_params',...
         'calculated','k','animal_spikesdir','already_sorted_998'};
      for thevar = varsToSave
         eval(sprintf('assignin(''base'',''%s'',%s)',thevar{:},thevar{:}));
      end
   end
end
% END Load .NEV File function


%% Load electrode function
function Xraw = loadelec(elec,nwaves,nSamplesPerPacket)

    % AZ20090126: now chooses random sample of spikes, size nwaves
%     [timestamps_wf, waveforms] = nevwaves(elec,[1 nwaves]);
    [T, waveforms] = nevwaves(elec,nwaves,nSamplesPerPacket);  %% Get the first nwaves

   Xraw = waveforms;

   if isempty(Xraw)
      set(findobj(gcf,'Tag','expnev1.txt.status'),'String','No data!  Electrode skipped!');
      return;
   end

    if(isempty(Xraw))
        warning('explorenev:process_elecs:loadelec:NoData',...
           'No data!  Electrode #%d skipped!',elec);
        return;
    end

end
% END Load electrode function


%% Load 998 thresholds, if any
function threshold = load998s(animaltype,iseries,animal_spikesdir,threshold,...
   already_sorted_998)
global DIRS
   if exist('already_sorted_998','var') && ~isempty(already_sorted_998)
      for elec = already_sorted_998(:,3)'-1000
         unit = load(sprintf('%s%u%s',[DIRS.spikes filesep animaltype filesep],...
            iseries,[filesep ...
           animal_spikesdir{already_sorted_998(...
                            already_sorted_998(:,3)==elec+1000,5)}]));
         threshold(elec) = unit.unit.gmm.threshold;
         clear unit
      end
   end
end
% END load 998 thresholds


%% Process loaded electrode data
function [fh,X,Xraw,iXraw,scaling_params] = process_elecs(...
   fh,A,X,Xraw,iXraw,threshold,nwaves,nSamplesPerPacket)

    set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Loading Electrode data... 0%');

    yMax = zeros(nnz(~isnan(A)),1);
    yMin = zeros(nnz(~isnan(A)),1);
    
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            elec = A(i,j);
            if ~isnan(elec)
               Xraw{elec} = loadelec(elec,nwaves,nSamplesPerPacket);   %% load electrode data
                
               % Re-align offset traces (anything with a global minimum ocurring more than
               %   3 samples after the mode of the minima for all traces is moved to have
               %   its minimum at the mode)
               % [xr,yr] = temporary variables (values & indices of minima)
               [yr,yr] = min(Xraw{elec},[],1);
               for a = find(yr-mode(yr)>3)
                  Xraw{elec}(:,a) = [Xraw{elec}(yr(a)-mode(yr):end,a)' zeros(1,yr(a)-mode(yr)-1)];
               end

               [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},fh);
                
               if ~isempty(X{elec})
                  yMax(elec) = max(max(X{elec}(:,1:end)));
                  yMin(elec) = min(min(X{elec}(:,1:end)));
               else
                  yMax(elec) =  1;
                  yMin(elec) = -1;
               end
            end
            set(findobj(gcf,'Tag','expnev1.txt.status'),'String',['Loading Electrode data... ',num2str((i-1)*10+j),'%']);
            drawnow;
        end
    end

    yMaxALL = max(yMax);
    yMinALL = min(yMin);
    ryALL   = yMaxALL-yMinALL;
    
    scaling_params.fitChange = false;
    if get(fh.expnev1.scaling.button.global,'Value') && ~get(fh.expnev1.scaling.button.fit,'Value')
%         scaling_params.fitChange = true;
        scaling_params.yAxisMin = (yMinALL-0.1*ryALL)*ones(size(yMin));
        scaling_params.yAxisMax = (yMaxALL+0.1*ryALL)*ones(size(yMax));
    elseif ~get(fh.expnev1.scaling.button.global,'Value') && get(fh.expnev1.scaling.button.fit,'Value')
%         scaling_params.fitChange = true;
        scaling_params.yAxisMin = yMin-0.1*(yMax-yMin);
        scaling_params.yAxisMax = yMax+0.1*(yMax-yMin);
    else
        warning('Invalid Scaling Option');
    end
    
    scaling_params.yMaxALL = yMaxALL;
    scaling_params.yMinALL = yMinALL;
    scaling_params.ryALL   = ryALL;
    scaling_params.yMax    = yMax;
    scaling_params.yMin    = yMin;

end
% END process_elecs function


%% Plot electrode axes function
function calculated = plotElecAxes(nevopen_outcome,calculated)
    if ~isempty(nevopen_outcome)
%         if get(hObject,'Value') ~= strcmp(get(findobj('Tag',['elec1',tagname{fignum}]),'Visible'),'on')
%             % does this case ever happen?
%             for a = 1:96
%                 eval(sprintf(['hObject = findobj(''Tag'',[''elec'',num2str(a),''',tagname{fignum},''']);']));
%                 set(hObject,'Visible',on_or_off{b});
%             end
%         end
        
%         if fignum == 2 && calculated(fignum).PCA == 0
%             for a = 1:96
%                 set(findobj('Tag',['elec',num2str(a),tagname{fignum}]),...
%                     'YLim',[-0.5 0.5])
%             end
%         end
    end
    calculated(1).Axes = 0;  % make sure this is always 0
end
% END Plot electrode axes function


%% Showorhide function
function [fh,calculated,f2] = showorhide(hObject,eventdata,fh,f2,nevopen_outcome,...
   scaling_params,tagname,calculated,A,d,X,iXraw,Xraw,xAxisInMsec,plotColors,k,...
   class_idx,threshold)
% How SHOWORHIDE works: First, test if any data has been loaded already. If
% yes, do nothing (return), if no move on.  Second, prepare strings to be
% executed: 
%   (#1) test if a given feature (e.g., d') has already been calculated, and
%   determine value of its checkbox.  if not calculated, and checkbox is
%   checked, run plotElec* and changed calculated flag to 1.  if already
%   calculated, turn visibility on or off (#2)
%   Also, run code (#1) when re-scaling trace graphs

object_suffix = get(hObject,'Tag');
object_suffix = object_suffix(max(strfind(object_suffix,'.'))+1:end);

fromBaseWS = false;
if ~exist('nevopen_outcome','var')
   nevopen_outcome = evalin('base','nevopen_outcome');
   fromBaseWS = true;
end

% Test whether any data has been loaded or not
if isempty(nevopen_outcome)
   return;
end

if fromBaseWS
   % LOAD VARS FROM BASE WORKSPACE
   varsToLoad = {'fh','f2','tagname','calculated','scaling_params'};
   switch object_suffix
      case {'D'; 'Num'}
         varsToLoad = [varsToLoad,'A','d'];
      case 'Axes'
         % do nothing
      case {'Traces'; 'PCA'}
         varsToLoad = [varsToLoad,'A','X','iXraw','Xraw','xAxisInMsec','plotColors',...
            'k','class_idx','threshold'];
      otherwise
         % do nothing
   end
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end
end
   
switch gcf
   case fh.expnev1.fig
      fignum = 1;
   case f2
      fignum = 2;
   otherwise
      fignum = 1;
end

set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Refreshing view, please wait...');
drawnow;

on_or_off = {'on';'off'};
if get(hObject,'Value')
   b = 1; % on_or_off->'on'
else
   b = 2; % on_or_off->'off';
end
        
% PREPARE CODE: find all handles tagged 'elec(#)plot(suffix)'.
% traces' handles: all children of 'elec(#)plot', except elec(#)plotD and
% elec(#)plotNum
switch object_suffix
   case {'D'; 'Num'}
      cmd = sprintf(['hObject = findobj(''Tag'',[''elec'',num2str(a),''',tagname{fignum},''',object_suffix]);']);
   case 'Axes'
      cmd = sprintf(['hObject = findobj(''Tag'',[''elec'',num2str(a),''',tagname{fignum},''']);']);
   case {'Traces'; 'PCA'}
      cmd = sprintf(['hObject = allchild(findobj(''Tag'',[''elec'',num2str(a),''',tagname{fignum},''']));\n',...
         'h = findobj(''Tag'',[''elec'',num2str(a),''',tagname{fignum},''',''D'']);\n',...
         'if ~isempty(h); hObject(find(hObject==h(1))) = []; end\n',...
         'h = findobj(''Tag'',[''elec'',num2str(a),''',tagname{fignum},''',''Num'']);\n',...
         'if ~isempty(h); hObject(find(hObject==h(1))) = []; end\n']);
   otherwise
      set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Error: Bad ''Show'' selection.');
end
        
% RUN CODE
%(#1)
if (eval(sprintf(['~calculated(fignum).',object_suffix])) && b == 1) || ...
      (strcmp(object_suffix,'Traces') && scaling_params.fitChange) % (b=1 is 'on')
   eval(sprintf(['calculated(fignum).',object_suffix, ' = 1;']));
   switch object_suffix
      case {'D'; 'Num'}
         eval(sprintf(['plotElec',object_suffix,'(fignum,tagname,A,d);']));
      case 'Axes'
         eval(sprintf(['calculated = plotElec',object_suffix,'(nevopen_outcome,calculated);']));
      case {'Traces'; 'PCA'}
         eval(sprintf(['plotElec',object_suffix,'(fignum,tagname,A,X,',...
            'iXraw,Xraw,xAxisInMsec,plotColors,k,class_idx,threshold,',...
            'scaling_params.yAxisMin,scaling_params.yAxisMax);']));
      otherwise
         set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Error: Bad ''Show'' selection.');
   end
elseif eval(sprintf(['calculated(fignum).',object_suffix]))% && b == 2
   
   %(#2)
   for a = 1:96
      eval(cmd);
      set(hObject,'Visible',on_or_off{b});
   end
   
   if b == 2
      eval(sprintf(['calculated(fignum).',object_suffix, ' = 1;']));
   end
   
end

set(findobj(gcf,'Tag','expnev1.txt.status'),'String','');

if fromBaseWS
   varsToSave = {'fh','f2','calculated'};
   for thevar = varsToSave
      eval(sprintf('assignin(''base'',''%s'',%s)',thevar{:},thevar{:}));
   end
end
end
% END showorhide function


%% ClusterButton callback function
function [fh,d,class_idx,elec,calculated,k,V,scaling_params,X,Xraw,iXraw] = ...
   ClusterButton_Callback(hObject,eventdata,...
   X,nwaves,nClusters,A,iXraw,fh,calculated,V,Xraw,threshold,nevopen_outcome,...
   xAxisInMsec,plotColors)

   tagname = evalin('base','tagname');
   f2      = evalin('base','f2'     );
   
   fromBaseWS = false;
   if ~exist('fh','var')
      % Get updated vars from base workspace
      fromBaseWS = true;
      varsToLoad = {'fh','X','nwaves','nClusters','A','iXraw',...
            'calculated','V','Xraw','threshold','nevopen_outcome',...
            'xAxisInMsec','plotColors'};
      for thevar = varsToLoad
         eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
      end
   end

%     if ~isempty(nevopen_outcome)
   nelecs    = nnz( ~isnan(A));
   d         = zeros(nelecs,1);
   class_idx =  cell(nelecs,4); %CLASS_IDX: set second dimension to be max nClusters possible
   elec      = [];

   % If figure 2 is open, notify user it will be redrawn
   if exist('f2','var') && ~isempty(f2)
      set(eval(sprintf(['fh.expnev',num2str(2),'.txt.status'])),'String','Refreshing view, please wait...');
      drawnow;
   end

   if size(X{1},2) ~= nwaves
      [fh,X,Xraw,iXraw,scaling_params] = process_elecs(...
      fh,A,X,Xraw,iXraw,threshold,nwaves);
   end

   [class_idx,d,fh,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
   class_idx,A,iXraw,fh,calculated,V);

   set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Plotting Electrode data...');

   calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated);
   drawnow;
   %         plotAllElecs;
   [fh,calculated,f2] = showorhide(fh.expnev1.show.box.Num   ,[],fh,f2,...
      nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,Xraw,...
      xAxisInMsec,plotColors,k,class_idx,threshold);
   [fh,calculated,f2] = showorhide(fh.expnev1.show.box.Traces,[],fh,f2,...
      nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,Xraw,...
      xAxisInMsec,plotColors,k,class_idx,threshold);
   [fh,calculated,f2] = showorhide(fh.expnev1.show.box.D     ,[],fh,f2,...
      nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,Xraw,...
      xAxisInMsec,plotColors,k,class_idx,threshold);
   [fh,calculated,f2] = showorhide(fh.expnev1.show.box.Axes  ,[],fh,f2,...
      nevopen_outcome,scaling_params,tagname,calculated,A,d,X,iXraw,Xraw,...
      xAxisInMsec,plotColors,k,class_idx,threshold);
   drawnow;
   set(findobj(gcf,'Tag','expnev1.txt.status'),'String','');

   % If figure 2 is open, redraw it
   if ~isempty(f2)
      pcafigure(hObject,eventdata);
   end
%     end

   % Update vars in base workspace
   if fromBaseWS
      varsToSave = {'fh','d','class_idx','elec','calculated','k','V',...
         'scaling_params','X','Xraw','iXraw'};
      for thevar = varsToSave
         eval(sprintf('assignin(''base'',''%s'',%s)',thevar{:},thevar{:}));
      end
   end
end
% END ClusterButton callback function


%% ClusterPopupN callback function
function [fh,nClusters,k,nwaves,d,class_idx,elec] = ClusterPopupN_Callback(...
   hObject,eventdata,fh,nevopen_outcome,X,nwaves,nClusters,A,iXraw,calculated,...
   V,Xraw,threshold,xAxisInMsec,plotColors)

fromBaseWS = false;
if ~exist('fh','var')
   % Get updated vars from base workspace
   fromBaseWS = true;
   varsToLoad = {'fh'};
   nevopen_outcome = evalin('base','nevopen_outcome');
   if ~isempty(nevopen_outcome)
      varsToLoad = [varsToLoad,'X','nwaves','nClusters','A','iXraw',...
         'calculated','V','Xraw','threshold','xAxisInMsec','plotColors'];
   end
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end
end

nClusters = get(fh.expnev1.cluster.popup.n,'String');
nClusters = nClusters{get(fh.expnev1.cluster.popup.n,'Value')};
nClusters = str2double(nClusters);
try k = nClusters*ones(96,1); catch, k = zeros(96,1); end

nwaves = get(fh.expnev1.cluster.popup.nwaves,'String');
nwaves = str2double(nwaves{get(fh.expnev1.cluster.popup.nwaves,'Value')});% # of waves
if nClusters == 1
   set(fh.expnev1.show.box.D,'Value',0);
   set(fh.expnev1.show.box.D,'Enable','off');
else
   %         set(fh.expnev1.show.box.D,'Value',1);
   set(fh.expnev1.show.box.D,'Enable','on');
end

if ~isempty(nevopen_outcome)
   [fh,d,class_idx,elec,calculated,k,V,scaling_params,X,Xraw,iXraw] = ...
      ClusterButton_Callback(hObject,eventdata,...
      X,nwaves,nClusters,A,iXraw,fh,calculated,V,Xraw,threshold,...
      nevopen_outcome,xAxisInMsec,plotColors);
end

% Update vars in base workspace
if fromBaseWS
   varsToSave = {'fh', 'nevopen_outcome'};
   if ~isempty(nevopen_outcome)
      varsToSave = [varsToSave,'d','class_idx','elec'];
   end
   for thevar = varsToSave
      eval(sprintf('assignin(''base'',''%s'',%s)',thevar{:},thevar{:}));
   end
end
end
% END ClusterPopupN callback function


%% editThreshold callback function
function editThreshold_Callback(hObject, eventdata)
   % LOAD VARS FROM BASE WORKSPACE
   varsToLoad = {'fh','A','X','iXraw','Xraw','scaling_params','d','nClusters',...
      'class_idx','calculated','V','tagname','f2'};
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end
   
   threshold = str2double(get(fh.expnev1.edit.thresh,'String'))*ones(96,1);
   
   for i = 1:size(A,1)
      for j = 1:size(A,2)
         elec = A(i,j);
         if ~isnan(elec)
            [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},fh);
            
            if ~isempty(X{elec})
               scaling_params.yMax(elec) = max(max(X{elec}(:,1:end)));
               scaling_params.yMin(elec) = min(min(X{elec}(:,1:end)));
            else
               scaling_params.yMax(elec) =  1;
               scaling_params.yMin(elec) = -1;
            end
         end
         set(findobj(gcf,'Tag','expnev1.txt.status'),'String',...
            sprintf('Rethresholding data... %i%%',...
            100*sub2ind(size(A'),j,i)/numel(A)));
         drawnow;
      end
   end
   
   [class_idx,d,fh,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
     class_idx,A,iXraw,fh,calculated,V);
   
   calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated);
   drawnow;
   
   % SAVE VARS TO BASE WORKSPACE
   varsToSave = {'fh','X','iXraw','scaling_params','d','class_idx','k','V',...
      'calculated'};
   for thevar = varsToSave
      eval(sprintf('assignin(''base'',''%s'',%s)',thevar{:},thevar{:}));
   end
   
   %         plotAllElecs;
   showorhide(fh.expnev1.show.box.Num);
   showorhide(fh.expnev1.show.box.Traces);
   showorhide(fh.expnev1.show.box.D);
   showorhide(fh.expnev1.show.box.Axes);
   drawnow;
   set(findobj(gcf,'Tag','expnev1.txt.status'),'String','');

   % If figure 2 is open, redraw it
   if ~isempty(f2)
      pcafigure(hObject,eventdata);
   end
   set(gcf,'Pointer','arrow');
end
% END editThreshold callback function


%% save998 callback function
function save998_Callback(hObject, eventdata)
   % LOAD VARS FROM BASE WORKSPACE
   varsToLoad = {'A','animaltype','iseries','iexp','threshold'};
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end
   
   set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Saving 998s...');
   drawnow;

%    threshold = str2double(get(fh.expnev1.edit.thresh,'String'))*ones(96,1);
   ichans = unique(A(~isnan(A)))';
   create998(animaltype,iseries,iexp,ichans,threshold);
   
   set(findobj(gcf,'Tag','expnev1.txt.status'),'String','');
end
% END save998 callback function


%% Scaling Buttons function
function scaling_button_press(ignoreme, eventdata)
   % LOAD VARS FROM BASE WORKSPACE
   varsToLoad = {'scaling_params','nevopen_outcome','fh'};
   for thevar = varsToLoad
      eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
   end
    % sp = scaling_params as defined in the process_elecs function
    sp = scaling_params;
    if ~isempty(nevopen_outcome)
        if strcmp(get(eventdata.NewValue,'String'),'Fit')
            sp.fitChange = true;
            sp.yAxisMin = sp.yMin-0.1*(sp.yMax-sp.yMin);
            sp.yAxisMax = sp.yMax+0.1*(sp.yMax-sp.yMin);
            % TODO: hide y-axis ticks
        elseif strcmp(get(eventdata.NewValue,'String'),'Global')
            sp.fitChange = true;
            sp.yAxisMin = (sp.yMinALL-0.1*sp.ryALL)*ones(size(sp.yMin));
            sp.yAxisMax = (sp.yMaxALL+0.1*sp.ryALL)*ones(size(sp.yMax));
            % TODO: show y-axis ticks
        else
            warning('Invalid Scaling Option');
        end

        assignin('base','scaling_params',sp);

        set(findobj(gcf,'Tag','expnev1.txt.status'),'String','Plotting Electrode data...');
        drawnow;
%         plotAllElecs;
%         showorhide(fh.expnev1.show.box.Num);

        showorhide(fh.expnev1.show.box.Traces);
        
        showorhide(fh.expnev1.show.box.D);
        showorhide(fh.expnev1.show.box.Axes);
        drawnow;
        set(findobj(gcf,'Tag','expnev1.txt.status'),'String','');
        
        sp.fitChange = false;
        assignin('base','scaling_params',sp);
    end
end
% END Global Scaling Button function


%% Second Figure (PCA graphs) function
function pcafigure(hObject,eventdata)
    if ~isempty(f2) && strcmp(get(hObject,'Tag'),'expnev1.PCABox')
        calculated = elecSubplotDestroyChildren(2,tagname,A,nClusters,calculated);
        cancel_f2;
    else %if isempty(f2) && strcmp(get(hObject,'Tag'),'expnev1.PCABox') || strcmp(get(hObject,'Tag'),'ClusterButton')
        % Initialize Figure
        if strcmp(get(hObject,'Tag'),'expnev1.PCABox')
            f2 = figure('Position',[520 380 672 420],'Name','explorenev6-PCA',...
                'Toolbar','figure','MenuBar','none','NumberTitle','off');
            % Create the 'status' text area
            f2h.status = uicontrol('Style','text','Tag','status','Units','normalized',...
                'Position',[0.8 0 0.2 0.075],'BackgroundColor',get(f2,'Color'));
            f2h.show.panel = uipanel('Title','Show','Tag','show.panel',...%'Visible','off',...
                'Units','normalized','Position',[0.8 0.225 0.2 0.15],...
                'BackgroundColor',get(f2,'Color'));
            f2h.show.box.Num = uicontrol('Parent',f2h.show.panel,'String','Elec #',...
                'Style','checkbox','Tag','show.box.Num','Value',0,...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[11 4 60 23],...
                'callback',@showorhide);
            f2h.show.box.Axes = uicontrol('Parent',f2h.show.panel,'String','Axes',...
                'Style','checkbox','Tag','show.box.Axes','Value',1,...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[71 24 60 23],...
                'callback',@showorhide);
            f2h.show.box.PCA = uicontrol('Parent',f2h.show.panel,'String','Clusters',...
                'Style','checkbox','Tag','show.box.PCA','Value',1,...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[71 4 60 23],...
                'callback',@showorhide);
            f2h.show.box.D = uicontrol('Parent',f2h.show.panel,'String','d''',...
                'Style','checkbox','Tag','show.box.D','Value',1,'Enable','on',...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[11 24 40 23],...
                'callback',@showorhide);

            hObject = fh.expnev1.PCABox; figure(f2);
            calculated = elecSubplotCreate(         2,tagname,A,  calculated);
        else
            hObject = fh.expnev1.PCABox; figure(f2);
            calculated = elecSubplotDestroyChildren(2,tagname,A,nClusters,calculated);
        end
        
        if ~isempty(nevopen_outcome)
            showorhide(f2h.show.box.Num);
            plotElecPCA(2);
            showorhide(f2h.show.box.D);
            showorhide(f2h.show.box.Axes);
        else
            set(f2h.status,'String','No .nev file selected.');
        end
        drawnow;
    end
end
% END Second Figure (PCA graphs) function


%% Plot electrode traces function
function plotElecPCA(fignum,hObject,eventdata)%,hObject,eventdata,handles)
    %% TODO / BUG: (can't be fixed?  tried using uistack.  it seems that
    %% regardless of which electrode is plotted, after the 20th plot, the
    %% whole figure refreshes, & won't let you move the text/number objects
    %% behind the plotted points)
    for i = 1:10
        for j = 1:10
            elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
            if ~isnan(elec)
                h = findobj('Tag',['elec',num2str(elec),tagname{fignum}]);
%                 set(gcf,'CurrentAxes',h(1));
                subplot('Position',get(h,'Position'))
                
                for a = 1:k(elec)
                    if ~isempty(class_idx{elec,a})
                        scatter(V{elec}(class_idx{elec,a},1),V{elec}(class_idx{elec,a},2),...
                            5,plotColors{a});
                    end
                    xr = range(V{elec}(:,1)); yr = range(V{elec}(:,2));
                    axis([min(V{elec}(:,1))-0.10*xr max(V{elec}(:,1))+0.15*xr ...
                          min(V{elec}(:,2))-0.15*yr max(V{elec}(:,2))+0.10*yr  ])
                    hold on, plot(0,0,'w.','markersize',15); plot(0,0,'r.','markersize',7);
                end

            end
        end
    end
end
% END Plot electrode traces function


% %% Load vars from base workspace function
% function loadFromBaseWS(varsToLoad)
% for thevar = varsToLoad
%    eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
% end
% end
% % End Load vars from base workspace function
% 
% 
% %% Save vars to base workspace function
% function saveToBaseWS(varsToSave)
% for thevar = varsToSave
%    eval(sprintf('%s = evalin(''base'',''%s'');',thevar{:},thevar{:}));
% end
% end
% % End Save vars to base workspace function