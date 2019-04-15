function explorenev6(varargin)

%% HINT: If you can't close a window because of a CloseRequestFcn error, run:
%%       delete(get(0,'Children'))

%% TODO: Fix Figure 2 so that Elec #s are below clusters
global X Xraw iXraw threshold d nClusters class_idx A calculated k V ...
   xAxisInMsec plotColors yAxisMin yAxisMax

%% Begin Main Figure Main Function
% toggle: 'Menu'/'Toolbar' = 'none'/'figure'
if ~isempty(findobj(get(0,'Children'),'Tag','fig.explorenev'))
   f1h.fig.expnev1 = findobj(get(0,'Children'),'Tag','fig.explorenev');
   set(0,'CurrentFigure',f1h.fig.expnev1);clf;
else
   f1h.fig.expnev1 = figure('Position',[520 380 672 420],'Name','explorenev6',...
      'Toolbar','none','MenuBar','none','NumberTitle','off',...
      'Tag','fig.explorenev');
end
set(gcf,'units','normalized','outerposition',[0 0.05 1 0.95]);
% Create a variable to store the handle of the new GUI figure
% f2 becomes non-empty when the new figure is created.
f2 = []; f2h = [];

%% Initialize GUI objects
% Create the 'load' pushbutton
f1h.loadbutton = uicontrol(f1h.fig.expnev1,'Style', 'pushbutton','String','Load',...
    'Units','normalized','Tag','loadbutton','Position',[0.82 0.93 0.06 0.05],...
    'callback',@load_button_press);
% Create the 'status' text area
f1h.txt.status = uicontrol('Style','text','Tag','txt.status','Units','normalized',...
    'Position',[0.8 0 0.2 0.075],'BackgroundColor',get(f1h.fig.expnev1,'Color'));
% Create the 'Scaling' panel, containing scaling buttons
f1h.scalingpanel = uibuttongroup('Title','Scaling','Tag','scalingpanel',...
    'Units','normalized','Position',[0.8 0.075 0.2 0.15],...
    'BackgroundColor',get(f1h.fig.expnev1,'Color'),'SelectionChangeFcn',@scaling_button_press);
% Create the 'Global' scaling button
f1h.scalingButtonGlobal = uicontrol('Parent',f1h.scalingpanel,'String','Global',...
    'Style','radiobutton','Tag','scalingButtonGlobal','Value',0,...
    'Units','pixels','BackgroundColor',get(f1h.fig.expnev1,'Color'),'Position',[11 24 87 23]);
% Create the 'Fit'    scaling button
f1h.scalingButtonFit    = uicontrol('Parent',f1h.scalingpanel,'String','Fit',...
    'Style','radiobutton','Tag','scalingButtonFit'   ,'Value',1,...
    'Units','pixels','BackgroundColor',get(f1h.fig.expnev1,'Color'),'Position',[11 4 87 23]);
% Create the 'Show' panel, containing show boxes
f1h.showpanel = uipanel('Title','Show','Tag','showpanel',...%'Visible','off',...
    'Units','normalized','Position',[0.8 0.225 0.2 0.15],...
    'BackgroundColor',get(f1h.fig.expnev1,'Color'));
f1h.showBoxD = uicontrol('Parent',f1h.showpanel,'String','d''',...
    'Style','checkbox','Tag','showBoxD','Value',1,'Enable','off',...
    'Units','pixels','BackgroundColor',get(f1h.fig.expnev1,'Color'),'Position',[11 24 40 23],...
    'callback',@showorhide);  %@showBoxD_Callback);
f1h.showBoxNum = uicontrol('Parent',f1h.showpanel,'String','Elec #',...
    'Style','checkbox','Tag','showBoxNum','Value',1,...
    'Units','pixels','BackgroundColor',get(f1h.fig.expnev1,'Color'),'Position',[11 4 60 23],...
    'callback',@showorhide);
f1h.showBoxAxes = uicontrol('Parent',f1h.showpanel,'String','Axes',...
    'Style','checkbox','Tag','showBoxAxes','Value',1,...
    'Units','pixels','BackgroundColor',get(f1h.fig.expnev1,'Color'),'Position',[71 24 60 23],...
    'callback',@showorhide);
f1h.showBoxTraces = uicontrol('Parent',f1h.showpanel,'String','Traces',...
    'Style','checkbox','Tag','showBoxTraces','Value',1,...
    'Units','pixels','BackgroundColor',get(f1h.fig.expnev1,'Color'),'Position',[71 4 60 23],...
    'callback',@showorhide);
% Create the 'Clustering' panel
f1h.ClusteringPanel = uipanel('Title','Clustering','Tag','ClusteringPanel',...
    'Units','normalized','Position',[0.8 0.375 0.2 0.15],...%'Visible','off',...
    'BackgroundColor',get(f1h.fig.expnev1,'Color'));
% f1h.textClusters = uicontrol('Parent',f1h.ClusteringPanel,'Style','text',...
%     'Tag','textClusters',...
%     'Units','normalized','String','Clusters','Position',[0.43 0.35 0.5 0.5],...
%     'HorizontalAlignment','left','BackgroundColor',get(f1h.fig.expnev1,'Color'));
f1h.ClusteringNButton = uicontrol('Parent',f1h.ClusteringPanel,'Style','pushbutton',...
    'Tag','ClusteringNButton',...
    'Units','normalized','String','Cluster!','Position',[0.525 0.525 0.4 0.45],...
    'HorizontalAlignment','left','callback',@ClusteringNButton_Callback);
f1h.ClusteringNPopup = uicontrol('Parent',f1h.ClusteringPanel,'Style','popupmenu',...
    'Tag','ClusteringNPopup',...%'BackgroundColor',get(f1h.fig.expnev1,'Color'),...
    'String',{'1';'2';'3';'4';'best'},'Units','normalized',...
    'Position',[0.08 0.525 0.375 0.45],'callback',@ClusteringNPopup_Callback);
f1h.ClusteringNWavesPopup = uicontrol('Parent',f1h.ClusteringPanel,'Style',...
    'popupmenu','Tag','ClusteringNWavesPopup',...%'BackgroundColor',get(f1h.fig.expnev1,'Color'),...
    'String',{'100';'500';'1000';'2500'},'Units','normalized',...
    'Position',[0.08 0.10 0.375 0.45],'callback',@ClusteringNPopup_Callback);
f1h.ClusteringNWavesText = uicontrol('Parent',f1h.ClusteringPanel,'Style',...
    'text','Tag','ClusteringNWavesText','Units','normalized',...
    'Position',[0.525 0 0.4 0.45],'BackgroundColor',get(f1h.fig.expnev1,'Color'),...
    'String','# Waves');

f1h.PCABox = uicontrol('String','PCA (New Window)',...
    'Style','checkbox','Tag','PCABox','Value',0,'Enable','on',...
    'Units','normalized','BackgroundColor',get(f1h.fig.expnev1,'Color'),...
    'Position',[0.81 0.525 0.2 0.05],...
    'callback',@pcafigure);  %@showBoxD_Callback);
% f1h.SortnevLinksBox = uicontrol('String','Sortnev links',...
%     'Style','checkbox','Tag','PCABox','Value',1,'Enable','on',...
%     'Units','normalized','BackgroundColor',get(f1h.fig.expnev1,'Color'),...
%     'Position',[0.81 0.575 0.2 0.05],...
%     'callback',@pcafigure);  %@showBoxD_Callback);
 
% AZ20090220: Thresholding tool
% Create the 'Threshold' panel
% f1h.ThreshPanel = uipanel('Title','Threshold','Tag','ThreshPanel',...%'Visible','off',...
%     'Units','normalized','Position',[0.8 0.625 0.2 0.15],...
%     'BackgroundColor',get(f1h.fig.expnev1,'Color'));
f1h.txtThresh = uicontrol('Parent',f1h.fig.expnev1,'Units','normalized',...
   'HorizontalAlignment','left','BackgroundColor',get(f1h.fig.expnev1,'Color'),...
   'ListboxTop',0,'String','Threshold: +/-','Style','text','Tag','txtThresh',...
   'Position',[0.810 0.700 0.140 0.026]);
f1h.edThresh = uicontrol('Parent',f1h.fig.expnev1,'Units','normalized',...
   'BackgroundColor',[1 1 1],'FontSize',9,...
   'Callback',@edThresh_Callback,...%'KeyPressFcn',@shortcutKeys,...
   'ListboxTop',0,'String','32','Style','edit','Tag','edThresh',...
   'Position',[0.920 0.696 0.040 0.032]);
f1h.pshbtn.save998 = uicontrol('Parent',f1h.fig.expnev1,'Style','pushbutton',...
   'String','Save 998','Units','normalized','Tag','pshbtn.save998',...
   'Callback',@save998_Callback,...
   'Position',[0.890 0.930 0.085 0.050]);

% FOR subplots created in elecSubplotCreate
% Define a context menu; it is not attached to anything
f1h.menu.subplotfcns(1) = uicontextmenu('Tag','menu.subplotfcns(1)');
% Define the context menu items and install their callbacks
f1h.menu.subplotfcns(2) = uimenu(f1h.menu.subplotfcns(1),  ...
   'Label','Open in sortnev2','Tag','menu.subplotfcns(2)', ...
   'Callback', ...
   ['spname = get(get(gco,''Parent''),''Tag'');',...
    'elec = str2double(spname(5:end-4));',...
    'threshold(elec) = max(get(findobj(''Tag'',sprintf(''elec%dplotThresh'',elec)),''YData''));',...
    'fname = regexp(get(gcf,''Name''),''[\s,]'',''split'');',...
    'fname = fname(~strcmp(fname,''''));',...
    'sortnev2(fname{2},str2double(fname{end-1}),str2double(fname{end}),elec,threshold(elec));']);
%% TODO: figure out a way to avoid this quoted code absurdity
f1h.menu.subplotfcns(3) = uimenu(f1h.menu.subplotfcns(1),...
   'Label','Change Threshold','Tag','menu.subplotfcns(3)',...
   'Callback',{@editThresholdDialog, f1h});

%    ['spname=get(get(gco,''Parent''),''Tag'');',...
%     'elec = str2double(spname(5:6));',...
%     ' answer = inputdlg(sprintf(''New threshold: for elec %2.0f'',elec),sprintf(''Change Threshold (elec %f)'',elec),1,{''32''});',...
%     ' if ~isempty(answer{1});',...
%     '    threshold(elec) = str2double(answer{1});',...
%     '    [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},f1h);',...
%     ' end']);
% f1h.menu.subplotfcns(4) = uimenu(f1h.menu.subplotfcns(1),...
%    ...%'Callback', @clickEllipse_Callback,...
%    'Label','fprintf(''Electrode #%f'',elec)','Tag','menu.subplotfcns(4)');


% f1h.sidebar = uipanel('Title','Controls','Tag','ClusteringPanel',...%'Visible','off',...
%     'Units','normalized','Position',[0.8 0.375 0.2 0.095],...
%     'BackgroundColor',get(f1h.fig.expnev1,'Color'));


%% Initialize variables for common use across multiple nested functions
global DIRS serverName;
DIRS = struct('data',[],'spikes',[],'camera',[],'xfiles',[],'michigan',[],...
    'Cerebus',[],'stimInfo',[]);
serverName    = [];
serverDataDir = [];
SetDefaultDirs;
A = zeros(10,10);
a = []; b = []; h = []; m = []; n =[];
dColor    = []; d = [];
nevopen_outcome   = [];
yAxisMax  = zeros(96,1); yAxisMin  = zeros(96,1);
scaling_params    = [];  scaling_params.fit  = [];
on_or_off = {'on';'off'};

already_sorted    = []; animal_spikesdir  = []; already_sorted_to_load = [];
animaltype        = []; seriesdirs        = []; already_sorted_998     = [];
iseries           = []; iexp              = [];
nextUnitNum       = []; takenUnitNums     = [];
allseries         = []; parts             = [];

nClusters = get(f1h.ClusteringNPopup,'String');
nClusters = nClusters{get(f1h.ClusteringNPopup,'Value')};          % nClusters = # of classes/neurons
nClusters = str2double(nClusters);
try k = nClusters*ones(96,1); catch, k = zeros(96,1); end

nwaves = get(f1h.ClusteringNWavesPopup,'String');
nwaves = str2double(nwaves{get(f1h.ClusteringNWavesPopup,'Value')});% # of waves
SamplingRateInKHZ = 30;
xAxisInMsec = (1:48)/SamplingRateInKHZ;

elec      = [];
d         = zeros(96,1);
X         =  cell(96,1);
class_idx =  cell(96,4);  %CLASS_IDX: set second dimension to be max nClusters possible
V         =  cell(96,1);

Xraw      =  cell(96,1);
iXraw     =  cell(96,1);
Xrealigned = cell(96,1);
threshold = str2double(get(f1h.edThresh,'String'))*ones(96,1);

class_idx_temp    = cell(1,4);
d_temp    = []; k_temp    = [];
deltat    = [];

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
dMin = []; dMax = []; yMaxALL = []; yMinALL = []; ryALL = [];
tagname{1} = 'plot';
tagname{2} = 'pca';

% for i = 1:10
%     for j = 1:10
%         f1h.plot{i,j} = uicontrol('Style','text','Tag',['plot',num2str((i-1)*10+j)],'Units','normalized',...
%         'HorizontalAlignment','Center',...
%         'FontWeight','Bold','FontUnits','normalized','FontSize',1.25,...
%         'BackgroundColor',[rand(1,3)],...
%         'String',num2str(round(rand*100)),'ForegroundColor',[1 1-rand 1-rand],...
%         'Position',[Xoffset+(j-1)*xoffset 1-Yoffset-(i-1)*yoffset xsize ysize]);
%         drawnow;
%     end
% end


%% Setup GUI default values

set(f1h.ClusteringNPopup,'Value',2);
set(f1h.ClusteringNWavesPopup,'Value',3);
% Run ClusteringNPopup Callback once, to make sure d' box is appropirately
% enabled/disabled
ClusteringNPopup_Callback(f1h.ClusteringNPopup,[]);

fname = [];
% AZ20090502 from: LB 080808 to allow automatic file specification
if length(varargin) == 3
   animaltype = varargin{1};
   iseries    = varargin{2};
   iexp       = varargin{3};
   fdir       = [DIRS.Cerebus filesep animaltype filesep];
   fname      = sprintf('u%03d_%03d.nev',iseries,iexp);
   set(f1h.fig.expnev1,'Name',['explorenev6: ',animaltype,',',num2str(iseries),',',...
      num2str(iexp)]);
   
   load_button_press;
end

%% % End Main function (all that follows is Nested functions)


% %% CREATE (UI Controls) FUNCTIONS

%% Load .NEV File Function
function load_button_press(hObject, eventdata)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% IF RELOADING, CLEAR EXISTING DATA
if nevopen_outcome
   d         = zeros(96,1);
   X         =  cell(96,1);
   class_idx =  cell(96,4);  %CLASS_IDX: set second dimension to be max nClusters possible
   elec      = [];
   calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated);
end
    
[f1h,nevopen_outcome,animaltype,fname] = opennevfile(f1h,fname,fdir);
    

   a = regexp(fname,'[u_\.]');
    iseries = str2double(fname(a(1)+1:a(2)-1));
    iexp    = str2double(fname(a(2)+1:a(3)-1));
    i = strfind(fdir,'\');
    animaltype = fdir(i(end-1)+1:i(end)-1);

    [A,C] = UtahGetLayout(animaltype,iseries); A

    sortedUnitCheck;
    load998s;
   
    process_elecs;
%     set(f1h.showpanel,'Visible','on');
    [class_idx,d,f1h,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
     class_idx,A,iXraw,f1h,calculated,V);

    set(findobj(0,'Tag','txt.status'),'String','Plotting Electrode data...');
    drawnow;
%     plotAllElecs;
    calculated = elecSubplotCreate(         1,tagname,A,  calculated);
    showorhide(f1h.showBoxNum);
    showorhide(f1h.showBoxTraces);
    showorhide(f1h.showBoxD);
    showorhide(f1h.showBoxAxes);
    drawnow;
    set(f1h.fig.expnev1,'Name',sprintf('explorenev6: %s,%2d,%2d',animaltype,iseries,iexp));
    set(findobj(gcf,'Tag','txt.status'),'String','');
   set(gcf,'Pointer','arrow');
end
% END Load .NEV File function


%% Load electrode function
function Xraw = loadelec(elec,nwaves)

    % AZ20090126: now chooses random sample of spikes, size nwaves
%     [timestamps_wf, waveforms] = nevwaves(elec,[1 nwaves]);
    [T, waveforms] = nevwaves(elec,nwaves);  %% Get the first nwaves

   Xraw = waveforms;

   if isempty(Xraw)
      set(findobj(gcf,'Tag','txt.status'),'String','No data!  Electrode skipped!');
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
function load998s
   if exist('already_sorted_998','var') && ~isempty(already_sorted_998)
      for elec = already_sorted_998(:,3)'-1000
         unit = load([DIRS.spikes filesep animaltype filesep num2str(iseries) filesep ...
            animal_spikesdir{already_sorted_998(already_sorted_998(:,3)==elec+1000,5)}]);
         threshold(elec) = unit.unit.gmm.threshold;
         clear unit
      end
   end
end
% END load 998 thresholds


%% Process loaded electrode data
function process_elecs

    set(findobj(0,'Tag','txt.status'),'String','Loading Electrode data... 0%');

    for i = 1:size(A,1)
        for j = 1:size(A,2)
            elec = A(i,j);
            if ~isnan(elec)
               Xraw{elec} = loadelec(elec,nwaves);   %% load electrode data
                
               % Re-align offset traces (anything with a global minimum ocurring more than
               %   3 samples after the mode of the minima for all traces is moved to have
               %   its minimum at the mode)
               % [xr,yr] = temporary variables (values & indices of minima)
               [xr,yr] = min(Xraw{elec},[],1);
               for a = find(yr-mode(yr)>3)
                  Xraw{elec}(:,a) = [Xraw{elec}(yr(a)-mode(yr):end,a)' zeros(1,yr(a)-mode(yr)-1)];
               end

               [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},f1h);
                
               if ~isempty(X{elec})
                  yMax(elec) = max(max(X{elec}(:,1:end)));
                  yMin(elec) = min(min(X{elec}(:,1:end)));
               else
                  yMax(elec) =  1;
                  yMin(elec) = -1;
               end
            end
            set(findobj(0,'Tag','txt.status'),'String',['Loading Electrode data... ',num2str((i-1)*10+j),'%']);
            drawnow;
        end
    end

    yMaxALL = max(yMax);
    yMinALL = min(yMin);
    ryALL   = yMaxALL-yMinALL;

    scaling_params.yMaxALL = yMaxALL;
    scaling_params.yMinALL = yMinALL;
    scaling_params.ryALL   = ryALL;
    scaling_params.yMax    = yMax;
    scaling_params.yMin    = yMin;
    
    if get(f1h.scalingButtonGlobal,'Value') && ~get(f1h.scalingButtonFit,'Value')
        yAxisMin = (yMinALL-0.1*ryALL)*ones(size(yMin));
        yAxisMax = (yMaxALL+0.1*ryALL)*ones(size(yMax));
    elseif ~get(f1h.scalingButtonGlobal,'Value') && get(f1h.scalingButtonFit,'Value')
        yAxisMin = yMin-0.1*(yMax-yMin);
        yAxisMax = yMax+0.1*(yMax-yMin);
    else
        warning('Invalid Scaling Option');
    end

end
% END process_elecs function


%% Classify (using Gaussian Mixture Modeling) function
function [class_idx_temp,d_temp,best_i] = classify_gmm(d_temp,elec,Xelec,k_temp,class_idx_temp)
% nClusters = # of classes

    % CLASSIFY: GMM, then PCA->d'
    if k_temp > 1
        
        [U,S,V{elec}] = svds(Xelec,3); % SVD
        
        obj = [];
        obj = gmdistribution.fit(V{elec},k_temp,'Start','randSample');
%         dim = 1:size(V,2);

%         for j = dim
%             n = dim(find(dim~=j));
%             fit(j).gmm = gmdistribution(obj.mu(:,n),obj.Sigma(n,n,:),obj.PComponents);
% 
% %             figure(j); hold on;
% %             scatter(V(:,n(1)),V(:,n(2)),10,'.')
% %             conth(j) = ezcontour(@(x,y)pdf(fit(j).gmm,[x y]),xlim,ylim);
% 
%             Cn(j,:) = cluster(fit(j).gmm,V(:,n));
%         end

        C = cluster(obj,V{elec});

        for a = 1:k_temp
            class_idx_temp{a} = find(C==a);
        end
        
        % Calculate d' (for all permutations)

        [d_temp,best_i] = dprime(V{elec},class_idx_temp);

    elseif k_temp == 1 % i.e., if NOT clustering
        class_idx_temp{1} = (1:length(Xelec))';
        % make mean(d) = 0.5
        d_temp = 0.5;
    end

end
% END Classify (using Gaussian Mixture Modeling) function


% %% Plot all electrodes function
% function plotAllElecs
%     dMax    = max(d);
%     dMin    = min(d);
%     
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
%             if ~isnan(elec)
%                 subplot('Position',...
%                     [Xoffset+(j-1)*xoffset 1-Yoffset-(i-1)*yoffset xsize ysize],...
%                     'Tag',['elec',num2str(elec),'plot']);
%                 hold on;
%                 
%                 dColor = [1 1-(d(elec)-dMin)/(dMax-dMin) 1-(d(elec)-dMin)/(dMax-dMin)];
%                 text('Tag',['elec',num2str(elec),'plotNum'],'Position',[0.5 0.5],...
%                     'String',num2str(elec),'Units','normalized',...
%                     'HorizontalAlignment','Center','VerticalAlignment','Middle',...
%                     'FontWeight','Bold','FontUnits','normalized','BackgroundColor',...
%                     'none','FontSize',1.25,'Color',dColor,'Visible','off');
% 
%                 for a = 1:nClusters
%                     plot(xAxisInMsec,X{elec}(:,class_idx{elec,a}),'Color',plotColors{a})
%                 end
% 
%                 axis([xAxisInMsec(1) xAxisInMsec(end) yAxisMin(elec) yAxisMax(elec)]);
%                 set(gca, 'Color', get(gcf, 'Color'));
%                 if     j ==  1
%                     set(gca,'xtick',[])
%                 elseif i == 10
%                     set(gca,'ytick',[])
%                 else
%                     set(gca,'xtick',[],'ytick',[])
%                 end
% 
%                 text('Tag',['elec',num2str(elec),'plotD'],'Position',[0.8 0.2],...
%                     'String',num2str(d(elec),'%2.2g'),'Units','normalized',...
%                     'HorizontalAlignment','Center','VerticalAlignment','Middle',...
%                     'FontWeight','Bold','FontUnits','normalized','BackgroundColor',...
%                     [1 0.7 0.7],'FontSize',0.15,'Color',[1 0 0],'Visible','off');
%             end
% %             set(findobj(gcf,'Tag','txt.status'),'String',['Plotting Electrode data... ',num2str((i-1)*10+j),'%']);
% %             drawnow;
%         end
%     end
%     
%     % d', Elec #, Axes, Traces Visibility (based on button status)
%     hObject = {'showBoxD'; 'showBoxNum'; 'showBoxAxes'; 'showBoxTraces'};
%     for m = 1:4
%         eval(sprintf('%s',['showorhide(f1h.',hObject{m},');']));
%     end
% end
% % END Plot all electrodes function


% %% ELECSUBPLOTCREATE.M: Create Subplot for electrode function
% function calculated = elecSubplotCreate(fignum,tagname,A,calculated,elecs)
%     if nargin < 5
%        elecs = A';
%        elecs = elecs(~isnan(elecs))';
%     end
%     
%    for elec = elecs
% %     for i = 1:10
% %         for j = 1:10
% %             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
%             [i,j] = ind2sub(size(A),find(A==elec));
%             if ~isnan(elec)
%                 subplot('Position',...
%                     [Xoffset+(j-1)*xoffset 1-Yoffset-(i-1)*yoffset xsize ysize],...
%                     'Tag',['elec',num2str(elec),tagname{fignum}],'DrawMode','fast');
%                 hold on;
%                 
%                 set(gca, 'Color', get(gcf, 'Color'));
%                 if     j ==  1% && fignum == 1
%                     set(gca,'xtick',[])
%                 elseif i == 10% && fignum == 1
%                     set(gca,'ytick',[])
%                 else
%                     set(gca,'xtick',[],'ytick',[])
%                 end
%             end
% %         end
%     end
%     
%     calculated(fignum).D      = 0;
%     calculated(fignum).Num    = 0;
%     calculated(fignum).Traces = 0;
%     calculated(fignum).Axes   = 0;  % make sure this is always 0
% end
% % END Create Subplot for electrode function


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


% %% Plot electrode Sortnev Checkboxes function
% function plotElecSortnevBoxes(fignum,hObject,eventdata)
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
%             if ~isnan(elec)
%                 h = findobj('Tag',['elec',num2str(elec),tagname{fignum}]);
%                 set(gcf,'CurrentAxes',h(1));
% 
%                 % TODO: Create box which opens SORTNEV
% %                 f1h.SortnevBox(elec) = uicontrol('String','PCA (New Window)',...
% %                     'Style','checkbox','Tag',['SortnevBox',num2str(elec)],'Value',0,'Enable','on',...
% %                     'Units','normalized','BackgroundColor',get(f1h.fig.expnev1,'Color'),...
% %                     'Position',[0.81 0.525 0.2 0.05],'Parent',h(1),...
% %                     'callback',@pcafigure);
%             end
%         end
%     end
% end
% % END Plot electrode Sortnev Checkboxes function


%% Showorhide function
function showorhide(hObject,eventdata)
% How SHOWORHIDE works: First, test if any data has been loaded already. If
% yes, do nothing (return), if no move on.  Second, prepare strings to be
% executed: 
%   (#1) test if a given feature (e.g., d') has already been calculated, and
%   determine value of its checkbox.  if not calculated, and checkbox is
%   checked, run plotElec* and changed calculated flag to 1.  if already
%   calculated, turn visibility on or off (#2)
%   Also, run code (#1) when re-scaling trace graphs
    switch gcf
       case f1h.fig.expnev1
          fignum = 1;
       case f2
          fignum = 2;
       otherwise
          fignum = 1;
    end

    if ~isempty(nevopen_outcome)  % Test whether any data has been loaded or not
        set(findobj(gcf,'Tag','txt.status'),'String','Refreshing view, please wait...');
        drawnow;

        if get(hObject,'Value')
            b = 1; % on_or_off->'on'
        else
            b = 2; % on_or_off->'off';
        end
        
        % PREPARE CODE: find all handles tagged 'elec(#)plot(suffix)'.
        % traces' handles: all children of 'elec(#)plot', except elec(#)plotD and
        % elec(#)plotNum
        object_suffix = get(hObject,'Tag');
        object_suffix = object_suffix(8:end);
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
                set(findobj(gcf,'Tag','txt.status'),'String','Error: Bad ''Show'' selection.');
        end
        
        % RUN CODE
        %(#1)
        if (eval(sprintf(['~calculated(fignum).',object_suffix])) && b == 1) || ~isempty(scaling_params.fit) % (b=1 is 'on')
            eval(sprintf(['calculated(fignum).',object_suffix, ' = 1;']));
            switch object_suffix
               case {'D'; 'Num'}
                   eval(sprintf(['plotElec',object_suffix,'(fignum,tagname,A,d);']));
               case 'Axes'
                   eval(sprintf(['calculated = plotElec',object_suffix,'(nevopen_outcome,calculated);']));
               case {'Traces'; 'PCA'}
                   eval(sprintf(['plotElec',object_suffix,'(fignum,tagname,A,X,',...
                      'iXraw,Xraw,xAxisInMsec,plotColors,k,class_idx,threshold,',...
                      'yAxisMin,yAxisMax);']));
               otherwise
                   set(findobj(gcf,'Tag','txt.status'),'String','Error: Bad ''Show'' selection.');
            end
        elseif eval(sprintf(['calculated(fignum).',object_suffix]))% && b == 2
%             'showorhide(hObject);\n',...

            %(#2)
            for a = 1:96
                eval(cmd);
                set(hObject,'Visible',on_or_off{b});
            end
            
            if b == 2
                eval(sprintf(['calculated(fignum).',object_suffix, ' = 1;']));
            end

        end

%         drawnow;
        set(findobj(gcf,'Tag','txt.status'),'String','');
    end
end
% END showorhide function

% %% showBoxD callback function
% function showBoxD_Callback(hObject, eventdata, handles)
% % hObject    handle to showBoxD (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
%     if nClusters > 1
%         b = 1; % on_or_off->'on';
%     elseif nClusters == 1
%         b = 2; % on_or_off->'off'
%     end
%     
%     set(hObject,'Enable',on_or_off{b});
%     showorhide(hObject);
% end
% % END showBoxD callback function


%% ClusteringNButton callback function
function ClusteringNButton_Callback(hObject, eventdata)
    if ~isempty(nevopen_outcome)
        d         = zeros(96,1);
        class_idx =  cell(96,4);  %CLASS_IDX: set second dimension to be max nClusters possible
        elec      = [];
        
        % If figure 2 is open, notify user it will be redrawn
        if ~isempty(f2)
            set(eval(sprintf(['f',num2str(2),'h.txt.status'])),'String','Refreshing view, please wait...');
            drawnow;
        end
        
        if size(X{1},2) ~= nwaves
            process_elecs;
        end

        [class_idx,d,f1h,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
         class_idx,A,iXraw,f1h,calculated,V);

        set(findobj(gcf,'Tag','txt.status'),'String','Plotting Electrode data...');
        
        calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated);
        drawnow;
%         plotAllElecs;
        showorhide(f1h.showBoxNum);
        showorhide(f1h.showBoxTraces);
        showorhide(f1h.showBoxD);
        showorhide(f1h.showBoxAxes);
        drawnow;
        set(findobj(gcf,'Tag','txt.status'),'String','');
        
        % If figure 2 is open, redraw it
        if ~isempty(f2)
            pcafigure(hObject,eventdata);
        end
    end
end
% END ClusteringNButton callback function


%% ClusteringNPopup callback function
function ClusteringNPopup_Callback(hObject, eventdata)
    nClusters = get(f1h.ClusteringNPopup,'String');
    nClusters = nClusters{get(f1h.ClusteringNPopup,'Value')};
    nClusters = str2double(nClusters);
    try k = nClusters*ones(96,1); catch, k = zeros(96,1); end
    
    nwaves = get(f1h.ClusteringNWavesPopup,'String');
    nwaves = str2double(nwaves{get(f1h.ClusteringNWavesPopup,'Value')});% # of waves
    if nClusters == 1
        set(f1h.showBoxD,'Value',0);
        set(f1h.showBoxD,'Enable','off');
    else
%         set(f1h.showBoxD,'Value',1);
        set(f1h.showBoxD,'Enable','on');
    end
    
    ClusteringNButton_Callback(hObject, eventdata);
end
% END ClusteringNPopup callback function


%% edThresh callback function
function edThresh_Callback(hObject, eventdata)
   
   threshold = str2double(get(f1h.edThresh,'String'))*ones(96,1);
   
   for i = 1:size(A,1)
      for j = 1:size(A,2)
         elec = A(i,j);
         if ~isnan(elec)
            [X{elec},iXraw{elec}] = changeThreshold(threshold(elec),Xraw{elec},f1h);
            
            if ~isempty(X{elec})
               yMax(elec) = max(max(X{elec}(:,1:end)));
               yMin(elec) = min(min(X{elec}(:,1:end)));
            else
               yMax(elec) =  1;
               yMin(elec) = -1;
            end
         end
         set(findobj(gcf,'Tag','txt.status'),'String',...
            sprintf('Rethresholding data... %i%%',...
            100*sub2ind(size(A'),j,i)/numel(A)));
         drawnow;
      end
   end
   
   [class_idx,d,f1h,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
     class_idx,A,iXraw,f1h,calculated,V);
   
   calculated = elecSubplotDestroyChildren(1,tagname,A,nClusters,calculated);
   drawnow;
   %         plotAllElecs;
   showorhide(f1h.showBoxNum);
   showorhide(f1h.showBoxTraces);
   showorhide(f1h.showBoxD);
   showorhide(f1h.showBoxAxes);
   drawnow;
   set(findobj(gcf,'Tag','txt.status'),'String','');

   % If figure 2 is open, redraw it
   if ~isempty(f2)
      pcafigure(hObject,eventdata);
   end
end
% END edThresh callback function


%% save998 callback function
function save998_Callback(hObject, eventdata)
   set(findobj(gcf,'Tag','txt.status'),'String','Saving 998s...');
   drawnow;

%    threshold = str2double(get(f1h.edThresh,'String'))*ones(96,1);
   ichans = unique(A(~isnan(A)))';
   create998(animaltype,iseries,iexp,ichans,threshold);
   
   set(findobj(gcf,'Tag','txt.status'),'String','');
end
% END save998 callback function


%% Scaling Buttons function
function scaling_button_press(hObject, eventdata)
    % sp = scaling_params as defined in the process_elecs function
    sp = scaling_params;
    if ~isempty(nevopen_outcome)
        if strcmp(get(eventdata.NewValue,'String'),'Fit')
            yAxisMin = sp.yMin-0.1*(sp.yMax-sp.yMin);
            yAxisMax = sp.yMax+0.1*(sp.yMax-sp.yMin);
            % TODO: hide y-axis ticks
        elseif strcmp(get(eventdata.NewValue,'String'),'Global')
            yAxisMin = (sp.yMinALL-0.1*sp.ryALL)*ones(size(sp.yMin));
            yAxisMax = (sp.yMaxALL+0.1*sp.ryALL)*ones(size(sp.yMax));
            % TODO: show y-axis ticks
        else
            warning('Invalid Scaling Option');
        end


        set(findobj(gcf,'Tag','txt.status'),'String','Plotting Electrode data...');
        drawnow;
%         plotAllElecs;
%         showorhide(f1h.showBoxNum);

        scaling_params.fit = 1;
        showorhide(f1h.showBoxTraces);
        scaling_params.fit = [];
        
        showorhide(f1h.showBoxD);
        showorhide(f1h.showBoxAxes);
        drawnow;
        set(findobj(gcf,'Tag','txt.status'),'String','');
    end
end
% END Global Scaling Button function


%% Second Figure (PCA graphs) function
function pcafigure(hObject,eventdata)
    if ~isempty(f2) && strcmp(get(hObject,'Tag'),'PCABox')
        calculated = elecSubplotDestroyChildren(2,tagname,A,nClusters,calculated);
        cancel_f2;
    else %if isempty(f2) && strcmp(get(hObject,'Tag'),'PCABox') || strcmp(get(hObject,'Tag'),'ClusteringNButton')
        % Initialize Figure
        if strcmp(get(hObject,'Tag'),'PCABox')
            f2 = figure('Position',[520 380 672 420],'Name','explorenev6-PCA',...
                'Toolbar','figure','MenuBar','none','NumberTitle','off');
            % Create the 'status' text area
            f2h.status = uicontrol('Style','text','Tag','status','Units','normalized',...
                'Position',[0.8 0 0.2 0.075],'BackgroundColor',get(f2,'Color'));
            f2h.showpanel = uipanel('Title','Show','Tag','showpanel',...%'Visible','off',...
                'Units','normalized','Position',[0.8 0.225 0.2 0.15],...
                'BackgroundColor',get(f2,'Color'));
            f2h.showBoxNum = uicontrol('Parent',f2h.showpanel,'String','Elec #',...
                'Style','checkbox','Tag','showBoxNum','Value',0,...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[11 4 60 23],...
                'callback',@showorhide);
            f2h.showBoxAxes = uicontrol('Parent',f2h.showpanel,'String','Axes',...
                'Style','checkbox','Tag','showBoxAxes','Value',1,...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[71 24 60 23],...
                'callback',@showorhide);
            f2h.showBoxPCA = uicontrol('Parent',f2h.showpanel,'String','Clusters',...
                'Style','checkbox','Tag','showBoxPCA','Value',1,...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[71 4 60 23],...
                'callback',@showorhide);
            f2h.showBoxD = uicontrol('Parent',f2h.showpanel,'String','d''',...
                'Style','checkbox','Tag','showBoxD','Value',1,'Enable','on',...
                'Units','pixels','BackgroundColor',get(f2,'Color'),'Position',[11 24 40 23],...
                'callback',@showorhide);

            hObject = f1h.PCABox; figure(f2);
            calculated = elecSubplotCreate(         2,tagname,A,  calculated);
        else
            hObject = f1h.PCABox; figure(f2);
            calculated = elecSubplotDestroyChildren(2,tagname,A,nClusters,calculated);
        end
        
        if ~isempty(nevopen_outcome)
            showorhide(f2h.showBoxNum);
            plotElecPCA(2);
            showorhide(f2h.showBoxD);
            showorhide(f2h.showBoxAxes);
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

end