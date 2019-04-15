function sortnev2(varargin)
%      SORTNEV opens a gui for spike sorting of utah array data
%
%      SORTNEV(myAnimal,mySeries,myExpt) opens the gui and loads the associated .nev file
%      after the user hits the 'load' button
%
%      ############ modified by AB ################
%      Each time a unit is sorted from a channel, save it. A .nevsorted file will
%      be created with the electrode no. To continue, reload the .nev and
%      go to the next electrode of interest.
%
%      SORTNEV, by itself, creates a new SORTNEV or raises the existing
%      singleton*.
%
%      H = SORTNEV returns the handle to a new SORTNEV or the handle to
%      the existing singleton*.
%
%      SORTNEV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SORTNEV.M with the given input arguments.
%
%      SORTNEV('Property','Value',...) creates a new SORTNEV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sortnev_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sortnev_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Change log
%   Created by Dario
%   07-2008 AB modified
%   08-2008 LB modified help text
%              added option for loading specified .nev file
%   08-2008 LB made saving more flexible
%   08-2008 LB added a warning when loading new data in the presence of
%       sorted units
%   12-2008 AZ added default directory/file (local), handles cancelling
%       loading a file gracefully. Set usertxt to username logged into
%       computer.  added clustering status message.  made prettier.
%   01-2009 AZ added d' calculation, which treats all non-selected points
%       for a given unit as noise.  added 'unc' buttons to select all
%       remaining unclassified points
%   01-2009 AZ ported code to SORTNEV2: a non-guide version of sortnev

% Edit the above text to modify the response to help sortnev

%% HINT: If you can't close a window because of a CloseRequestFcn error, run:
%%       delete(get(0,'Children'))
global f1h

col    = {[1 0 0] [0 0 1] [0 1 0] [1 0 1] [0.5 0.5 0.5]};

% % Initialize figure window
% f1h = sortnev2_fig_init(f1h,col);
% % sortnev2_fig_init;


% Create a variable to store the handle of the new GUI figure
% f2 becomes non-empty when the new figure is created.
f1h.fig.f2 = [];

%% Begin Main Figure Main Function
% toggle: 'Menu'/'Toolbar' = 'none'/'figure'
f1h.fig.f1 = figure('Position',[520 170 654 783],'Name','sortnev2',...
    'Toolbar','none','MenuBar','none','NumberTitle','off',...
    'Color',get(0,'defaultUicontrolBackgroundColor'),...
    'CloseRequestFcn', @closefirstfigure);

% Main Controls (Top of window)
f1h.pshbtn.loadnev = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@loadnev_Callback,'KeyPressFcn',@shortcutKeys,'Visible','off',...
'ListboxTop',0,'String','Load','Tag','pshbtn.loadnev',...
'Position',[0.180 0.961 0.106 0.029]);
f1h.pshbtn.save = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@save_Callback,'KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','Save','Tag','pshbtn.save',...
'Position',[0.073 0.961 0.093 0.029]);
f1h.pshbtn.reload = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@loadelec,'KeyPressFcn',@shortcutKeys,...
'Position',[0.180 0.961 0.106 0.029],...
'String','Reload','Tag','pshbtn.elecnum');

f1h.txt.elecnum = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'HorizontalAlignment','center','KeyPressFcn',@shortcutKeys,...
'Position',[0.073 0.928 0.093 0.022],...
'String','Electrode #','Style','text','Tag','pshbtn.elecnum');
f1h.pshbtn.prevelec = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@elecrement_Callback,'KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','<','Tag','pshbtn.prevelec',...
'Position',[0.180 0.928 0.030 0.028]);
f1h.txtbox.elecnum = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'BackgroundColor',[1 1 1],'FontSize',9,...
'Callback',@elecnum_Callback,...%'KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','1','Style','edit','Tag','txtbox.elecnum',...
'Position',[0.218 0.928 0.030 0.028]);
f1h.pshbtn.nextelec = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@elecrement_Callback,'KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','>','Tag','pshbtn.nextelec',...
'Position',[0.256 0.928 0.030 0.028]);

f1h.txt.username = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'HorizontalAlignment','right','KeyPressFcn',@shortcutKeys,...
'Position',[0.297 0.961 0.058 0.022],...
'String','User:','Style','text','Tag','txt.username');
f1h.txtbox.username = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'HorizontalAlignment','center',...%'KeyPressFcn',@shortcutKeys,...
'Position',[0.364 0.961 0.101 0.029],...
'String','User:','Style','edit','Tag','txtbox.username');


% f1h.txt.cluster = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
% 'String','Cluster','Tag','txt.cluster',...
% 'FontSize',10,'FontWeight','bold','HorizontalAlignment','left',...
% 'ListboxTop',0,'String','Cluster:','Style','text','Tag','txt.cluster',...
% 'Position',[0.567 0.942 0.075 0.029]);
f1h.pshbtn.addClearEllipses = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'String','Add Ellipse','Tag','pshbtn.addClearEllipses',...
'Position',[0.567 0.961 0.115 0.029],...
'Callback',@addClearEllipses_Callback,'KeyPressFcn',@shortcutKeys);

f1h.panel.clusterAuto = uipanel('Parent',f1h.fig.f1,'Units','normalized',...
'Tag','panel.clusterAuto','Title','remove','Visible','on',...
'BorderType','line','HighlightColor',[0.7 0.7 0.7],...
'Position',[0.695 0.956 0.112 0.037]);
f1h.pshbtn.clusterAuto = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'String','Auto','Tag','pshbtn.clusterAuto',...
'Position',[0.700 0.961 0.050 0.029],...
'Callback',@clusterAuto_Callback,'KeyPressFcn',@shortcutKeys);
f1h.pop.nclusters = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Position',[0.755 0.961 0.045 0.027],'BackgroundColor',[1 1 1],...
'String',{ '1'; '2'; '3'; '4' },'Style','popupmenu','Value',2,'Tag','pop.nclusters',...
'KeyPressFcn',@shortcutKeys);
% f1h.txt.centroids = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
% 'FontSize',10,'FontWeight','bold','HorizontalAlignment','left',...
% 'KeyPressFcn',@shortcutKeys,...
% 'Position',[0.815 0.928 0.100 0.029],...
% 'String','Centroids','Style','text','Tag','txt.centroids');

f1h.pshbtn.showall = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@showall_Callback,'Enable','on','KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','Show All','Tag','pshbtn.showall',...
'Position',[0.860 0.961 0.077 0.029]);
f1h.pshbtn.showpca = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@showpca_Callback,'Enable','on','KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','PCA','Tag','pshbtn.showpca',...
'Position',[0.945 0.961 0.040 0.029]);

f1h.txt.status = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'FontSize',10,'FontWeight','bold','HorizontalAlignment','left',...
'ListboxTop',0,'String','','Style','text','Tag','txt.status',...
'Position',[0.075 0.890 0.390 0.026]);
% f1h.txt.nevfile = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
% 'FontSize',10,'FontWeight','bold','HorizontalAlignment','left',...
% 'ListboxTop',0,'String','','Style','text','Tag','txt.nevfile',...
% 'Position',[0.075 0.965 0.875 0.026]);
f1h.txt.nchan = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'FontSize',10,'FontWeight','bold','HorizontalAlignment','left',...
'ListboxTop',0,'String','','Style','text','Tag','txt.nchan',...
'Position',[0.075 0.861 0.390 0.026]);


% Unclassified plot (left side)
f1h.axes.unclassified = axes('Parent',f1h.fig.f1,'Units','normalized',...
'Position',[0.075 0.506 0.396 0.287],'Color',get(f1h.fig.f1,'Color'),...
'XColor',get(0,'defaultaxesXColor'),'YColor',get(0,'defaultaxesYColor'),...
'Tag','unclassified','Visible','on');


% PCA plots on right side
f1h.axes.pc(1) = axes('Parent',f1h.fig.f1,'Units','normalized',...
'Position',[0.567 0.640 0.388 0.275],'Color',get(f1h.fig.f1,'Color'),'Box','off',...
'XColor',get(0,'defaultaxesXColor'),'YColor',get(0,'defaultaxesYColor'),...
'Tag','pc(1).ax','Visible','on','XColor',[0.4 0   0],'YColor',[0 0.4 0  ]);
f1h.axes.pc(2) = axes('Parent',f1h.fig.f1,'Units','normalized',...
'Position',[0.567 0.335 0.388 0.275],'Color',get(f1h.fig.f1,'Color'),'Box','off',...
'XColor',get(0,'defaultaxesXColor'),'YColor',get(0,'defaultaxesYColor'),...
'Tag','pc(2).ax','Visible','on','XColor',[0.4 0   0],'YColor',[0 0   0.4]);
f1h.axes.pc(3) = axes('Parent',f1h.fig.f1,'Units','normalized',...
'Position',[0.567 0.030 0.388 0.275],'Color',get(f1h.fig.f1,'Color'),'Box','off',...
'XColor',get(0,'defaultaxesXColor'),'YColor',get(0,'defaultaxesYColor'),...
'Tag','pc(3).ax','Visible','on','XColor',[0   0.4 0],'YColor',[0 0   0.4]);


%% UNIT CONTROLS
f1h.panel.unitcontrols = uipanel('Parent',f1h.fig.f1,'Units','normalized',...
'Tag','panel.unitcontrols','Title','Sorting','Visible','on',...
'Position',[0.045 0.040 0.450 0.425]);
f1h.txt.unitnum       = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','Unit #',     'Style','text','Tag','txt.unitnum',...
'Position',[0.025 0.945 0.250 0.045],'BackgroundColor','r');
f1h.txt.unitdprime    = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','d'' =',     'Style','text','Tag','txt.unitdprime',...
'Position',[0.025 0.880 0.250 0.045]);
f1h.txt.saveunit      = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','SAVE',     'Style','text','Tag','txt.saveunit',...
'Position',[0.025 0.815 0.250 0.045]);
f1h.txt.showtrace     = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','Traces','Style','text','Tag','txt.showtrace',...
'Position',[0.025 0.750 0.250 0.045]);
f1h.txt.meanwave      = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','Mean','Style','text','Tag','txt.meanwave',...
'Position',[0.025 0.685 0.250 0.045]);
f1h.txt.bringtotop    = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','Bring to Top',  'Style','text','Tag','txt.bringtotop',...
'Position',[0.025 0.620 0.250 0.045]);
f1h.btngrp.bringtotop = uibuttongroup('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'BorderType','none','visible','on','Tag','btngrp.bringtotop',...
'SelectionChangeFcn',@highlightPointTrace,...
'Position',[0.295 0.620 0.680 0.045]);
f1h.txt.oldMeanWaves  = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
'FontSize',9,'FontWeight','normal','HorizontalAlignment','left',...
'ListboxTop',0,'String','Old Means','Style','text','Tag','txt.oldMeanWaves',...
'Position',[0.025 0.555 0.250 0.045]);

for i = 1:5
   f1h.txtbox.unitnum(i)     = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
   'FontSize',9,'Enable','off',...
   'ListboxTop',0,'String','','Style','edit','Tag',['txtbox.unitnum(',num2str(i),')'],...
   'Position',[0.155+i*0.140 0.945 0.120 0.050],'BackgroundColor',col{i});

   f1h.txt.unitdprimeVal(i)  = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
   'FontSize',9,'FontWeight','normal','HorizontalAlignment','center',...   
   'ListboxTop',0,'String','',        'Style','text',...
   'Tag',['txt.unitdprimeVal(',num2str(i),')'],...
   'Position',[0.155+i*0.140 0.880 0.120 0.050]);

   f1h.chkbx.unitsave(i)     = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
   'Style','checkbox','Tag',['chkbx.unitsave(',num2str(i),')'],'Value',0,'Enable','off',...
   'Position',[0.190+i*0.140 0.815 0.050 0.045],'BackgroundColor',get(f1h.fig.f1,'Color'));

   f1h.chkbx.showtraces(i)   = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
   'Style','checkbox','Tag',['chkbx.showtraces(',num2str(i),')'],'Value',0,'Enable','off',...
   'Position',[0.190+i*0.140 0.750 0.050 0.045],'BackgroundColor',get(f1h.fig.f1,'Color'),...
   'Callback',@showTraces);

   f1h.chkbx.meanwave(i)     = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
   'Style','checkbox','Tag',['chkbx.meanwave(',num2str(i),')'],'Value',0,'Enable','off',...
   'Position',[0.190+i*0.140 0.685 0.050 0.045],'BackgroundColor',get(f1h.fig.f1,'Color'),...
   'Callback',@showMeanTrace);

   f1h.radio.bringtotop(i)   = uicontrol('Parent',f1h.btngrp.bringtotop, 'Units','normalized',...
   'Style','radio',   'Tag',['radio.bringtotop(',num2str(i),')'],'Value',0,'Enable','off',...
   'Position',[(((115-50-6)/560)/2)+((137.5+2.5)/680)*(i-1) 0.000 50/680 1.000]);

   f1h.chkbx.oldMeanWaves(i) = uicontrol('Parent',f1h.panel.unitcontrols,'Units','normalized',...
   'Style','checkbox','Tag',['chkbx.oldMeanWaves(',num2str(i),')'],'Value',0,'Enable','off',...
   'Position',[0.190+i*0.140 0.555 0.050 0.045],'BackgroundColor',get(f1h.fig.f1,'Color'),...
   'Callback',@showMeanTrace);
end

% AZ20090220: Thresholding tool
f1h.txt.thresh = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Visible','on',...
'ListboxTop',0,'String','Threshold: +/-','Style','text','Tag','txt.thresh',...
'Position',[0.180 0.008 0.140 0.026]);
f1h.txtbox.thresh = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'BackgroundColor',[1 1 1],'FontSize',9,'Visible','on',...
'Callback',@loadelec,...%'KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','32','Style','edit','Tag','txtbox.thresh',...
'Position',[0.320 0.008 0.040 0.028]);
f1h.pshbtn.save998 = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
'Callback',@save998_Callback,'KeyPressFcn',@shortcutKeys,...
'ListboxTop',0,'String','Save 998','Tag','pshbtn.save998',...
'Position',[0.370 0.008 0.093 0.029]);

% f1h.chkbx.fixedsc = uicontrol('Parent',f1h.fig.f1,'Units','normalized','Enable','off',...
% 'String','Fixed Scaling','Style','radiobutton','Value',1,'Tag','chkbx.fixedsc',...
% 'Position',[0.038 0.011 0.132 0.019],'KeyPressFcn',@shortcutKeys);%,...
% f1h.chkbx.interact = uicontrol('Parent',f1h.fig.f1,'Units','normalized',...
% 'String','Interactive Plotting','Style','checkbox','Value',0,'Tag','chkbx.interact',...
% 'Position',[0.200 0.011 0.200 0.019],'Callback',@interactivePlotToggle,...
% 'KeyPressFcn',@shortcutKeys);

% f1h.plotInteractive = false;
f1h.highlights = [];

% FOR ELLIPSES CREATED IN plotGMMellipses.m
% Define a context menu; it is not attached to anything
f1h.menu.pcaellipse(1) = uicontextmenu;
% Define the context menu items and install their callbacks
f1h.menu.pcaellipse(2) = uimenu(f1h.menu.pcaellipse(1),...
   'Callback', @clickEllipse_Callback,...
   'Label','edit ellipse');% #',...
             %regexp(get(gco),'[0-9]','match','once') ]);
f1h.menu.pcaellipse(3) = uimenu(f1h.menu.pcaellipse(1),...
   'Callback', @clickEllipse_Callback,...
   'Label','remove ellipse');


%% CREATE FUNCTIONS: username_CreateFcn (+ elecnum_CreateFcn + nclusters_CreateFcn)
% AZ 20081216: Sets username to username logged into computer
OSName = computer;
switch OSName
    case {'PCWIN','PCWIN64'}
        UserName = getenv('UserName');
        
        set(f1h.txtbox.elecnum,'BackgroundColor','white'); % from elecnum_CreateFcn
        
        if isequal(get(f1h.txtbox.username,'BackgroundColor'), ...
                     get(0,'defaultUicontrolBackgroundColor'))
                   set(f1h.txtbox.username,'BackgroundColor','white');
        end
        
        %% nclusters_CreateFcn
        if isequal(get(f1h.pop.nclusters,'BackgroundColor'), ...
                   get(0,'defaultUicontrolBackgroundColor'))
                   set(f1h.pop.nclusters,'BackgroundColor','white');
        end
    case {'GLNX86','GLNXA64','MACI','SOL64'}
        UserName = getenv('USER');
        
        set(f1h.txtbox.elecnum,'BackgroundColor', ...
         get(0,'defaultUicontrolBackgroundColor')); % from elecnum_CreateFcn
    otherwise
        %Leave usrtxt at default value of 'dario'
        UserName = 'dario';
end
set(f1h.txtbox.username,'String',UserName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %% KEYPRESS FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(f1h.fig.f1,'KeyPressFcn',@shortcutKeys);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             %% GUI MAIN FUNCTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION
global DIRS serverName pepNEV
DIRS = struct('data',[],'spikes',[],'camera',[],'xfiles',[],'michigan',[],...
    'Cerebus',[],'stimInfo',[]);
serverDataDir = [];
SetDefaultDirs;

nevsorted         = [];
mwaves            = []; fname             = [];
elec              = []; already_sorted    = []; already_sorted_to_load = [];
animal_spikesdir  = []; animaltype        = []; already_sorted_998     = [];

iseries           = []; iexp              = [];
seriesdirs        = [];

sortedUnitsFound = false;
currentStatusObj = [];

S = []; T = []; U = []; V = []; X = []; C = []; k = []; d = [];
Xraw = [];  iXraw = []; Xrealigned  = [];

prev              = []; gmm               = struct;
% gmm.obj = gmdistribution(0,1,1);
%uidx              = [];
sidx              = cell(1,5);
f1h.plots         = []; unit              = struct;
nevopen_outcome   = []; threshold         = [];
idx_tograph       = []; sidx_tograph      = cell(1,5);
idx_tograph_only  = []; %sidx_tograph_only = cell(1,5);
%fdialog           = []; fdialogh          = [];

nelecs = [];
SamplingRateInKHZ = 30;
nSamplesPerPacket = 48;
xAxisInMsec = (1:nSamplesPerPacket)/SamplingRateInKHZ;
threshold = 32;

sortedUnitsFound = [];

% Temporary variables
xr = []; yr = []; ry = [];
My = []; my = [];
a  = []; b  = []; i  = []; j = [];
allseries   = []; parts = [];

% AZ20090502 from: LB 080808 to allow automatic file specification
if length(varargin) >= 3
   animaltype = varargin{1};
   iseries    = varargin{2};
   iexp       = varargin{3};
   if length(varargin) >= 4 % Set elec #
      elec = varargin{4};
      set(f1h.txtbox.elecnum,'String',num2str(elec));
      if length(varargin) == 5 % Set thresholds
         threshold = varargin{5};
         set(f1h.txtbox.thresh,'String',num2str(threshold));
      end
   end
   
   fdir  = [DIRS.Cerebus filesep animaltype filesep];
   fname = sprintf('u%03d_%03d.nev',iseries,iexp);
   
%    set(f1h.txt.nevfile, 'String', [fdir fname]);
   set(f1h.fig.f1,'Name',['sortnev2: ',animaltype,', Series ',num2str(iseries),...
      ', Experiment ',num2str(iexp)]);
   loadnev_Callback;
end

set(f1h.txtbox.thresh,'String',num2str(threshold));

% varargout = {gmm};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %% CALLBACK FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOADNEV_CALLBACK CALLBACK FUNCTION
function loadnev_Callback(hObject, eventdata)
    % hObject    handle to loadnev_Callback (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

%     global fname nev elec X sidx nevsorted  hfile mwaves DIRS already_sorted animal_spikesdir;

[f1h,nevopen_outcome,animaltype,fname,nevsorted,mwaves,SamplingRateInKHZ,...
   nelecs,nSamplesPerPacket] = opennevfile(f1h,fname,fdir,nevsorted);

xAxisInMsec = (1:nSamplesPerPacket)/SamplingRateInKHZ;

figure(f1h.fig.f1);
sortedUnitCheck;

loadelec;   %% load electrode data
end
% END LOADNEV_CALLBACK function


%% ELECNUM_CALLBACK CALLBACK FUNCTION
function elecnum_Callback(hObject, eventdata)
% hObject    handle to elecnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%     [sidx,gmm,nevsorted,mwaves] = applyandsave(sidx,gmm,nevsorted,mwaves);  %% apply rules and save

   elecNew = str2double(get(f1h.txtbox.elecnum,'String'));  %% get new electrode
   if ~isnan(elecNew)
      elec = fix(elecNew);
      set(f1h.txtbox.elecnum,'String',num2str(elec));
   else % Display error if given a non-integer input, replace with previous value
      set(f1h.txtbox.elecnum,'String',num2str(elec));
      
      angry;
      return;
   end

   for i = 1:4
      set(f1h.chkbx.unitsave(i),'Value',0);
   end

   if nevopen_outcome
      loadelec;
   end % end if nevopen_outcome
end
% END ELECNUM_CALLBACK CALLBACK FUNCTION


%% NEXTELEC CALLBACK FUNCTION
% --- Executes on button press in elecrement.
function elecrement_Callback(hObject, eventdata)
% hObject    handle to elecrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%     [sidx,gmm,nevsorted,mwaves] = applyandsave(sidx,gmm,nevsorted,mwaves);  %% apply rules and save

%     elec = elec+1;
   set(f1h.txt.nchan,'String',sprintf(['Total sorted units = ' num2str(length(nevsorted))]));
   if     gco == f1h.pshbtn.nextelec && elec < 96
      set(f1h.txtbox.elecnum,'String',num2str(elec+1));
   elseif gco == f1h.pshbtn.prevelec && elec > 1
      set(f1h.txtbox.elecnum,'String',num2str(elec-1));
   end

   elecnum_Callback(hObject,eventdata);
%     for i = 1:4
%         set(f1h.chkbx.unitsave(i),'Value',0);
%     end
% 
%     loadelec;

end
% END NEXTELEC CALLBACK FUNCTION


%% CLUSTERAUTO_CALLBACK CALLBACK FUNCTION
function clusterAuto_Callback(hObject, eventdata)
% hObject    handle to clusterAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   gmm = [];
   [f1h,gmm,d,sidx,sidx_tograph,C,S,U,V,Xrealigned] = cluster_gmm3sd( ...
    f1h,gmm,d,sidx,sidx_tograph,  S,U,V,Xrealigned,                   ...
    iXraw,idx_tograph,xAxisInMsec,SamplingRateInKHZ,col,threshold,prev    );
end
% END CLUSTERAUTO_CALLBACK CALLBACK FUNCTION


%% ADDCLEARELLIPSES_CALLBACK CALLBACK FUNCTION
function addClearEllipses_Callback(hObject, eventdata)
% hObject    handle to addClearEllipses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   if ~isempty(sidx{4}) && strcmp(get(f1h.pshbtn.addClearEllipses,'String'),'Clear Ellipses')
      clearEllipses;
      gmm = struct; k = []; sidx = cell(1,5); sidx_tograph = cell(1,5);
%       set(f1h.pshbtn.addClearEllipses,'String','Add Ellipse');

      f1h = sortnev2_controls_update(false,f1h);
      [Xrealigned,idx_tograph,idx_tograph_only,sidx_tograph,sidx,f1h] = ...
          sortnev2_gui_update(Xrealigned,V,iXraw,sidx_tograph,gmm,prev,xAxisInMsec,col,...
          sidx,d,threshold,f1h,'Ellipses Cleared.');
      return;
   end
   
   [f1h,gmm,k] = createPCellipse(f1h,gmm,col);

   if isempty(gmm) || isempty(k)
      return;
   end
   
   if k == 4
      set(f1h.pshbtn.addClearEllipses,'String','Clear Ellipses');
   end
      
	set(f1h.pop.nclusters,'Value',size(gmm.mu,1)); drawnow;
   
   [f1h,gmm,d,sidx,sidx_tograph,C,S,U,V,Xrealigned] = cluster_gmm3sd( ...
    f1h,gmm,d,sidx,sidx_tograph,  S,U,V,Xrealigned,                   ...
    iXraw,idx_tograph,xAxisInMsec,SamplingRateInKHZ,col,threshold,prev    );
end
% END ADDCLEARELLIPSES_CALLBACK CALLBACK FUNCTION;


%% --- Executes on button press in save.
function save_Callback(hObject, eventdata)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(f1h.fig.f1,'Pointer','watch'); drawnow;

%% PRE-SAVE CHECKS
% AZ 2009-03-17
% If a save unit's d' is higher than current d', ask user to confirm save
if exist('prev','var') && sum(size(prev)) > 1 && isfield(prev(1),'replace_dprime')
   
   answer = questdlg(['Would you like to replace the previously sorted units,', ...
      ' or create new ones?  Note: old units will be moved to ''..', ...
      filesep,'archive',filesep,'''.'], ...
      'sortnev2: Replace previously sorted units?','Replace','Cancel','Replace');

   if strcmp(answer,'Replace') && max(d) < max(prev(1).replace_dprime)
      answer = questdlg({['A previously saved unit for this electrode has a higher d'' (',...
         num2str(max(prev(1).replace_dprime)),').'];...
         'Are you sure you want to replace it?'},...
         'sortnev2: Replace previously saved unit with higher d''?',...
         'Replace','Cancel','Replace');
   end
   
   if strcmp(answer,'Replace')
      % Check if you will be overwriting a unit another user saved
      for i = 1:size(prev,2)
         if prev(i).unit.iexp == iexp && ~isempty(prev(i).unit.source) && ...
                                          ~strcmp(prev(i).unit.source,UserName)
            prevUser = prev(i).unit.source;
         end
      end
      if exist('prevUser','var') && ~isempty(prevUser)
         answer = questdlg({['The unit you are replacing was saved by another user (',...
            prev(1).unit.source,').'];...
            'Are you sure you want to replace it?'},...
            ['sortnev2: Replace unit saved by ',prev(1).unit.source,'?'],...
            'Replace','Cancel','Replace');
      end
%       prev(1).replace        = currUnitNum;
   end
   if strcmp(answer,'Cancel')
      set(findobj(gcf,'Tag','txt.status'),'String','Save cancelled.');
      set(f1h.fig.f1,'Pointer','arrow');
      drawnow;
      return;
   end
end
      
[sidx,gmm,nevsorted,mwaves] = applyandsave(sidx,gmm,nevsorted,mwaves);  %%

anywave =      get(f1h.chkbx.unitsave(1),'Value') | ...
               get(f1h.chkbx.unitsave(2),'Value') | ...
               get(f1h.chkbx.unitsave(3),'Value') | ...
               get(f1h.chkbx.unitsave(4),'Value') ;  % at least one!

if ~anywave || isempty(nevsorted)
    set(findobj(gcf,'Tag','txt.status'),'String','Noting to save! First select unit.');
    set(f1h.fig.f1,'Pointer','arrow');
    return;
end

gmm.threshold = threshold(1);
%% SAVE UNIT(S)
units = nev2unit(animaltype,elec,iseries,iexp,gmm,nevsorted,mwaves,UserName);
% Save 998
saveSupraThreshold(animaltype,elec,iseries,iexp,gmm,UserName,T,...
   Xrealigned,iXraw,threshold);

%% POST-SAVE VARIABLE UPDATING
% Reset unit checkboxes
for i = 1:4
    set(f1h.chkbx.unitsave(i),'Value',0);
end

% AZ20090331: TODO: FIX ANY PROBLEMS ARISING WHEN WE STOP USING THE CODE BELOW
% AZ20090211: Instead of polling directories after each save, manually
% insert saved units into already_sorted & animal_spikesdir variables
if exist('prev','var') && sum(size(prev)) > 1 && ...
      ( ~isfield(prev(1),'replace') || isfield(prev(1),'loaded') )
   % DID WE ADD PREVIOUSLY SORTED UNITS, USING SAME OR DIFF PARAMS
   %   if so, 
   % xr,yr,ry = temporary variables (loaded,saved,comparison)
   xr = already_sorted_to_load(:,1:4);
   yr = [];
   if ~isempty(xr)
      yr = [ones(nnz(gmm.icell),1)*[iseries iexp 1000+elec] ...
         gmm.icell(nonzeros(gmm.icell))'];
      
      % Put units from this exact experiment first:
      i = find(xr(:,2) == iexp)';
      %Reorder
      xr = xr([i setdiff(1:size(xr,1),i)],:);
      
      if size(xr,1) == size(yr,1)
         ry = (xr == yr);
      else
         ry = [1 zeros(1,3)];
      end
   end
end

if exist('ry','var') && all(sum(ry,2)==4)
   % do nothing (we have replaced existing units)
   set(findobj(gcf,'Tag','txt.status'),'String','Re-sorted data saved!');
   set(f1h.fig.f1,'Pointer','arrow');
   return;
else
   % if creating new units, add them to our list
   i = find(any(ry==0,2));
%    if ~isempty(i)
   already_sorted = [ [ones(length(i),1)*[iseries iexp 1000+elec] ...
      gmm.icell(i)' length(animal_spikesdir)+(1:length(i))']; ...
      already_sorted ];
%    else
   i = setdiff(find(gmm.icell),i);
   already_sorted = [ [ones(length(i),1)*[iseries iexp 1000+elec] ...
      gmm.icell(i)' length(animal_spikesdir)+(1:length(i))']; ...
      already_sorted ];
%    end

   nonempty_icells = find(gmm.icell);
   for i = 1:nnz(gmm.icell)
      name = regexp(UnitGetFilename(animaltype,iseries,iexp,elec+1000,...
         gmm.icell(nonempty_icells(i))),filesep,'split');
      animal_spikesdir{end+i} = name{3};
   end
   % now, don't need to call
   % sortedUnitCheck;

   set(findobj(gcf,'Tag','txt.status'),'String','Re-sorted data saved!');
   set(f1h.fig.f1,'Pointer','arrow');         
end
end


%% --- Executes on button press in save.
function save998_Callback(hObject, eventdata)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(f1h.fig.f1,'Pointer','watch'); drawnow;

saveSupraThreshold(animaltype,elec,iseries,iexp,gmm,UserName,T,...
   Xrealigned,iXraw,threshold);

set(findobj(gcf,'Tag','txt.status'),'String','Thresholded MUA data saved!');
set(f1h.fig.f1,'Pointer','arrow');         
end


function showall_Callback(hObject, eventdata)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   set(f1h.fig.f1,'Pointer','watch'); drawnow;
   xr =               get(f1h.pop.nclusters,'String');  %temporary var
   k  = str2double(xr{get(f1h.pop.nclusters,'Value' )});

   save_Debug;
   set(f1h.fig.f1,'Pointer','arrow');
end

function showpca_Callback(hObject, eventdata)
   set(f1h.fig.f1,'Pointer','watch'); drawnow;
   f1h = PCA_Debug(f1h,xAxisInMsec,U,S);
   set(f1h.fig.f1,'Pointer','arrow');
end


% AZ 20090602: avoid variable scope issues by calling external function here
function clickEllipse_Callback(hObject,eventdata)
   switch hObject 
      case f1h.menu.pcaellipse(2), action =   'edit';
      case f1h.menu.pcaellipse(3), action = 'remove';
      otherwise
         error('sortnev2:clickEllipse:BadSelection','Bad Right Click Selection');
   end
   
   [f1h,gmm,d,sidx,sidx_tograph,C,S,U,V,Xrealigned] = clickEllipse( ...
    f1h,gmm,d,sidx,sidx_tograph,  S,U,V,Xrealigned,                 ...
    iXraw,idx_tograph,xAxisInMsec,SamplingRateInKHZ,col,threshold,prev,action);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          %% NON-CALLBACK FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOADELEC FUNCTION
function loadelec(hObject,eventdata)

   elec = str2double(get(f1h.txtbox.elecnum,'String'));  %% get new electrode
%    xAxisInMsec = (0:47)/SamplingRateInKHZ;
   nwaves = [];%AZ20090220% = 2500;
   % nwaves = [1 2500];

   gmm = [];
   sidx         = cell(1,5);
   sidx_tograph = cell(1,5);
   prev         = struct([]);

   %AZ20090107: consolidate code, clear dprime text
   set(f1h.fig.f1,'Pointer','watch'); drawnow;

   f1h = sortnev2_controls_update(false,f1h);

   set(f1h.fig.f1,'CurrentAxes',f1h.axes.pc(1)       ); cla;
   set(f1h.fig.f1,'CurrentAxes',f1h.axes.pc(2)       ); cla;
   set(f1h.fig.f1,'CurrentAxes',f1h.axes.pc(3)       ); cla;
   set(f1h.fig.f1,'CurrentAxes',f1h.axes.unclassified); cla; hold on;
   for field = {'plots'; 'means';}% 'highlights'}
      if       isfield(f1h,field)
         f1h = rmfield(f1h,field);
      end
   end
   clearHighlights;         

   set(findobj(gcf,'Tag','txt.status'),'String','Loading data... Please wait!');
   drawnow;

   %% nev = loadnev(fname,elec,'all');  %% load chan

   % skip re-loading/thresholding if reload button pushed
   if ~(exist('hObject','var') && (hObject == f1h.pshbtn.reload || hObject == f1h.txtbox.thresh))
      nevopen_outcome = nevopen(fname);
      if ~isempty(nwaves)
         [T, waveforms] = nevwaves(elec,nwaves);  %% Get the first nwaves
      else
         [T, waveforms] = nevwaves(elec);             %% Get ALL waveforms
      end
%       nevclose;
      fclose(pepNEV.index.fid);

      N = size(waveforms,2); % N waveforms
      Xraw = waveforms;

   %    %% remove big transients
   %    try
   %       p = prctile(X',[2 98])';
   %       idx = find(sum((X < repmat(p(:,1),[1 size(X,2)])) + (X > repmat(p(:,2),[1 size(X,2)]))) > 20);
   %       X(:,idx)  = [];
   %    catch
   %    end

      if isempty(Xraw)
         set(findobj(gcf,'Tag','txt.status'),'String','No data!  Electrode skipped!');
         set(f1h.fig.f1,'Pointer','arrow');
         return;
      else
         iXraw      = false(size(Xraw,2),1);
         Xrealigned = zeros(size(Xraw));
         
         %% AZ20090724: change threshold to 998's
         %% TODO: resave 998 when saving other units
         if ~isempty(already_sorted_998(already_sorted_998(:,3)==elec+1000,5))
            unit = load([DIRS.spikes filesep animaltype filesep num2str(iseries) filesep ...
               animal_spikesdir{already_sorted_998(already_sorted_998(:,3)==elec+1000,5)}]);
            threshold = unit.unit.gmm.threshold;
            set(f1h.txtbox.thresh,'String',num2str(threshold));
            clear unit
         end
         
         rethresh;
      end
   else
      %TODO: make less permissive
      if threshold ~= str2double(get(f1h.txtbox.thresh,'String')) || isempty(Xrealigned)
         rethresh;
      end
   end % end load/threshold only if reload button not just pushed

   if isempty(Xrealigned)
      set(findobj(gcf,'Tag','txt.status'),'String','Threshold too high!');
      return;
   end

   %AZ20090204: LOAD SAVED PARAMETERS
   already_sorted_to_load = already_sorted(already_sorted(:,3)==elec+1000,:);
   sorted_series  = already_sorted_to_load(:,1);
%    currUnitNum    = already_sorted_to_load(:,4);
   prevSavedUnits = already_sorted_to_load(:,5);
   
%    % Put units from this exact experiment first:
%    unit.iexp      = already_sorted_to_load(:,2);
%    i = find(unit.iexp == iexp)';

%    %Reorder
%    if i
%       sorted_series  =  sorted_series([i setdiff(1:size(prevSavedUnits,1),i)],:);
%       currUnitNum    =    currUnitNum([i setdiff(1:size(prevSavedUnits,1),i)],:);
%       prevSavedUnits = prevSavedUnits([i setdiff(1:size(prevSavedUnits,1),i)],:);
%    end
   
   % AZ 20090406
   if ~exist('hObject','var') || (exist('hObject','var') && hObject ~= f1h.txtbox.thresh)
       [Xrealigned,V,xAxisInMsec,sidx,d,threshold,f1h,iXraw,gmm,prev,S,U,C,k] = ...
          prevUnitProcess(Xrealigned,gmm,prev,xAxisInMsec,Xraw,iXraw,nevopen_outcome,...
                          SamplingRateInKHZ,already_sorted_to_load,animaltype,...
                          sorted_series,animal_spikesdir,prevSavedUnits,iseries,...
                          iexp,threshold,k,f1h);
   end
   
   [Xrealigned,idx_tograph,idx_tograph_only,sidx_tograph,sidx,f1h] = ...
       sortnev2_gui_update(Xrealigned,V,iXraw,sidx_tograph,gmm,prev,xAxisInMsec,col,...
       sidx,d,threshold,f1h,'Saved parameters loaded successfully.');
end
% END LOADELEC function


function [sidx,gmm,nevsorted,mwaves] = applyandsave(sidx,gmm,nevsorted,mwaves)

   %AZ20081218% blen = 50000;            %% block length in waves...

   xr =               get(f1h.pop.nclusters,'String');  % temporary var
   k  = str2double(xr{get(f1h.pop.nclusters,'Value' )});
   anywave =      get(f1h.chkbx.unitsave(1),'Value') | ...
                  get(f1h.chkbx.unitsave(2),'Value') | ...
                  get(f1h.chkbx.unitsave(3),'Value') | ...
                  get(f1h.chkbx.unitsave(4),'Value') ;  % at least one!

	if anywave && ~isempty(Xrealigned) %AZ20081217 don't accidentally add nothing to nevsorted/mwaves

      %AZ20090107: Change text only if actually applying and saving...
      set(findobj(gcf,'Tag','txt.status'),'String',sprintf('Re-sorting in progress...  Please wait!'));

      if ~isempty(T)   %AZ20081218

% AZ20090129: DEBUG %
% save_Debug;
% AZ20090129: END DEBUG %

         sidxorg = cell(1,5);
         for i = 1:4
            bol = get(f1h.chkbx.unitsave(i),'Value');
            if bol  %&& ~isempty(spoly{i})) %AZ20090204
               sidxorg{i} = T(C==i);  %AZ20090203
            end
         end

      end

      j = 1;
      gmm.icell = zeros(1,k);
      gmm.dprime = zeros(1,k);
        for i = 1:4
            if ~isempty(sidxorg{i})
                nevsorted{i} = {elec sort(sidxorg{i})};
                mwaves{i} =  mean(Xrealigned(:,sidx{i}),2);
%                 gmm.icell(i) = 1;
                gmm.S     = S;             %AZ20090204: save matrices   to reconstruct PCA
                gmm.U     = U;             %AZ20090204: save matrices   to reconstruct PCA
                gmm.dprime(i) = d(i);%AZ20090211: save dprime for a selected unit
%                 if exist('prev','var') && sum(size(prev)) > 1 && isfield(prev(1),'loaded')
%                     gmm.icell(i) = prev(j).loaded(i);
%                     j = j + 1;
%                 elseif exist('prev','var') && sum(size(prev)) > 1 && ...
%                         isfield(prev(1),'replace') && j <= length(prev.replace)
%                     gmm.icell(i) = prev(1).replace(j);
%                     j = j + 1;
%                 else
%                     gmm.icell(i) = j;
%                     j = j + 1;
%                 end
% AZ20090505
%                gmm.icell(i) = str2double(get(f1h.txtbox.unitnum(i),'String'));
               gmm.icell(i) = j;
               j = j + 1;
            end
         end

      %AZ20090107: Change text only if actually applying and saving...
      set(findobj(gcf,'Tag','txt.status'),'String','Re-sorting done!');
      pause(0.5); %AZ20090107: sped up (was 1)
	end

    drawnow;
end


%% SELECTTHREREST2 function (obsolete)
% function selectTheRest2(hObject,eventdata)
%    % SELECTTHEREST2
%    % Selects all unclassified point as unit # i
%    % change log:
%    %   created by Andrew Zaharia 2009-01-07
%    %   AZ 2009-01-29: calculates convex hull, puts points in spoly
%    %   AZ 2009-02-24: adapted for sortnev2
% 
%    j = get(hObject,'Tag');
%    j = str2double(j(13));
%    set(findobj(gcf,'Tag','txt.status'),'String',['Putting all unclassified points in ch',num2str(j)])
% 
%    set(f1h.pop.nclusters,'Value',max(howManyClasses(sidx(1:4)),1));
%    k = get(f1h.pop.nclusters,'Value');
% 
%    sidx{j} = sidx{5};
%    sidx{5} = [];
% 
%    if isempty(vertcat(sidx{setdiff(1:4,j)}))
%       gmm.obj = gmdistribution.fit(V,1,'Start','randSample');
%       C = cluster(gmm.obj,V);
%       for i = 1:4
%          sidx{i} = find(C==i);
%       end
%    end
% %    dim = 1:size(V,2);
%    
%    if k == 1
%       d = zeros(1,4);
%    else
%       for i = 1:4
%          d(i) = dprime(V,sidx{i},vertcat(sidx{setdiff(1:4,i)}));
%       end
%    end
%    
%    sidx_tograph = sidx;
%    for i = 1:5
%       sidx_tograph{i}(setdiff(sidx{i},idx_tograph)) = [];
%    end
% 
%    f1h = plotPCA(        f1h,           V,sidx_tograph,gmm,                 col,sidx);
%    f1h = plotWaveforms(  f1h,Xrealigned,  sidx_tograph,    prev,xAxisInMsec,col,sidx,threshold);
% %    f1h = plotsortedunits(f1h,Xrealigned,  sidx_tograph,    prev,xAxisInMsec,col,sidx,d);
%    
%    % disable 'unc' buttons
%    for i = 1:4
% %        set(f1h.pshbtn.unit(i).unc,'Enable','Off');
%    end
% 
%    set(findobj(gcf,'Tag','txt.status'),'String','');
% end


%% Original figure's closing function
%% HINT: If you can't close a window because of a CloseRequestFcn error, run:
%%       delete(get(0,'Children'))
function closefirstfigure(hObject, eventdata, handles)
    % If the new figure window the new button created is open ...
    if ~isempty(f1h.fig.f2)
        % Close it
        delete(f1h.fig.f2);
    end
    % Close the original figure
    delete(f1h.fig.f1);
end
% End of original figure's closing function

end