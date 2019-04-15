function [f1h,nevopen_outcome,animaltype,fname,nevsorted,mwaves,SamplingRateInKHZ,...
   nelecs,nSamplesPerPacket] = opennevfile(f1h,fname,fdir,nevsorted)
% OPENNEVFILE.M: Checks to see you have the right nevfile, or asks you to 
% pick one, then opens it

global DIRS

if exist('nevsorted','var') && ~isempty(nevsorted)
   mychoice = questdlg(['You have sorted units which you will lose ',...
      'when loading a new file. Do you want to continue loading?'], ...
      'Loading a nev file', 'No');
   switch mychoice
      case {'No', 'Cancel'}
         return;
   end
end
% set(gcbf,'WindowButtonDownFcn','zoom(gcbf,''down'')');  %% set to zoom...

if isfield(f1h,'chkbx') && isfield(f1h.chkbx,'unitsave')
   for i = 1:4
      set(f1h.chkbx.unitsave(i),'Value',0);
   end
end

%     fname = get(f1h.txt.nevfile, 'String');
if isempty(fname)
   %AZ 20081209: was
   %[fname,p] = uigetfile('*.nev','Select a NEV file');
   % added default directory/file, outcome
   [fname,fdir,outcome] = uigetfile('*.nev','Select a NEV file',[DIRS.Cerebus ...
      filesep 'CATZ077' filesep 'u004_022.nev']);
   %     [fname,p,outcome] = uigetfile('*.nev','Select a NEV file','..\..\Data\Cerebus\');
   if ~outcome
      set(findobj(gcf,'Tag','txt.status'),'String','No nev file selected.  Please try again.');
      return;
   end
%    from explorenev6 (this seems unnecessary)
else
   if strcmpi(fname(1:2),'C:')
      [fdir fname ext] = fileparts(fname);
      fname = [fname ext];
      clear ext;
   end
%    fdir  = [DIRS.Cerebus filesep animaltype filesep];
%    fname = sprintf('u%03d_%03d.nev',series,iexp);
end

nevsorted = {};
mwaves = {};  %% mean waveforms
% from sortnev2: unnecessary?
% elec = str2double(get(f1h.txtbox.elecnum,'String'));
% %     set(f1h.pshbtn.nextelec,'String','Electrode #');
% set(f1h.txtbox.elecnum,'String',num2str(elec));

set(gcf,'Pointer','watch');
set(findobj(gcf,'Tag','txt.status'),'String','Opening nev file...');
drawnow;
fprintf('Opening nev file %s\n', fname);

% ADDED by AB to work locally (07-08-08)
%---------------------------------------
if fname
   animaltype = fdir(end-strfind(fliplr(fdir),fliplr(DIRS.Cerebus))+3:end-1);
   if isempty(animaltype)
      if fdir(end) == filesep
         fdir(end) = [];
      end
      animaltype = fdir(find(fdir==filesep,1,'last')+1:end);
%       animaltype = fdir(17:end);
   end
   name = fname(1:regexp(fname,'\.nev')-1);
%    name = name{1};
end

xst = exist(['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'], 'file');
if ~xst
   nevfile = dir([fdir fname]);
   if isempty(nevfile)
      error('%s\\%s is not a valid nev file.',animaltype,fname);
   end
   reply = questdlg(['Do you want to copy the source nev-file (',...
      num2str(round(nevfile.bytes/(1024^2))),'MB) to C:\WINDOWS\Temp?'],...
      'sortnev2: Copy to Temp folder?', 'Yes');
   if strcmp(reply,'Cancel')% exit
      set(findobj(gcf,'Tag','txt.status'),'String','File open canceled.');
      set(f1h.fig.f1,'Pointer','arrow');
      nevopen_outcome = [];
      drawnow;
      return;
   end
   if strcmpi(reply,'Yes')
      mkdir(  ['C:\WINDOWS\Temp\',animaltype,filesep]);
      copyfile([fdir fname], ...
         ['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'],'f');
      fname = ['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'];
% from explorenev6: unnecessary?
%    elseif strcmpi(reply,'No')
%       fname = [fdir fname];
   end
   %AZ20080422: If file is in Temp, don't ask user, just use it
else%if xst && ~strcmpi(fname,['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'])
   % reply = input(...
   % 'Found a local copy in the Temp of C:\WINDOWS, do you want to use it? y/n [y]: ', 's');
   % if isempty(reply), reply = 'y'; end
   % if strcmpi(reply,'y')
   fname = ['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'];
   % end
end%
%--------------------------------------

[nevopen_outcome,SamplingRateInKHZ,nelecs,nSamplesPerPacket] = nevopen(fname);

if ~nevopen_outcome
%    nevclose;
   fclose(pepNEV.index.fid);
   set(findobj(gcf,'Tag','txt.status'),'String','Could not open file!');
   set(gcf,'Pointer','arrow');
   return;
   %     else
   %         set(f1h.txt.nevfile,'String',fname);
end