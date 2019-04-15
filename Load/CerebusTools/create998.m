function create998(animaltype,iseries,iexp,ichans,threshold)
% CREATE998.M: Loads an experiment, then creates a MUA 998 unit for 
% supra-threshold traces
% Usage: 
% create998(animaltype,iseries,iexp,ichans,threshold)
% or
% create998(animaltype,iseries,iexp,ichans,[old_iseries old_iexp])
% 
% Examples:
% create998('CATZ077',4,26)
% No threshold value assumes default value of 32.
% No ichans argument assumes you want to create a 998 for all electrodes.
% Alternatively, you can specify which electrodes to create 998's for:
% create998('CATZ077',4,26,1     )
% create998('CATZ077',4,26,1:9,32)
% Or specify the threshold only:
% create998('CATZ077',4,26,[] ,32)
% 
% Also, you can pass a 2-D vector to threshold, specifying an earlier
% [series experiment] with saved 998s that you want to use for this iseries
% and iexp:
% create998('CATZ077',4,26,1:9,[4 25])
% create998('CATZ077',4,26,[] ,[4 25])
% 
% 
% Created by AZ 2009-07-20
% 
% SEE ALSO: CHANGETHRESHOLD, SAVESUPRATHRESHOLD, SORTNEV2
SetDefaultDirs;

if nargin < 5 || isempty(threshold)
   threshold = 32; %default value
end
if nargin < 4 || isempty(ichans)
   ichans = 1:96;
end

% if threshold is a 2-D vector, replicate thresholds from saved 998 units
% from a specified [series experiment]
switch size(threshold,2)
   case 1     
      switch size(threshold,1)
         case 1 % expand threshold to 96-D vector
            threshold   = threshold*ones(96,1);
         case 96
            % do nothing
         otherwise
            error('create998:BadInput','if threshold argument is 1-D, it should have 1 or 96 elements');
      end
   case 2     % load thresholds from 998s
      old_iseries = threshold(1);
      old_iexp    = threshold(2);
      threshold   = zeros(96,1);
      
      for elec = ichans
         old_file = [DIRS.spikes filesep animaltype filesep num2str(old_iseries) ...
            filesep animaltype '_s' sprintf('%02d',old_iseries) '_e' ...
            sprintf('%02d',old_iexp) '_c' num2str(elec+1000) '_u998.mat'];
         if exist(old_file,'file')
            unit = load(old_file);
            threshold(elec) = unit.unit.gmm.threshold;
            clear unit
         else
            fprintf('998 unit does not exist for series %d exp %d elec %2d\n',...
               [old_iseries old_iexp elec]);
            threshold(elec) = 0;
         end
      end
      
   otherwise
      error('create998:BadInput','threshold argument should be a 1- or 2-D vector');
end

SetDefaultDirs;

fdir  = [DIRS.Cerebus filesep animaltype filesep];
fname = sprintf('u%03d_%03d.nev',iseries,iexp);

%% From SORTNEV2
OSName = computer;
switch OSName
   case {'PCWIN','PCWIN64'}
      UserName = getenv('UserName');
   case {'GLNX86','GLNXA64','MACI','SOL64'}
      UserName = getenv('USER');
   otherwise
      %Leave usrtxt at default value of 'dario'
      UserName = 'dario';
end
% end SORTNEV2

%% from LOADNEV
if fname
   animaltype = fdir(end-strfind(fliplr(fdir),fliplr(DIRS.Cerebus))+3:end-1);
   name = regexp(fname,'\.nev','split');
   name = name{1};
end

xst = exist(['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'], 'file');
if ~xst
   nevfile = dir([fdir fname]);
   reply = questdlg(['Do you want to copy the source nev-file (',...
      num2str(round(nevfile.bytes/(1024^2))),'MB) to C:\WINDOWS\Temp?'],...
      'sortnev2: Copy to Temp folder?', 'Yes');
   if strcmp(reply,'Cancel')% exit
      warning('create998:FileOpenCanceled','File open canceled.');
      return;
   end
   if strcmpi(reply,'Yes')
      mkdir(['C:\WINDOWS\Temp\',animaltype,filesep]);
      copyfile([fdir fname], ['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'],'f');
      fname = ['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'];
   end
   %AZ20080422: If file is in Temp, don't ask user, just use it
else
   fname = ['C:\WINDOWS\Temp\',animaltype,filesep,name,'.nev'];
end
nevopen_outcome = nevopen(fname);

if ~nevopen_outcome
   error('create998:FileOpenError','Could not open file!');
end
% End LOADNEV

for elec = ichans
   if threshold(elec)
      % Test if file already exists, and if so, is threshold the same
      existingUnit = sprintf('%s_s%02d_e%02d_c%d_u%03d.mat',...
         [DIRS.spikes filesep animaltype filesep num2str(iseries) filesep animaltype],...
         iseries,iexp,elec+1000,998);
      xst = exist(existingUnit,'file');
      if xst
         load(existingUnit);
      end
      
      if ~xst || unit.gmm.threshold ~= threshold(elec)
         
         fprintf('Elec #%02d: Loading data... (threshold will be set to %2d)\n',...
            [elec threshold(elec)])
         
         %% from LOADELEC
         nwaves = [];%AZ20090220% = 2500;
         % nwaves = [1 2500];
         
         gmm = [];
         if ~isempty(nwaves)
            [T, waveforms] = nevwaves(elec,nwaves);  %% Get the first nwaves
         else
            [T, waveforms] = nevwaves(elec);             %% Get ALL waveforms
         end
         
%          N = size(waveforms,2); % N waveforms
         Xraw = waveforms;
         
         if isempty(Xraw)
            error('create998:NoData','No data!  Electrode skipped!');
         end
         
         [Xrealigned,iXraw] = changeThreshold(threshold(elec),Xraw,[],nevopen_outcome);
         
         % RESIZE ( take only samples (7:34) (t = -0.1:1/30:0.8) )
         if size(Xrealigned,1) == 48
            Xrealigned = Xrealigned(7:34,:);
         end
         
         if isempty(Xrealigned)
            fprintf('Elec #%02d: Threshold (%2d) too high!\n',elec,threshold(elec))
            % End LOADELEC
         else
            %% SAVE 998
            saveSupraThreshold(animaltype,elec,iseries,iexp,gmm,UserName,T,...
               Xrealigned,iXraw,threshold(elec));
         end
      else
         fprintf('Elec #%02d: 998 unit with same threshold (%2d) already exists\n',...
            elec,threshold(elec));
      end % end if unit already exists
   end
end % end for all elecs

fprintf('**************** ¡¡¡DONE!!! ****************\n')