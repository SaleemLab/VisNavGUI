% SORTEDUNITCHECK.M: check for previously sorted units
% Works by polling 'series' subdirectories of a given animal, pulls info
% from filenames and stores them in a matrix, 'already_sorted'
% Variables needing prior initialization: DIRS, animaltype, iseries (optional)

% CHECK FOR PREVIOUSLY SORTED UNITS
sortedUnitsFound = false;
currentStatusObj = findobj('-regexp','Tag','txt.status','-and','Parent',gcf);

while ~sortedUnitsFound
   if exist('fh','var')
      set(currentStatusObj,'String','Checking for sorted units...'); drawnow;
   end
   % already_sorted   = [];
   animal_spikesdir = [];

   if ~exist('iseries','var') || isempty(iseries) || exist('fh','var')
      seriesdirs = dir([DIRS.spikes filesep animaltype]);

      if isempty(seriesdirs)
         break;
      end
      % don't poll 'archive' directory
      if strcmp(seriesdirs(end).name,'archive')
          seriesdirs = seriesdirs(1:end-1);
      end
   else
      seriesdirs(3).name = num2str(iseries);
   end

   % Run 'dir' on all 'series' directories
   for i = 3:length(seriesdirs)
      if exist('fh','var') % for sortnev2 gui
         set(fh.expnev1.txt.status,'String',sprintf('Checking for sorted units... %2.0f%%',... 
            (i-2)/(length(seriesdirs)-2)*100/3)); drawnow; 
      end
      allseries = str2double(seriesdirs(i).name); % to make sure we have a (numeric) series dir
      animal_spikesdir     = [ animal_spikesdir; ...
         what([DIRS.spikes filesep animaltype filesep num2str(allseries)]) ];
   %   animal_spikesdir(1).name = [DIRS.spikes,filesep,animaltype,filesep,num2str(allseries)];
   end
   % Use only the .mat files
   animal_spikesdir = vertcat(animal_spikesdir(:).mat);
   % Split filenames into series, exp, elec#, and unit info

   % Matrix with cols: iseries | iexp | elec | unit | dir_index
   already_sorted   = zeros(size(animal_spikesdir,1),5); % pre-allocate
   parts = cell(size(animal_spikesdir,1),4);

   a = regexprep(animal_spikesdir(:),'^CAT[Z]{0,1}[0-9]{2,4}_s|[\.mat]{4}$','');
   parts = regexp(a,'_[a-z]','split');
   for i = 1:size(parts,1)
   %     [pathstr, name, ext, versn] = fileparts(animal_spikesdir{i});
   %     if strcmp(ext,'.mat')
   %     set(fh.expnev1.txt.status,'String',['Checking for sorted units... ' ...
   %         num2str(50 + ((i-2)/(length(animal_spikesdir)-2)*100/2)) '%']); drawnow;
   %         parts = regexp(animal_spikesdir{i},'_[a-z]','split');

           if parts{i}{1}(1) <= '9'
   %             parts{i}{5} = parts{i}{5}(1:length(parts{i}{5})-4);
               already_sorted(i,:) = [ str2double(parts{i}{1}) ...
                                       str2double(parts{i}{2}) ...
                                       str2double(parts{i}{3}) ...
                                       str2double(parts{i}{4}) ...
                                       i];
           end
   %     end
      if exist('fh','var') && floor(i/2000) == i/2000
         set(currentStatusObj,'String',sprintf('Checking for sorted units... %2.0f%%',... 
            ( 1 + (i-2)*2/(length(animal_spikesdir)-2) )*100/3)); drawnow;
      end
   end

   % % Matrix with cols: iseries | iexp | elec | unit | dir_index
   % already_sorted   = zeros(size(animal_spikesdir,1),4); % pre-allocate
   % a = regexprep(animal_spikesdir(:),'^CAT[Z]{0,1}[0-9]{2,4}_s|[\.mat]{4}$','');
   % clear b;
   % for i = 1:size(animal_spikesdir,1)
   %    if a{i}(1) <= '9'
   %       b(i,:) = regexp(a{i},'_[a-z]','split');
   %       if floor(i/1000) == i/1000
   %          set(fh.expnev1.txt.status,'String',['Checking for sorted units... ' ... 
   %             num2str(( 1 + (i-2)*2/(length(animal_spikesdir)-2) )*100/3) '%']); drawnow;
   %       end
   %    end
   % end
   % already_sorted = [str2double(b) (1:size(animal_spikesdir,1))'];
   % 
   % delete extra entries
   already_sorted = already_sorted( ~isnan(already_sorted(:,4)),:);
   already_sorted = already_sorted(logical(already_sorted(:,4)),:);

   % Load any saved 998 thresholds
   if exist('iseries','var')
      already_sorted_998 = already_sorted(sum([already_sorted(:,1) == iseries ...
                                               already_sorted(:,2) == iexp    ...
                                               already_sorted(:,4) == 998        ],2)==3,:);
   end

   % Don't care about 999 or 998 (?) units
   already_sorted = already_sorted(already_sorted(:,4) < 100   ,:);

   % already_sorted = flipud(already_sorted); % AZ 2009-02-18: prioritize later experiments

   % AZ 20090406: Don't care about fancy numbering anymore
   % takenUnitNums = unique(already_sorted(:,4));
   % takenUnitNums = sort(takenUnitNums);
   % a = 1:max(takenUnitNums);
   % b = zeros(1,a(end));
   % b(takenUnitNums) = a(takenUnitNums);
   % % nextUnitNum   = max(already_sorted(:,4)) + 1;
   % nextUnitNum   = min(find(a-b));
   % if isempty(nextUnitNum)
   %     nextUnitNum = max(takenUnitNums)+1;
   % end

   % % Don't care about other files
   % already_sorted = already_sorted(find(already_sorted(:,1)==iseries),:);
   % already_sorted = already_sorted(find(already_sorted(:,2)==iexp   ),:);
   sortedUnitsFound = true;
end