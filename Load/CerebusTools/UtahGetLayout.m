function [arrayLayout, cableEntry] = UtahGetLayout(animal, series, arraySerialNumber)
% UtahGetLayout - gives the physical arrangement of Utah array channels in an NxN matrix
%
% [matrix, site] = UtahGetLayout(animal, iseries)
% 
% [matrix, site] = UtahGetLayout(animal, iseries, arraySerialNumber)
%
% It returns a matrix with the layout of the channels, and the site where
% the cable is connected to the array
%
% If the argument arraySerialNumber is not provided, the function might
% prompt the user to specify the serial number of the array that was used in a
% particular series. Here's a list of the available options:
% 
% Array Serial Number       Used in 
% 1017-0010                 10x10 array, cat87 s11-13, cat86 s10-s13, cat84 s1-s9, cat81 s4
% 1017-0011                 10x10 array, cat83 s5, cat82 all series, cat81 s3
% 1017-0012                 10x10 array, cat88 all series, cat87 s1-10, cat86 s1-s9, cat85 s4-s10, cat84 s10-s16, cat83 s6-s12, cat80-cat75 all series
% 1017-0019                 10x10 array, cat74-cat72 all series
% 1017-0047                 5x5 array, 0.5 mm long, only used in mice so far
% 1017-0048                 4x5 array, 0.5 mm long, only used in mice so far
%
%
% 2008-Feb-21 SK
% 
% 2008-May SK added check if we can access the machine on which the cerebus data reside
% 2008-Apr LB added the use of DIRS
% 2008-Aug LB added check for variable cableEntry (old saving format was without cableEntry)
% 2010-Feb AZ added 1017-0048 array

global DIRS;

nos = {'1017-0010','1017-0011','1017-0012','1017-0019','1017-0047','1017-0048'};
extNos = {'1017-0010   (cat81, cat84, cat86-87)', '1017-0011   (cat81-83)', ...
   '1017-0012   (cat75-80, cat83-88)', '1017-0019   (cat72-74)', ...
   '1017-0047   (5x5, 0.5 mm long)',  '1017-0048   (4x5, 0.5 mm long)'};

if nargin < 2
    error('Invalid number of arguments.');
end

if nargin == 2
    arrayNo = [];
end

if nargin == 3 && ~ismember(arraySerialNumber, nos)
    error('Invalid Array Serial Number');
elseif nargin == 3 && ismember(arraySerialNumber, nos)
    arrayNo = arraySerialNumber;
end

dataDir = strcat(DIRS.Cerebus, '\', animal, '\');
matFileName = strcat(dataDir, sprintf('layoutSeries%03d.mat', series)); 

% check if we have access to the machine
tmp = ls(DIRS.Cerebus);
if isempty(tmp)
    error('Cannot access %s', DIRS.Cerebus);
end

% check if a file with the layout already exists
if exist(matFileName, 'file')
    load(matFileName);
    if ~exist('cableEntry', 'var')
        arrayNo = [];
    else
        return;
    end
end

% if the *.mat file is not found in the data directory and the serial number is not provided as an argument,
% get the serial number from user

if isempty(arrayNo)
    [index, ok] = listdlg('PromptString', 'Select a serial number for the Utah array:',...
                    'ListSize', [300 75], 'SelectionMode','single', 'ListString', extNos);
    if ~ok
        error('No Array Serial Number specified.');
    else
        arrayNo = nos{index};
    end
end

% determine arrangement of banks
switch arrayNo
    
    case '1017-0010'  

        A = [60 70 29 39 89 99 48 58 7 17 67 77 26 36 86 96 45 ...
            55 4 14 64 74 23 33 83 93 42 52 11 21 71 81];

        B = [50 40 19 9 79 69 38 28 98 88 57 47 16 6 76 66 35 ...
            25 95 85 54 44 13 3 73 63 32 22 92 82 61 51];

        C = [20 30 80 90 49 59 8 18 68 78 27 37 87 97 46 56 5 ...
            15 65 75 24 34 84 94 43 53 2 12 62 72 31 41];

        eNo = [C B A]; % reflects how the connectors are plugged in into the preamp
        cableRow = 6; cableCol = 10;
        utahSites = flipud(reshape(1:100,10,10)'); % site numbering as specified in the data sheet
        arrayLayout = ones(size(utahSites)) * NaN;
    
    

    case '1017-0011'  % used in cats 81-82

        A = [60 70 29 39 89 99 48 58 7 17 67 77 26 36 86 96 45 ...
            55 4 14 64 74 23 33 83 93 42 52 11 21 71 81];

        B = [50 40 19 9 79 69 38 28 98 88 57 47 16 6 76 66 35 ...
            25 95 85 54 44 13 3 73 63 32 22 92 82 61 51];

        C = [20 30 80 90 49 59 8 18 68 78 27 37 87 97 46 56 5 ...
            15 65 75 24 34 84 94 43 53 2 12 62 72 31 41];

        eNo = [C B A]; % reflects how the connectors are plugged in into the preamp
        cableRow = 6; cableCol = 10;
        utahSites = flipud(reshape(1:100,10,10)'); % site numbering as specified in the data sheet
        arrayLayout = ones(size(utahSites)) * NaN;


    case '1017-0012'  % used in cats 75-79

        A = [60 70 29 39 89 99 48 58 7 17 67 77 26 36 86 96 45 ...
            55 4 14 64 74 23 33 83 93 42 52 11 21 71 81];
        B = [50 40 19 9 79 69 38 28 98 88 57 47 16 6 76 66 35 ...
            25 95 85 54 44 13 3 73 63 32 22 92 82 61 51];
        C = [20 30 80 90 100 59 8 18 68 78 27 37 87 97 46 56 5 ...
            15 65 75 24 34 84 94 43 53 2 12 62 72 31 41];

        eNo = [C B A]; % reflects how the connectors are plugged in into the preamp
        cableRow = 6; cableCol = 10;
        utahSites = flipud(reshape(1:100,10,10)'); % site numbering as specified in the data sheet
        arrayLayout = ones(size(utahSites)) * NaN;


    case '1017-0019' % used in cat72, cat73, cat74)

        A = [60 70 29 39 89 99 48 58 7 17 67 77 26 36 86 96 45 ...
            55 4 14 64 74 23 33 83 93 42 52 11 21 71 81];
        B = [50 40 19 9 79 69 38 28 98 88 57 47 16 6 76 66 35 ...
            25 95 85 54 44 13 3 73 63 32 22 92 82 61 51];
        C = [20 30 80 90 49 59 8 18 68 78 27 37 87 97 46 56 5 ...
            15 65 75 24 34 84 94 43 53 2 12 62 72 31 41];

        eNo = [C B A];
        cableRow = 6; cableCol = 10;
        utahSites = flipud(reshape(1:100,10,10)'); % site numbering as specified in the data sheet
        arrayLayout = ones(size(utahSites)) * NaN;


    case '1017-0047' % used in mice

        eNo = 1:25; % AZ 2010-02: officially, should be 1:32 
        cableRow = 3; cableCol = 1;
        utahSites = reshape(1:25,5,5); % site numbering as specified in the data sheet
        arrayLayout = ones(size(utahSites)) * NaN;
        
    case '1017-0048' % used in mice

        eNo = 1:20; % AZ 2010-02: officially, should be 1:32 
        cableRow = 3; cableCol = 1;
        utahSites = reshape(1:20,4,5); % site numbering as specified in the data sheet
        arrayLayout = ones(size(utahSites)) * NaN;

    otherwise
        error('Unknown Array Number: %s.', arrayNo);
end


% for ieNo = 1 : length(eNo)
%     [row, col] = find(utahSites == eNo(ieNo));
%     arrayLayout(row, col) = ieNo;
% end
% AZ20090604: Vectorized above code below
arrayLayout(eNo) = 1:size(eNo,2);
arrayLayout      = rot90(arrayLayout,1);

cableEntry = arrayLayout(cableRow, cableCol);


% make sure target directory exists (this is helpful for online analyses,
% in which case the data might not have been copied to zeon yet)
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

save(matFileName, 'arrayLayout', 'cableEntry');

