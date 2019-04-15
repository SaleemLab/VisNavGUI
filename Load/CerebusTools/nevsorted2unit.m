function unit = nevsorted2unit(inputFile,nextUnitNum)
% nevsorted2unit converts .nevsorted files to units
%
% unit = nevsorted2unit(inputFile)
% inputFile is a file with .nevsorted extension created by sortnev
%
% 07-08-08 AB created
% 2008-08 LB removed input of animal and parses inputFile instead, which
%   contains the animal name
% 2008-08 LB corrected some bugs in the numbering of electrodes and units
% 2008-12 AZ changed unit.source output to be the user's name instead of
% 'NEV file'
% 2009-02 AZ allows extra argument to increment unit.icell if units have
% already been sorted for given animal

if nargin < 2
    nextUnitNum = 1;
end

SetDefaultDirs;

load(inputFile,'-mat');

if isempty(resort.sorted)
    fprintf('<nevsorted2unit> Noting to save. Please select the units to save\n');
    unit = [];
    return;
end
    
[pathstr, name, ext, versn] = fileparts(inputFile); 
[iseries iexpt electrode] = strread(name,'u%f %f %f','delimiter','_');
aidx = strfind(pathstr, 'CATZ'); % assume it's a cat
if isempty(aidx) % it's a mouse
    aidx = regexp(pathstr, 'M\d\d\d\d\d\d'); % yymmdd
end
slashidx = strfind(pathstr, filesep);
animal = pathstr(aidx:slashidx(find(slashidx > aidx, 1))-1);
    
p = ProtocolLoad(animal,iseries,iexpt);

%% Reorganize sync times

SyncTimes = resort.sync.timestamps/30; % double( tt(ii==0) )/30;

nrepeats = floor( length(SyncTimes)/2 / p.nstim );

if nrepeats ~= p.nrepeats
    fprintf('Error! Protocol has %d repeats, Cerebus data has %d repeats. Have to stop here\n',p.nrepeats, nrepeats);
    return;
end

TT1 = zeros(p.nstim,nrepeats);
TT2 = zeros(p.nstim,nrepeats);
for istim = 1:p.nstim
    for irepeat = 1:nrepeats
        ipres = p.seqnums(istim,irepeat);
        TT1(istim,irepeat) = SyncTimes(2*ipres-1);
        TT2(istim,irepeat) = SyncTimes(2*ipres+0);
    end
end

stimdurs = (TT2 - TT1); % in milliseconds

%% Rearrange spike times and save unit structures

nsingUnits = length(resort.sorted);
if nsingUnits > 1
    eNums = [resort.sorted{:}];
    eNums = [eNums{1:2:end}];
    mye = unique(eNums);
    uNums = nan(size(eNums));
    for ie = 1 : length(mye)
        idx = find(eNums == mye(ie));
        uNums(idx) = 1:length(idx);
    end
else
    eNums = resort.sorted{1}{1};
    uNums = 1;
end

%AZ 20081217: Get username from usertxt field (for saving in unit.source)
UserName = get(findobj('Tag','usertxt'),'String');

for isingUnit = 1:nsingUnits

    TheseTimes = resort.sorted{isingUnit}{2}/30; %double( tt(ii==ie) )/30;
    SpikeTimes{isingUnit} = cell(p.nstim, nrepeats);
    for istim = 1:p.nstim
        for irepeat = 1:nrepeats
            jj = (TheseTimes>TT1(istim,irepeat)) & (TheseTimes<=TT2(istim,irepeat));
            SpikeTimes{isingUnit}{istim,irepeat} = ((TheseTimes( jj )-TT1(istim,irepeat))./1000)'; % convert from milliseconds to seconds
        end
    end

    unit = UnitCreate;
    unit.animal = p.animal;
    unit.ichan = eNums(isingUnit)+1000; % electrode number
    %unit.icell =  isingUnit;
    unit.icell = uNums(isingUnit)+nextUnitNum-1;
    unit.iseries = p.iseries;
    unit.iexp =  p.iexp;
    unit.stimdurs = stimdurs./1000;
    unit.timestamp = datestr(now);
    unit.spiketimes = SpikeTimes{isingUnit};
    unit.datatype = 'spiketimes';
%AZ 20081217    unit.source = 'NEV file';
    unit.source = UserName;  %AZ 20081217: save user name from usertxt field
    unit.gmm = resort.gmm;

    UnitSave(unit,DIRS.spikes);

end

%chanlist = resort.sorted{1}{1}+1000;

return


figs = AnalyzeRingach('catz085',6,15, 1 );



