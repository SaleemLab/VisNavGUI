function unitOUT = nev2unit(animaltype,elec,iseries,iexp,gmm,nevsorted,mwaves,UserName)
% nev2unit saves sorted .nev data directly to units
%
% 2009-02-05 AZ Created
% 2009-07-10 AZ Now outputs all units saved in one variable

% global nevsorted mwaves gmm h pepNEV DIRS animaltype iseries iexp elec
global DIRS pepNEV

for isingUnit = find(gmm.icell)

if isempty(nevsorted)
    set(h.status,'String','Nothing to save! First select unit.');
    return;
end

p = ProtocolLoad(animaltype,iseries,iexp);

%% Reorganize sync times

SyncTimes = pepNEV.sync.timestamps/30; % double( tt(ii==0) )/30;

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

% %AZ 20081217: Get username from usertxt field (for saving in unit.source)
% UserName = get(findobj('Tag','usertxt'),'String');

TheseTimes = nevsorted{isingUnit}{2}/30; %double( tt(ii==ie) )/30;
SpikeTimes{isingUnit} = cell(p.nstim, nrepeats);
for istim = 1:p.nstim
    for irepeat = 1:nrepeats
        jj = (TheseTimes>TT1(istim,irepeat)) & (TheseTimes<=TT2(istim,irepeat));
        SpikeTimes{isingUnit}{istim,irepeat} = ((TheseTimes( jj )-TT1(istim,irepeat))./1000)'; % convert from milliseconds to seconds
    end
end

unit = UnitCreate;
unit.animal = p.animal;
unit.ichan = elec+1000; % electrode number
%unit.icell =  isingUnit;
unit.icell = gmm.icell(isingUnit);
unit.iseries = p.iseries;
unit.iexp =  p.iexp;
unit.stimdurs = stimdurs./1000;
unit.timestamp = datestr(now);
unit.spiketimes = SpikeTimes{isingUnit};
unit.datatype = 'spiketimes';
%AZ 20081217    unit.source = 'NEV file';
unit.source = UserName;  %AZ 20081217: save user name from usertxt field
unit.gmm = gmm;
unit.gmm.dprime = gmm.dprime(isingUnit); %AZ20090211: save dprime for selected icell only

unit.prototype = mwaves{isingUnit}; %AZ

UnitSave(unit,DIRS.spikes);

unitOUT(isingUnit) = unit;

end % end for isingUnit



end %END FUNCTION NEV2UNIT