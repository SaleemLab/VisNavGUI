function Unit998FromNEV(animal,iseries,iexp)
% Unit998FromNEV
%
% Unit998FromNEV works on animal, iseries, and iexp that are currently in
% PICK.
%
% Unit998FromNEV(animal,iseries,iexp) lets you specify animal, iseries, and
% iexp.
%
% EXAMPLE: 
% SetDefaultDirs
% Unit998FromNEV('M110320TS',4,7)
%
% 2011-05 MC made it from earlier code
% 2011-05 ND made it work with saving by experiment or saving by repeat

global PICK pepNEV DIRS

if nargin < 3
    animal  = PICK.animal;
    iseries = PICK.iseries;
    iexp    = PICK.iexp;
end

p = ProtocolLoad(animal, iseries, iexp);

%% get layout (much of this code should be encapsulated in layout object)

[arrayLayout, ~, nshanks] = MichiganGetLayout(animal, iseries);
nchans = length(arrayLayout);
nChansPerShank = floor(nchans/nshanks); % an effort to drop the extra channels
nchans = nChansPerShank*nshanks;
ChanTags = zeros(nChansPerShank,nshanks);
if nshanks==1
    ChanTags(:) = 1:nchans;
else
    for ishank = 1:nshanks
        ChanTags(:,ishank) = 100*ishank + (1:nChansPerShank);
    end
end

% unscramble them
NewChanTags = ChanTags;
for ichan = 1:nchans
    NewChanTags(ichan) = ChanTags(arrayLayout==ichan);
end
ChanTags = NewChanTags;

%%

% this works for data saved repeat by repeat
NevDir = fullfile(DIRS.Cerebus, animal, num2str(iseries), num2str(iexp));
% NevFiles = dir( fullfile(NevDir, sprintf('%s_%d_%d_*.nev',animal,iseries,iexp) ) );
NevFiles = dir( fullfile(NevDir, sprintf('%s_%d_%d*.nev',animal,iseries,iexp) ) );
ExptDate = NevFiles(1).date;
nNevFiles = length(NevFiles);
% if nNevFiles == 1 and p.nrepeats > 1, then saved by experiment
if nNevFiles == 1
    if p.nrepeats > 1
        byRptFlag  = 0;
        byExptFlag = 1;
        nrepeats = p.nrepeats;
    else % p.nrepeats ==1
        % don't really know ==> could be saving by repeat or by experiment
        byRptFlag  = 0;
        byExptFlag = 0;
    end
else % many NevFiles ==> saving by repeat
    byRptFlag  = 1;
    byExptFlag = 0;
    nrepeats = min(nNevFiles,p.nrepeats); % this only is good if data has been saved by repeat
end


%% copy files to local directory so things are faster

if byRptFlag
    fprintf('Copying NEV files to local directory');
    TempDir = 'C:\WINDOWS\Temp\';
    for irepeat = 1:nrepeats
        fprintf('.');
        NevFileName = sprintf('%s_%d_%d_%d.nev',animal,iseries,iexp,irepeat);
        copyfile(fullfile(NevDir,NevFileName), fullfile(TempDir,NevFileName),'f');
    end
    fprintf(' done\n');
elseif byExptFlag
    fprintf('Copying NEV files to local directory');
    TempDir = 'C:\WINDOWS\Temp\';
    fprintf('.');
    NevFileName = sprintf('%s_%d_%d.nev',animal,iseries,iexp);
    copyfile(fullfile(NevDir,NevFileName), fullfile(TempDir,NevFileName),'f');
    fprintf(' done\n');
end


%% load the data

SpikeTimes = cell(p.nstim, nrepeats, nchans);
stimdurs = zeros(p.nstim,nrepeats);
% TotNSpikes = zeros(nchans,1);
% WaveForms = zeros(48,100,nchans); % this assumes the waveforms are 48 samples long...

if byRptFlag
    for irepeat = 1:nrepeats
        
        NevFileName = sprintf('%s_%d_%d_%d.nev',animal,iseries,iexp,irepeat);
        nevopen_outcome = nevopen(fullfile(TempDir,NevFileName)); % NevDir
        
        if ~nevopen_outcome, error('create998:FileOpenError','Could not open file!'); end
        
        fprintf('Loading spikes for repeat %d/%d -- ',irepeat,nrepeats);
        
        SyncTimes = pepNEV.sync.timestamps/30;
        
        TT1 = zeros(p.nstim,1);
        TT2 = zeros(p.nstim,1);
        for istim = 1:p.nstim
            ipres = p.seqnums(istim,irepeat)-p.nstim*(irepeat-1);
            TT1(istim,1) = SyncTimes(2*ipres-1);
            TT2(istim,1) = SyncTimes(2*ipres+0);
        end
        
        stimdurs(:,irepeat) = (TT2 - TT1); % in milliseconds    for ichan = 1:nchans
        if any(stimdurs(:,irepeat)==0), break; end
        
        for ichan = 1:nchans
            
            tt = nevwaves(ichan);             % was [tt, vv] = nevwaves(ichan);
            
            nspikes = length(tt);
            fprintf('%d ; ',nspikes);
            
            if nspikes > 0
                
                % AverageWaveForm(:,ichan) = AverageWaveForm(:,ichan) + sum(waveforms,2);
                % TotNSpikes(ichan) =  TotNSpikes(ichan) + nspikes;
                
                TheseTimes = tt/30;
                
                for istim = 1:p.nstim
                    jj = (TheseTimes>TT1(istim)) & (TheseTimes<=TT2(istim));
                    SpikeTimes{istim,irepeat,ichan} = ((TheseTimes( jj )-TT1(istim))./1000)'; % convert from milliseconds to seconds
                end
            end
        end
        fprintf('\n');
        
        nevclose;
    end
elseif byExptFlag
    for irepeat = 1:nrepeats
        
        NevFileName = sprintf('%s_%d_%d.nev',animal,iseries,iexp);
        nevopen_outcome = nevopen(fullfile(TempDir,NevFileName)); % NevDir
        
        if ~nevopen_outcome, error('create998:FileOpenError','Could not open file!'); end
        
        fprintf('Loading spikes for repeat %d/%d -- ',irepeat,nrepeats);
        
        SyncTimes = pepNEV.sync.timestamps/30;
        
        TT1 = zeros(p.nstim,1);
        TT2 = zeros(p.nstim,1);
        for istim = 1:p.nstim
            ipres = p.seqnums(istim,irepeat);
            TT1(istim,1) = SyncTimes(2*ipres-1);
            TT2(istim,1) = SyncTimes(2*ipres+0);
        end
        
        stimdurs(:,irepeat) = (TT2 - TT1); % in milliseconds    for ichan = 1:nchans
        if any(stimdurs(:,irepeat)==0), break; end
        
        for ichan = 1:nchans
            
            tt = nevwaves(ichan);             % was [tt, vv] = nevwaves(ichan);
            
            nspikes = length(tt);
            fprintf('%d ; ',nspikes);
            
            if nspikes > 0
                
                % AverageWaveForm(:,ichan) = AverageWaveForm(:,ichan) + sum(waveforms,2);
                % TotNSpikes(ichan) =  TotNSpikes(ichan) + nspikes;
                
                TheseTimes = tt/30;
                
                for istim = 1:p.nstim
                    jj = (TheseTimes>TT1(istim)) & (TheseTimes<=TT2(istim));
                    SpikeTimes{istim,irepeat,ichan} = ((TheseTimes( jj )-TT1(istim))./1000)'; % convert from milliseconds to seconds
                end
            end
        end
        fprintf('\n');
        
        nevclose;
    end
end

%% Save the units

for ichan = 1:nchans
    % UnitFileName = UnitGetFilename( animal, iseries, iexp, ChanTags(ichan), 998 );
    % if exist(UnitFileName,'file'), error('create998:FUnitExistsError','Unit 998 already exists!'); end

    unit = UnitCreate;
    unit.animal     = p.animal;
    unit.iseries    = p.iseries;
    unit.iexp       = p.iexp;
    unit.ichan      = ChanTags(ichan); % electrode number
    unit.icell      = 998;
    unit.stimdurs   = stimdurs./1000;
    unit.timestamp  = ExptDate;
    unit.spiketimes = cell(p.nstim,nrepeats);
    for istim = 1:p.nstim
        for irepeat = 1:nrepeats
            unit.spiketimes{istim, irepeat} = SpikeTimes{istim, irepeat, ichan};
        end
    end
    
    unit.datatype = 'spiketimes';
    unit.source = 'Original NEV file';

    % unit.prototype = AverageWaveForm(:,ichan)/ TotNSpikes(ichan);

    UnitSave(unit,DIRS.spikes);

end % for ichan

fprintf('**************** ¡¡¡DONE!!! ****************\n');

