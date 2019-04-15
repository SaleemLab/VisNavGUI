function chanlist = UnitGetCerebus(animal, iseries, iexp, strOverwrite, ChunkLength)
% UnitGetCerebus takes Cerebus multiunit data and turns it into Unit data
%
% chanlist = UnitGetCerebus(animal, iseries, iexp) loads Cerebus multiunit data from
% NEV files and saves them into 'unit' data. 
%
% Channel names start at 1001, and unit numbers are 999 (the standard for
% multiunit activity).
%
% This function relies on a global DIRS to know where the relevant
% directories are. Do SetDefaultDirs to set this global to default values. 
%
% UnutGetCerebus with no arguments uses the global PICK (i.e. it assumes
% you are running ExptPick and choosing experiments with it).
%
% UnitGetCerebus(animal, iseries, iexp) lets you specify manually.
%
% UnitGetCerebus(animal, iseries, iexp, 'overwrite') forces an overwrite of
% existing unit files (DEFAULT: '', which means do not overwrite).
%
% UnitGetCerebus(animal, iseries, iexp, '', ChunkLength) lets you specify
% the length of data chunks to read from the NEV file (DEFAULT: 1000000).
%
% During an experiment, it is likely that you will want to do:
% DIRS.Cerebus = '\\zanna\cerebus_data'. 
%
%
%   Example:  
% 
% SetDefaultDirs;
% UnitGetCerebus('CATZ072',5,4])

% 2007 Ian Nauhaus created as 'SendToZAnalyze'
% 2008-02 LB added CerebusDir as input argument
% 2008-02 SK added target directory as input argument, renamed handling of directories
% 2008-02 LB streamlined # of arguments and documentation
% 2008-02 MC reorganized and renamed as 'UnitGetCerebus'
% 2008-02 MC added 'sites' output and renamed the channels from 1001 on.
% 2008-02 MC major changes: reads NEV files, saves all channels
% 2008-02 MC checks first to see if the files exist 
% 2008-02 MC added 'source' field to unit data
% 2008-02 MC major change: introduced memory mapping of the NEV file
% 2008-03 LB bug fix: stimdurs and spiketimes are expressed in sec (was millisec before)
% 2008-03 MC makes it call UnitCreate
% 2008-03 LB bug fix (l. 172): replaced unit.spiketimes = {SpikeTimes{ie}} by unit.spiketimes = SpikeTimes{ie};
% 2008-06 LB bug fix: unit.spiketimes must have format [1, nspikes] to be compatible with existing functions (was [nspikes, 1])
% 2010-03 MC added support for PICK

global DIRS PICK

if nargin<1
    if isempty(PICK), error('You have to define PICK and make it global (or run ExptPick)'); end
    animal = PICK.animal;
    iseries = PICK.iseries;
    iexp = PICK.iexp;
end

if nargin<4
    strOverwrite = '';
end

if nargin<5
    ChunkLength = 1000000; % one million spikes at a time
end

p = ProtocolLoad(animal,iseries,iexp);

%% See if the files already exist

if ~strcmp( strOverwrite, 'overwrite' )
    filename = UnitGetFilename( animal,iseries,iexp, 1001, 999 );
    if exist(fullfile(DIRS.spikes,filename),'file')
        chanlist = 1001:1096;
        for ichan = chanlist
            filename = UnitGetFilename( animal,iseries,iexp, ichan, 999 );
            if ~exist(fullfile(DIRS.spikes,filename),'file')
                chanlist(ichan) = NaN;
            end
        end
        chanlist(isnan(chanlist)) = [];
        fprintf('\nFile %s already exists. \nRun UnitGetCerebus with option ''overwrite'' if you wish to overwrite.\n',filename);
        return
    end
end

%% Read the NEV file

chanlist = []; % start pessimistic
    
filename = fullfile( DIRS.Cerebus,animal,sprintf('u%03d_%03d.nev',iseries,iexp) );

if ~exist(filename,'file')
    fprintf('Did not find file %s.\n',filename);
    return;
end;

fid = fopen(filename, 'r', 'ieee-le');
if fid == -1, 
    fprintf('could not open file %s.\n',filename);
    return;
end;

fseek(fid, 12, 'bof');
offset = fread(fid, 1, 'uint32');
packetlen = fread(fid, 1, 'uint32');

fseek(fid,0,'eof');
totnspikes = floor( (ftell(fid)-offset)/packetlen );

fclose(fid);

fprintf('NEV file has %2.1f million spikes.\n',totnspikes/(10^6));

ii = zeros(totnspikes,1,'uint8');
tt = zeros(totnspikes,1,'uint32');

nloops = ceil(totnspikes/ChunkLength); 

fprintf('Will load them in %d chunks of %2.1f million spikes\n',nloops,ChunkLength/(10^6));
for iloop = 1:nloops
    LoopOffset = ChunkLength*(iloop-1);
    nRemaining = totnspikes - LoopOffset;
    nReadsThisLoop = min( nRemaining, ChunkLength );
    fprintf('Chunk %d of %d',iloop,nloops);
    tic;
    m = memmapfile(filename,            ...
        'Format', {   'uint32' [1 1] 'TimeStamp'; ...
        'uint16' [1 1] 'ChanName'; ...
        'int16' [48 1] 'Waveform'; ...
        'int16' [1 1]  'Padding'},  ...
        'Offset', offset+LoopOffset*packetlen, 'Repeat', nReadsThisLoop);
    dd = m.Data;
    ii(LoopOffset+(1:nReadsThisLoop)) = [dd.ChanName];
    tt(LoopOffset+(1:nReadsThisLoop)) = [dd.TimeStamp];
    t = toc;
    fprintf(' (took %2.1f s: should be done in %2.1f min)\n',t,t*(nloops-iloop)/60);    
end

%% Reorganize sync times

SyncTimes = double( tt(ii==0) )/30;

nrepeats = floor( length(SyncTimes)/2 / p.nstim );

if nrepeats > p.nrepeats
    fprintf('Error! protocol file has information only for %d of %d repeats. Have to stop here\n',p.nrepeats, nrepeats);
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

% Andrea stop here when you are looking for durations...

%% Rearrange spike times and save unit structures

ne = double(max(ii));

SpikeTimes = cell(ne,1); 
for ie = 1:ne
    fprintf('Electrode %d of %d -- ', ie, ne);
    TheseTimes = double( tt(ii==ie) )/30;
    SpikeTimes{ie} = cell(p.nstim, nrepeats);
    for istim = 1:p.nstim
        for irepeat = 1:nrepeats
            jj = (TheseTimes>TT1(istim,irepeat)) & (TheseTimes<=TT2(istim,irepeat));
            SpikeTimes{ie}{istim,irepeat} = ((TheseTimes( jj )-TT1(istim,irepeat))./1000)'; % convert from milliseconds to seconds
        end
    end
    drawnow
    unit = UnitCreate;
    unit.animal = animal;
    unit.ichan = ie+1000;
    unit.icell =  999;
    unit.iseries = iseries;
    unit.iexp =  iexp;
    unit.stimdurs = stimdurs./1000;
    unit.timestamp = datestr(now); 
    unit.spiketimes = SpikeTimes{ie};
    unit.datatype = 'spiketimes';
    unit.source = 'NEV file'; 
    
    UnitSave(unit,DIRS.spikes);
end

chanlist = (1:ne)+1000;

return


% % % sanity checks 
% 
% SpikeTimes = cell(ne,1); % here it would make sense to go into stimuli and repeats
% for ie = 1:ne
%     SpikeTimes{ie} = tt(ii==ie)/30;
% end
% 
% % nevopen(fullfile( DIRS.Cerebus,animal,sprintf('u%03d_%03d.nev',iseries,iexp)));
% % ie = 10;
% % [tt2, ww2] = nevwaves(ie);
% tt2 = tt2/30;
% 
% length(tt2)
% length(SpikeTimes{ie})
% 
% figure; plot( tt2, SpikeTimes{ie} )
% 
% e_ii = find(ii==ie);
% figure; 
% subplot(2,1,1);
% plot( ww(e_ii(1),:) );
% subplot(2,1,2);
% plot(ww2(:,1))



% 
% % nevroot = [CerebusDir '\' animal '\u' sprintf('%03d',iseries) '_'  sprintf('%03d',iexp)];
% nevroot = fullfile( CerebusDir,animal,sprintf('u%03d_%03d',iseries,iexp) );
% nev = load([nevroot '.mat'],'resp'); 
% 
% % should load the nev files, not the mat files!!!
% 
% switch MUAflag
%     case 1
%         idchan = 3;
%     case 2
%         idchan = 4;
%     case 3
%         idchan = 5;
%     otherwise
%         idchan = [];
% end
% 
% [Nstims Nreps] = size(nev.resp);
% 
% Stimes = cell(length(chans),1);
% stimdurs = NaN*zeros(Nstims,Nreps);
% 
% for stim = 1:Nstims
%     for rep = 1:Nreps
% 
%         if length(nev.resp{stim,rep}{151,2}) < 2
%             fprintf('<SendtoZAnalyze> Missing start/stop time for stim %d, repeat %d\n', stim, rep);
%             continue;
%         end
% 
%         Tstart = nev.resp{stim,rep}{151,2}(1);  %Trial start in sec.
%         Tstop = nev.resp{stim,rep}{151,2}(2);   %Trial stop in sec.
% 
%         for ichan = 1:length(chans)
%             chan = chans(ichan);
%             if MUAflag == -1
%                 Stimes{ichan}{stim,rep} = [...
%                     nev.resp{stim,rep}{chan,2}; ...
%                     nev.resp{stim,rep}{chan,3}; ...
%                     nev.resp{stim,rep}{chan,4}; ...
%                     nev.resp{stim,rep}{chan,5}];
%                 Stimes{ichan}{stim,rep} = sort(Stimes{ichan}{stim,rep}) - Tstart;
%             else
%                 Stimes{ichan}{stim,rep} = nev.resp{stim,rep}{chan,idchan} - Tstart;
%             end
%             if ichan == 1
%                 stimdurs(stim,rep) = Tstop - Tstart;
%             end
%         end
% 
%     end
% end
% 
% %% Make unit structures and save them
% 
% for ichan = 1:length(chans)
%     chan = chans(ichan);
%     unit = struct('animal', animal, 'ichan', chan+1000, 'icell', 999,'iseries', iseries,'iexp', iexp,...
%         'stimdurs',stimdurs,'timestamp','','prototype',[],'neighborhood',[],'spiketimes',{Stimes{ichan}},...
%         'datatype','spiketimes'); %#ok<NASGU>
%     UnitSave(unit,SpikeDir);
% end
% 
% chanlist = chans+1000;

