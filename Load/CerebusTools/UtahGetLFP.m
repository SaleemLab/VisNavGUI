function [] = UtahGetLFP(animal, iseries, iexp, RepeatFlag, nRepeats, prestimTime)
% UtahGetLFP converts the LFPs from the Utah array to a VDAQ.tensor averaged across
% repeats
%
% UtahGetLFP
% loads the data files for the experiment specified by the global PICK and
% updates the global variable VDAQ
%
% UtahGetLFP(ExptTag)
% loads the data specified by ExptTag.animal, ExptTag.iseries, ExptTag.iexp
%
% UtahGetLFP(animal, iseries, iexp)
% lets you specify animal, iseries, iexp (DEFAULT: the values in the global
% variable PICK).
%
% UtahGetLFP(animal, iseries, iexp, 'keeprepeats')
% lets you keep the individual repeats
% DEFAULT: average across repeats.
%
% UtahGetLFP(animal, iseries, iexp, '', nRepeats)
% lets you specify the first n repeats to take
% DEFAULT: use all repeats.
%
% UtahGetLFP(animal, iseries, iexp, '', nRepeats, prestimTime)
% allows you to specify a time window in ms before and after the stimulus duration 
% DEFAULT: 0
% if nRepeats is set to 0, all repeats are used
%
%
% to make sense of this, look at toolbox BrainPics
%
% 2008-04 SK and LB created it
% 2008-05 LB converted alldata to single (to save memory)
% 2008-08 LB added the option to limit the # of repeats
% 2009-11 MS added the option to include LFP before and after stimulus duration

global pepNEV
global DIRS
global PICK
global VDAQ

if nargin < 6, prestimTime = 0; end
if nargin < 5, nRepeats = 0; end;
if nargin < 4, RepeatFlag = ''; end
if nargin < 2
    ExptTag = animal; % the first argument was ExptInfo
    animal  = ExptTag.animal;
    iseries = ExptTag.iseries;
    iexp 	= ExptTag.iexp;
end
if nargin < 1
    if isempty(PICK)
        error('<UtahGetLFP> Could not find global variable PICK');
    end
    animal = PICK.animal;
    iseries = PICK.iseries;
    iexp = PICK.iexp;
end

%% load some information
p = ProtocolLoad(animal, iseries, iexp);
layout = UtahGetLayout(animal, iseries );
[nr,nc] = size(layout);

if ~nRepeats
    nRepeats = p.nrepeats;
end

fileName = fullfile(DIRS.Cerebus, animal, sprintf('u%03d_%03d', iseries, iexp));
nevopen([fileName '.nev']);

% sync times in units of clock ticks
stimOnTimes = pepNEV.sync.timestamps(1:2:end);
stimOffTimes = pepNEV.sync.timestamps(2:2:end);

% convert prestimTime to number of nev frames
if prestimTime
    nframes = stimOffTimes(1)-stimOnTimes(1);
    frameRate = nframes / p.pfiledurs(1);
    frameRate = round(frameRate/1000)*1000; % to get to a round number
    prestimTime = frameRate * prestimTime/1000; % prestimTime is in ms
end

% error checking
if diff([length(stimOnTimes) length(stimOffTimes) max(p.seqnums(:))]) ~= 0
    error('umatched number of start/stop sync times or number of trials for %s', fileName);
end

%% open the file containing the lfp data
fid = fopen([fileName '.ns3'], 'r', 'ieee-le');

% extract some information
fseek(fid, 24, 'bof');
infoVec = fread(fid,8,'char');
period = infoVec(1);
nchan =  infoVec(5);
enum = fread(fid, nchan, 'int32');  % electrode numbers: 384 (96*4)
filePtrPosition = ftell(fid); 

%% first run: find minimum number of samples
iTrials = 0;
fprintf('\nchecking number of samples ');
for istim = 1 : p.nstim
    fprintf('.');
    for irepeat =  1 : p.nrepeats
        iTrials = iTrials + 1;
        % find the appropriate segment
        t0 = stimOnTimes(p.seqnums(istim, irepeat)) - prestimTime;
        t1 = stimOffTimes(p.seqnums(istim, irepeat)) + prestimTime;
        % convert to units of samples
        T0 = ceil(t0/period);
        T1 = floor(t1/period);
        % position the file pointer
        fseek(fid,T0*nchan*2,'cof');
        r = fread(fid,nchan*(T1-T0+1),'int16=>int16');
        r = reshape(r,[nchan,T1-T0+1]);
        nSamples(iTrials) = size(r, 2); %#ok<AGROW>
        % reset the file pointer
        fseek(fid, filePtrPosition, 'bof');
    end
end
fprintf('\t\t\t done.\n');

nSamples = unique(nSamples);
minSamples = min(nSamples);

if size(nSamples,1) > 1 || size(nSamples,2) > 1
    fprintf('\nwarning: mismatch in number of samples between stimuli\n');
    fprintf('%d ', nSamples);
    fprintf('\n\n');
end 

%% get VDAQ tensor
VDAQ.tensor = {};
for istim = 1:p.nstim
    VDAQ.tensor{istim} = zeros(nr, nc, minSamples);
end

m = zeros(nr,nc, minSamples);

fprintf('Converting LFP to Tensor');
for istim = 1 : p.nstim
    fprintf('.');
    for irepeat =  1 : nRepeats
        % find the appropriate segment
        t0 = stimOnTimes(p.seqnums(istim, irepeat)) - prestimTime;
        t1 = stimOffTimes(p.seqnums(istim, irepeat)) + prestimTime;
        % convert to units of samples
        T0 = ceil(t0/period);
        T1 = floor(t1/period);
        % position the file pointer
        fseek(fid,T0*nchan*2,'cof');
        % read the data, convert to floats to normalize them later
        r = fread(fid,nchan*(T1-T0+1),'int16=>float');
        r = reshape(r,[nchan,T1-T0+1]);
        for ie = 1:length(enum)
            mye = enum(ie);
            [ir,ic] = find( layout == mye);
            if isempty(ir) || isempty(ic)
                % fprintf('\n\t<UtahGetLFP> Warning: electrode number %d not found in layout', mye);
                continue;
            end
            m(ir,ic,:) = r(mye,1:minSamples);
        end; clear r; 
        VDAQ.tensor{istim} = VDAQ.tensor{istim} + m/nRepeats; % divide by the # of repeats
        if strcmp(RepeatFlag,'keeprepeats')
            VDAQ.alldata{istim, irepeat} = single(m);
        end
        
        % repeats{irepeat} = r(:, 1:minSamples);
        % reset the file pointer
        fseek(fid, filePtrPosition, 'bof');
    end
end
fclose(fid);
fprintf('\t\t\t done.\n');

%% build the rest of the VDAQ structure
VDAQ.animal     = animal;
VDAQ.iseries    = iseries;
VDAQ.iexp       = iexp;
VDAQ.Frame0List = [];
VDAQ.ResizeFactor = 1;
VDAQ.nstim      = p.nstim;
VDAQ.nx         = nr;
VDAQ.ny         = nc;
VDAQ.transposed = 0;
VDAQ.binX       = 1;
VDAQ.binY       = 1;
VDAQ.MmPerCameraPix = 0.4; % assume 400 um spacing
VDAQ.nframes    = minSamples;
VDAQ.nrepeats   = nRepeats;
VDAQ.durs       = p.pfiledurs+2*prestimTime/1000; % prestimTime is in ms
if ~isempty(VDAQ.durs)
    VDAQ.FrameRate  = VDAQ.nframes/VDAQ.durs(1);
end

% fields that we probably don't need to populate:
% VDAQ.fileList = ;
% VDAQ.MeanIntensities = ;
% VDAQ.header = ;
% VDAQ.nsummedframes =;
