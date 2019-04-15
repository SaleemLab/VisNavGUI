function TensorFromCerebus(animal, iseries, iexp, Fs)
% TensorFromCerebus
%
% TensorFromCerebus(animal, iseries, iexp) creates a tensor and puts the
% result in the global variable VDAQ
% 
% TensorFromCerebus(animal, iseries, iexp, Fs) lets you specify the sampling
% frequency Fs in Hz (DEFAULT: 200)
%
% Example:
%
% SetDefaultDirs;
% global VDAQ
% TensorFromCerebus('CATZ072',5,4)
%

% 2008-02 Matteo Carandini
% 2008-02 MC updated so channel numbers don't have to start at 1

global DIRS
global VDAQ

if nargin<4
    Fs = 200;
end

% make the 'unit' files
chans = UnitGetCerebus(animal, iseries, iexp);
ChanOffset = chans(1)-1;

p = ProtocolLoad(animal,iseries,iexp);

layout = UtahGetLayout(animal, iseries );

[nr,nc] = size(layout);

rr = cell(nr,nc);
ee = cell(nr,nc);

for ichan = chans % 1:96
    [ir,ic] = find( layout == (ichan-ChanOffset) ); 
    u = UnitLoad(DIRS.spikes,animal, iseries, iexp,ichan,999);
    [rr{ir,ic}, ee{ir,ic}] = UnitGetRates( u, 1/Fs ); 
end

LengthMatrix = NaN*zeros(nr,nc);
for ir = 1:nr
    for ic = 1:nc
        if ~isempty(rr{ir,ic})
            LengthMatrix(ir,ic) = numel(rr{ir,ic}{1});
        end
    end
end
if nanstd(LengthMatrix(:))>0
    error('Problem with durations');
end
nt = nanmedian(LengthMatrix(:));

%% create the VDAQ structure 

VDAQ = [];

VDAQ.tensor = cell(p.nstim,1);
for istim = 1:p.nstim
    VDAQ.tensor{istim} = zeros(nr,nc,nt);
    for ir = 1:nr
        for ic = 1:nc
            if ~isempty( rr{ir,ic} )
                VDAQ.tensor{istim}(ir,ic,1:length(rr{ir,ic}{istim})) = rr{ir,ic}{istim};
            end
        end
    end
end

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
VDAQ.nframes    = nt;
VDAQ.nrepeats   = u.nrepeats;
VDAQ.durs       = u.stimdurs;
VDAQ.FrameRate  = VDAQ.nframes/VDAQ.durs(1);

% fields that we probably don't need to populate:
% VDAQ.fileList = ;
% VDAQ.MeanIntensities = ;
% VDAQ.header = ;
% VDAQ.nsummedframes =;



