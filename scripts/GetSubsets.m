function tidx = GetSubsets(S,contrast,gain,roomlength,outcome,speed_th,FnoBlanks,FnoAfterBlanks)
if nargin < 2 || isempty(contrast)
    contrast = 1:numel(S.contrastValues);
    gain = 1:numel(S.gainValues);
    roomlength = 1:numel(S.roomlengthValues);
    outcome = 1:numel(S.outcomeValues);
end
if nargin < 3 || isempty(gain)
    gain = 1:numel(S.gainValues);
    roomlength = 1:numel(S.roomlengthValues);
    outcome = 1:numel(S.outcomeValues);
end
if nargin < 4 || isempty(roomlength)
    roomlength = 1:numel(S.roomlengthValues);
    outcome = 1:numel(S.outcomeValues);
end
if nargin < 5 || isempty(outcome)
    outcome = 1:numel(S.outcomeValues);
end
if nargin < 6
    speed_th = S.SpeedThresh;
end
if nargin < 7
    FnoBlanks = true;
end
if nargin < 8
    FnoAfterBlanks = false;%true;%
end

tidx = false(size(S.smthBallSpd));
tidx(ismember(S.contrast,S.contrastValues(contrast)) & ismember(S.gain,S.gainValues(gain)) &...
     ismember(S.roomLength,S.roomlengthValues(roomlength)) & ismember(S.outcome,S.outcomeValues(outcome))) = true;

tidx = tidx & S.smthBallSpd > speed_th & ~isnan(S.smthBallSpd) & ~isnan(S.smthTrajSpd) & ~isnan(S.gain);

if FnoBlanks && isfield(S,'blanks')
    tidx = tidx & ~(S.blanks);
    for tt = -30:30
        tidx = tidx & circshift(~(S.blanks),tt);
    end
end
if FnoAfterBlanks && isfield(S,'afterblanks')
    tidx = tidx & ~(S.afterblanks);
end
end