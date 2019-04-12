function obj = Calculate1Dmaps(obj, varXname, Tsmthwin, Xbinsize, Xsmthwin, delayT, Fcircular)
size_th = 100;
if nargin < 5
    Xsmthwin = 2*Xbinsize;
end
if nargin < 6
    delayT = 0;
end
if nargin < 7
    Fcircular = true;
end
disp(['time window = ' num2str(Tsmthwin) ' ms']);

nbcont = numel(obj.SubsetVal.contrast) + 1;
nbgain = numel(obj.SubsetVal.gain);
nbroomlength = numel(obj.SubsetVal.roomlength);
nboutcome = numel(obj.SubsetVal.outcome);

spiketrain = circshift(obj.data.es.spikeTrain,[-delayT 0]);
sampleRate = mean(1./obj.data.es.sampleSize);
for icell = 1:size(spiketrain,2)
    spiketrain(:,icell) = smthInTime(spiketrain(:,icell), sampleRate, Tsmthwin, 'same', [], 'boxcarsum_centered');
end

if strcmp(varXname,'trajPercent') || strcmp(varXname,'spiketrajPercent')
    minvarXname = 0;
    maxvarXname = 100;
    numBins = round(maxvarXname/Xbinsize);
elseif strcmp(varXname,'trajPercentunwrapped')
    minvarXname = 0;
    maxvarXname = 200;
    numBins = round(maxvarXname/Xbinsize);
elseif strcmp(varXname,'LFPphase') || strcmp(varXname,'LFPphase2') || strcmp(varXname,'spikePhase')
    minvarXname = 0;
    maxvarXname = 360;
    numBins = round(maxvarXname/Xbinsize);
elseif strcmp(varXname,'gain')
    [gainval,~,ic] = unique(obj.data.es.(varXname)(~isnan(obj.data.es.(varXname))),'sorted');
    minvarXname = min(ic);
    maxvarXname = max(ic);
    numBins = numel(unique(ic));
else
    minvarXname = min(obj.data.es.(varXname)(:));
    maxvarXname = max(obj.data.es.(varXname)(:));
end
if Fcircular
    for k = 1:size(obj.data.es.(varXname),2)
        Xtemp = unwrap(obj.data.es.(varXname)(:,k)/maxvarXname*2*pi)*maxvarXname/(2*pi);
        Xtemp(~isnan(Xtemp)) = smthInTime(Xtemp(~isnan(Xtemp)), sampleRate, Tsmthwin, 'same', [], 'boxcar_centered');
        varX(:,k) = mod(Xtemp,maxvarXname);
    end
else
    for k = 1:size(obj.data.es.(varXname),2)
        Xtemp = obj.data.es.(varXname)(:,k);
        Xtemp(~isnan(Xtemp)) = smthInTime(Xtemp(~isnan(Xtemp)), sampleRate, Tsmthwin, 'same', [], 'boxcar_centered');
        varX(:,k) = Xtemp;
    end
end

if strcmp(varXname,'gain')
    for k = 1:size(obj.data.es.(varXname),2)
        [gainval,~,ic] = unique(varX(~isnan(varX(:,k)),k),'sorted');
        varX(~isnan(varX(:,k)),k) = ic;
    end
end

for k = 1:size(obj.data.es.(varXname),2)
    [varX(:,k), bins] = normalise1var(varX(:,k), numBins, [], [minvarXname maxvarXname]);
end

if size(varX,2) == 2
    varX_probe1 = repmat(varX(:,1),[1 sum(obj.CellInfo.Probe == 1)]);
    varX_probe2 = repmat(varX(:,2),[1 sum(obj.CellInfo.Probe == 2)]);
    varX = cat(2,varX_probe1,varX_probe2);
end

win = max(1,round(Tsmthwin/(1000/sampleRate)));
T = obj.data.es.sampleSize.*win;

gbase = find(obj.SubsetVal.gain == mode(obj.data.es.gain));
gainlist = 1:nbgain;
gainlist = [gbase gainlist(gainlist~=gbase)];
nShf = 100;%500;
obj.maps1d.(varXname) = [];
for c = 1:nbcont
    for g = gainlist
        for r = 1:nbroomlength
            for o = 1:nboutcome
                obj.maps1d.(varXname){c, g, r, o} = ToneDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize));
                obj.maps1d.(varXname){c, g, r, o}.sampleRate = obj.data.es.sampleRate;
                obj.maps1d.(varXname){c, g, r, o}.Fcircular = Fcircular;
                obj.maps1d.(varXname){c, g, r, o}.numBins = numBins;
                obj.maps1d.(varXname){c, g, r, o}.bins = bins;
                obj.maps1d.(varXname){c, g, r, o}.qthreshold = 1;
                obj.maps1d.(varXname){c, g, r, o}.Fdiscarditer = false;
                obj.maps1d.(varXname){c, g, r, o}.Fgoodcells = obj.CellInfo.Goodcluster;
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                if strcmp(varXname,'trajPercentunwrapped')
                    idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                elseif strcmp(varXname,'gain')
                    idx = false(size(obj.data.es.(varXname)));
                    for gg = gainlist
                        idx = idx | obj.getSubsets(contidx, gg, r, o, obj.SpeedThresh, true, true);
                    end
                else
                    idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                end
                if o == 3 && c == nbcont
                    if sum(idx) > size_th
                        if ~strcmp(varXname,'gain')
                            obj.maps1d.(varXname){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                        else
                            if g == gbase
                                obj.maps1d.(varXname){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                            end
                        end
                    end
                end
                
                obj.maps1d.([varXname '_Shf']){c, g, r, o} = ToneDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize));
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.sampleRate = obj.data.es.sampleRate;
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.Fcircular = Fcircular;
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.numBins = numBins;
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.bins = bins;
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.qthreshold = 1;
                if ~strcmp(varXname,'gain')
                    obj.maps1d.([varXname '_Shf']){c, g, r, o}.Fdiscarditer = true;
                else
                    obj.maps1d.([varXname '_Shf']){c, g, r, o}.Fdiscarditer = false;
                end
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.Fgoodcells = obj.CellInfo.Goodcluster;
                if ~strcmp(varXname,'trajPercentunwrapped')
                    obj.maps1d.([varXname '_Shf']){c, g, r, o}.Fshuffle = true;
                else
                    obj.maps1d.([varXname '_Shf']){c, g, r, o}.Fshufflehalf = true;
                end
                obj.maps1d.([varXname '_Shf']){c, g, r, o}.kfold = nShf;
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                if strcmp(varXname,'trajPercentunwrapped')
                    idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                elseif strcmp(varXname,'gain')
                    idx = false(size(obj.data.es.(varXname)));
                    for gg = gainlist
                        idx = idx | obj.getSubsets(contidx, gg, r, o, obj.SpeedThresh, true, true);
                    end
                else
                    idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                end
                if o == 3 && c == nbcont
                    if sum(idx) > size_th
                        if ~strcmp(varXname,'gain')
                            obj.maps1d.([varXname '_Shf']){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                        else
                            if g == gbase
                                obj.maps1d.([varXname '_Shf']){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                            end
                        end
                    end
                end
                
                if ~strcmp(varXname,'gain')
                    obj.maps1d.([varXname '_2fold']){c, g, r, o} = ToneDimMap('Xsmth_win', floor(Xsmthwin/Xbinsize));
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.sampleRate = obj.data.es.sampleRate;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.Fcircular = Fcircular;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.numBins = numBins;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.bins = bins;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.qthreshold = 1;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.Fgoodcells = obj.CellInfo.Goodcluster;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.kfold = 2;
                    obj.maps1d.([varXname '_2fold']){c, g, r, o}.Fdiscarditer = false;
                    if c > numel(obj.SubsetVal.contrast)
                        contidx = find(obj.SubsetVal.contrast>0);
                    else
                        contidx = c;
                    end
                    if strcmp(varXname,'trajPercentunwrapped')
                        idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                    elseif strcmp(varXname,'gain')
                        idx = false(size(obj.data.es.(varXname)));
                        for gg = gainlist
                            idx = idx | obj.getSubsets(contidx, gg, r, o, obj.SpeedThresh, true, true);
                        end
                    else
                        idx = obj.getSubsets(contidx, g, r, o, obj.SpeedThresh, true, true);
                    end
                    if o == 3 && c == nbcont
                        if sum(idx) > size_th
                            obj.maps1d.([varXname '_2fold']){c, g, r, o}.trainSpikeMap(varX(idx,:), spiketrain(idx,:), T(idx));
                        end
                    end
                end
            end
        end
    end
end
end