function Spk = createSpikefile_fromold(DIRname, animal, iseries, iexp, igroup, iaddinfo, SpikeSorterType)
if nargin < 7
    SpikeSorterType = 'klusta';
end

if numel(iaddinfo) == 1 && numel(igroup) > 1
    iaddinfotemp = iaddinfo{1};
    iaddinfo = cell(1,numel(igroup));
    for group_idx = 1:numel(igroup)
        iaddinfo{group_idx} = iaddinfotemp;
    end
end

%for loading kwik file clusters
nCells = 0;
idx = 1;
for igroup_idx = 1:length(igroup)
    try
        switch SpikeSorterType
            case 'klusta'
                temp = getKwikSpikes(DIRname, animal, iseries, iexp, igroup(igroup_idx),iaddinfo{igroup_idx});
                tempNumCells = length(temp);
                for k = 1:tempNumCells
                    temp(k).ProbeIDs = iaddinfo{igroup_idx};
                end
                if idx==1
                    chans=temp;
                    numCells = tempNumCells;
                else
                    for icell = 1:tempNumCells
                        chans(numCells+1) = temp(icell);
                        numCells = numCells + 1;
                    end
                end
            case 'kilosort'
                temp = getKiloSortSpikes(DIRname, animal, iseries, iexp, igroup(igroup_idx),iaddinfo{igroup_idx});
                tempNumCells = length(temp);
                for k = 1:tempNumCells
                    temp(k).ProbeIDs = iaddinfo{igroup_idx};
                end
                if idx==1
                    chans=temp;
                    numCells = tempNumCells;
                else
                    for icell = 1:tempNumCells
                        chans(numCells+1) = temp(icell);
                        numCells = numCells + 1;
                    end
                end
        end
        idx = idx + 1;
    catch
        warning(['no spikes sorted on shank ' num2str(igroup(igroup_idx))]);
    end
    
    % getting the spike Times
    for ichan = 1:length(chans)
        Spk.animal{nCells + ichan} = animal;
        Spk.series{nCells + ichan} = iseries;
        Spk.exp{nCells + ichan} = iexp;
        Spk.sampleRate = chans(ichan).sampleRate;
        Spk.spikeIDs{nCells + ichan} = chans(ichan).id;
        Spk.chanIDs{nCells + ichan} = chans(ichan).ichan;
        Spk.bestchan{nCells + ichan} = chans(ichan).bestchan;
        Spk.waveform{nCells + ichan} = chans(ichan).waveform;
        Spk.ProbeIDs{nCells + ichan} = chans(ichan).ProbeIDs;
        Spk.spikeTimes{nCells + ichan} = chans(ichan).spiketimes;
    end
    
    nCells = length(chans);
end

save([DIRname filesep 'Spkdata_' animal '_' num2str(iseries) '_' num2str(iexp)], 'Spk', '-v7.3'); 
end

