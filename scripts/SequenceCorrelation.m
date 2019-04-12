function SequenceCorrelation
nShf = 100;
repPval = sum(repmat(EXP.CellInfo.fieldAmp{size(EXP.CellInfo.fieldAmp,1), 2, 1, 3},[nShf 1])' >= EXP.CellInfo.fieldShfAmpiter{size(EXP.CellInfo.fieldShfAmpiter,1), 2, 1, 3},2)/nShf;
CA1cellidx = find(EXP.CellInfo.Goodcluster & ~EXP.CellInfo.Finterneuron & EXP.CellInfo.Probe==1 & repPval' < 0.01 & nanmean(EXP.data.es.spikeTrain,1)*60<4);
CA1spikes = EXP.data.es.spikeTrain(:,CA1cellidx);
CA1fieldPos = EXP.CellInfo.fieldPos{3,2,1,3}(CA1cellidx);
[~,CA1sorted] = sort(CA1fieldPos,'ascend');
CA1fieldPos = CA1fieldPos(CA1sorted);
CA1spikes = CA1spikes(:,CA1sorted);
rho = NaN(1,size(CA1spikes,1));
rho_Shf = NaN(100,size(CA1spikes,1));
for tt = 5:10000%(size(CA1spikes,1)-round(0.3/(1/60)))
    if sum(sum(CA1spikes((tt-4):tt-1,:),1)>0,2) < inf
        evttemp = CA1spikes(tt:(tt+round(0.5/(1/60))),:);
        spkSeq = [];
        for kk = 1:size(evttemp,1)
            spkID = find(evttemp(kk,:)>0);
            evttemp(:,spkID) = 0;
            spkSeq = [spkSeq spkID(:)'];
        end
        if numel(spkSeq) > 3%size(CA1spikes,2)/3
            for kcell = 0:numel(spkSeq)-1
                rankCorr = corr(spkSeq',circshift(sort(unique(spkSeq'),'ascend'),kcell),'type','Spearman');
                if abs(rankCorr) > abs(rho(tt)) || isnan(rho(tt))
                    rho(tt) = rankCorr;
                end
                for ishf = 1:100
                    rankCorr = corr(spkSeq(randperm(numel(spkSeq)))',circshift(sort(unique(spkSeq'),'ascend'),kcell),'type','Spearman');
                    if abs(rankCorr) > abs(rho_Shf(ishf,tt)) || isnan(rho_Shf(ishf,tt))
                        rho_Shf(ishf,tt) = rankCorr;
                    end
                end
            end
        end
    end
end
figure;plot(max(abs(rho),[],1))
hold on;plot(EXP.data.es.traj(1:10000)/100)
for icell = 1:size(CA1spikes,2)
    spkevt = find(CA1spikes(:,icell)>0);
    hold on;plot(spkevt,icell/size(CA1spikes,2)*ones(size(spkevt)),'o');
end

end