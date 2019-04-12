function fitthetadec
meandecErr = Prange/(2*pi)*circ_dist(2*pi/Prange*EXP.Bayes.MeanDecodedPosition0{selprobe,1}(tidx & spdidx & phsidx & goodidx),2*pi/Prange*X(tidx & spdidx & phsidx & goodidx));
phsval = round(mod(EXP.Bayes.LFPphase{1},360));
phsvaltemp = phsval(tidx & spdidx & phsidx & goodidx);
phsvaltemp = round(phsvaltemp/30);
phsvaltemp = phsvaltemp(randperm(numel(phsvaltemp)));
CVO = crossValPartition(ones(1,numel(phsvaltemp)),20);
meanpredErr = NaN(size(meandecErr));
improv = [];
for iter = 1:20
map = fast1Dmap(phsvaltemp(CVO.train{iter}),meandecErr(CVO.train{iter}),1,1,NaN,1);
meanpredErr(CVO.test{iter}) = map(phsvaltemp(CVO.test{iter})+1);
improv(iter) = Prange/(2*pi)*(circ_std(2*pi/Prange*meandecErr(CVO.test{iter})) - circ_std(circ_dist(2*pi/Prange*meanpredErr(CVO.test{iter}),2*pi/Prange*meandecErr(CVO.test{iter}))));
end