function [Errmean,Latopt,latcorrectionlist] = PopEquivalentLatency
expt = getExperimentList;
nanimal = numel(expt);
nProbe = 2;

Tsmthwin = 15;
Xsmthwin = 4;
SpeedThreshold = 5;
Tsmthwin_dec = 250;
nDecbins = 100;
nspeedbins = 5;
neyebins = 1;
nphsbins = 1;
cellstr = 'goodonly';
filesuffix_EXP = ['Twin' num2str(Tsmthwin) '_' 'Xwin' num2str(Xsmthwin) '_' 'spdth' num2str(SpeedThreshold) '_' 'Decwin' num2str(Tsmthwin_dec) '_' 'nDecbins' num2str(nDecbins) '_' num2str(nspeedbins) 'speedbins' '_' num2str(neyebins) 'eyebins' '_' num2str(nphsbins) 'thetabins' '_' cellstr];
area_str = 'CA1V1';
datadir = 'D:\DATA\batch';

icount = zeros(1,nProbe);
Errmean = cell(nProbe,3);
Latopt = cell(nProbe,3);

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        for iprobe = 1:nProbe
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                    dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
                    S = load([dDIRname filesep 'latency_' filesuffix_EXP '.mat']);
                    latcorrectionlist = S.latcorrectionlist;
                    icount(iprobe) = icount(iprobe) +1;
                    for g = [2 1 3]
                        Errmean{iprobe,g}(icount(iprobe),:) = getCircularAverage(squeeze(nanmean(S.Err{iprobe,g},3))',0,1);
%                         Latopttemp = find(diff(sign(Errmean{iprobe,g}(icount(iprobe),:)-Errmean{iprobe,2}(icount(iprobe),:)))~=0,1,'first');
                        [~,Latopttemp] = min(abs(Errmean{iprobe,g}(icount(iprobe),:)-Errmean{iprobe,2}(icount(iprobe),:)));
                        if ~isempty(Latopttemp)
                            Latopt{iprobe,g}(icount(iprobe)) = latcorrectionlist(Latopttemp)+0.5*mean(diff(latcorrectionlist));%latcorrectionlist(Latopttemp);%
                        else
                            Latopt{iprobe,g}(icount(iprobe)) = latcorrectionlist(end);
                        end
                    end
                end
            end
        end
    end
end
end