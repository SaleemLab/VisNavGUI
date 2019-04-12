function [relErr,Err,Err_Shf,nspdbinslist] = PopSpeednbins
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
relErr = cell(nProbe,3);
Err = cell(nProbe,3);
Err_Shf = cell(nProbe,3);

for ianimal = 1:nanimal
    animalname = expt(ianimal).animal;
    for iseries = 1:numel(expt(ianimal).series)
        disp([animalname num2str(expt(ianimal).series{iseries})]);
        for iprobe = 1:nProbe
            if ((iprobe == 1 && expt(ianimal).goodCA1{iseries} == 1) || (iprobe == 2 && expt(ianimal).goodV1{iseries} == 1))
                if ~isempty(strfind(area_str,expt(ianimal).area{iseries}))
                    dDIRname = [datadir  filesep animalname filesep num2str(expt(ianimal).series{iseries}) filesep 'processed'];
                    S = load([dDIRname filesep 'Speedbins_' filesuffix_EXP '.mat']);
                    nspdbinslist = S.nspdbinslist;
                    Prange = 100;
                    Xrange = 100;
                    icount(iprobe) = icount(iprobe) +1;
                    for g = [2 1 3]
                        for ibins = 1:numel(nspdbinslist)
                            if ~isempty(S.Xsmth0{iprobe,g,ibins})
                                meantraj = NaN(Xrange,1);
                                meantraj_Shf = NaN(Xrange,1);
                                for xx = 1:Xrange
                                    if sum(S.Xsmth0{iprobe,g,ibins}==xx)>0
                                        meantraj(xx) = mod(Prange + Prange/(2*pi)*circ_mean(2*pi/Prange*S.MeanDecodedPosition{iprobe,g,ibins}(S.Xsmth0{iprobe,g,ibins}==xx & ~isnan(S.MeanDecodedPosition{iprobe,g,ibins}))),Prange);
                                        meantraj_Shf(xx) = mod(Prange + Prange/(2*pi)*circ_mean(2*pi/Prange*S.MeanDecodedPosition_Shf{iprobe,g,ibins}(S.Xsmth0_Shf{iprobe,g,ibins}==xx & ~isnan(S.MeanDecodedPosition_Shf{iprobe,g,ibins}))),Prange);
                                    end
                                end
                                meandecErr2 = circ_dist(2*pi/Prange*S.MeanDecodedPosition{iprobe,g,ibins},2*pi/Prange*meantraj(S.Xsmth0{iprobe,g,ibins}));
                                Err{iprobe,g}(icount(iprobe),ibins) = Prange/(2*pi)*circ_std(meandecErr2);
                                meandecErr2_Shf = circ_dist(2*pi/Prange*S.MeanDecodedPosition_Shf{iprobe,g,ibins},2*pi/Prange*meantraj_Shf(S.Xsmth0_Shf{iprobe,g,ibins}));
                                Err_Shf{iprobe,g}(icount(iprobe),ibins) = Prange/(2*pi)*circ_std(meandecErr2_Shf);
                                relErr{iprobe,g}(icount(iprobe),ibins) = Err_Shf{iprobe,g}(icount(iprobe),ibins) - Err{iprobe,g}(icount(iprobe),ibins);
                            end
                        end
%                         relErr{iprobe,g}(icount(iprobe),:) = S.maxdecErrAll_Shf{iprobe,g} - S.maxdecErrAll{iprobe,g};
%                         Err{iprobe,g}(icount(iprobe),:) = S.maxdecErrAll{iprobe,g};
%                         Err_Shf{iprobe,g}(icount(iprobe),:) = S.maxdecErrAll_Shf{iprobe,g};
                    end
                end
            end
        end
    end
end
end