function createSynchFiles_Blackrock(nevDir)
global pepNEV

filepath_nev = ls([nevDir filesep '*.nev']);
%Will use the digital screen times from the photodiode as the aperiodic
%signal

filepath_ns5 = ls([nevDir filesep '*.ns5']);
%in ns5 recording:
% - channels #33 is PD 
% - channels #34 is licks 
% - channels #35 is rotary encoder


if size(filepath_nev, 1) > 1 || size(filepath_ns5, 1) > 1
    error('recordings should be split in different folders')
end
nsopen([nevDir filesep filepath_ns5(1,:)]);
nChannels = size(pepNEV.ns.Data.data,1);
PhotoDiode = pepNEV.ns.Data.data(nChannels-2,:)';
Licks = pepNEV.ns.Data.data(nChannels-1,:)';
RotEnc = pepNEV.ns.Data.data(nChannels,:)';
fulllength = length(pepNEV.ns.Data.data(1,:));

nevFileName = [nevDir filesep filepath_nev(1,:)];
openNEV(nevFileName,'nomat','overwrite');
nevopen(nevFileName);

screenTimes  = pepNEV.sync.timestamps;
screenTimes(find(diff(screenTimes)<100) + 1) = [];
Aperiodic = zeros(fulllength,1,'int16');
Aperiodic(screenTimes) = 1;

UpSampling = 1;

savedfilename = filepath_ns5(1:(strfind(filepath_ns5,'.')-1));
save([nevDir filesep 'Synch_BR_' savedfilename '.mat'], 'UpSampling', 'Aperiodic', 'PhotoDiode', 'Licks', 'RotEnc', '-v7.3');

nevclose;

clear global pepNEV

end