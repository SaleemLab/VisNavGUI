function saveSupraThreshold(animaltype,elec,iseries,iexp,gmm,UserName,T,...
   Xrealigned,iXraw,threshold)
% SAVESUPRATHRESHOLD.M: Save all events above a user-defined threshold as
% 998 units

gmm.icell     = 998;
gmm.dprime    = 0;
gmm.threshold = threshold(1);

nevsorted{1} = {[] T(iXraw)};
mwaves{1}    = mean(Xrealigned(:,iXraw),2);

nev2unit(animaltype,elec,iseries,iexp,gmm,nevsorted,mwaves,UserName);

end