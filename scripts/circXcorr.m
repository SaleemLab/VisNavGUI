function Xcorr_circ = circXcorr(map,mapbase)

mapcorr = map - nanmean(map);
mapcorr = mapcorr./sqrt(sum(mapcorr.^2));
mapbasecorr = mapbase - nanmean(mapbase);
mapbasecorr = mapbasecorr./sqrt(sum(mapbasecorr.^2));
n = numel(mapbase);
Xcorr_circ = zeros(1,numel(mapbase));

mapXcorr_circ = mapbasecorr(:)*(mapcorr(:)')/sqrt(sum(mapcorr.^2)*sum(mapbasecorr.^2));
Xcorr_circ(1) = sum(diag(mapXcorr_circ));
for xshift = 1:n-1
    Xcorr_circ(xshift+1) = sum(diag(mapXcorr_circ,xshift)) + sum(diag(mapXcorr_circ,xshift-n));
end
Xcorr_circ = circshift(Xcorr_circ,floor(n/2)-1);

% xshiftlim = floor(numel(mapbase)/2);
% ishift = 0;
% for xshift = -xshiftlim:xshiftlim-1
%     ishift = ishift + 1;
%     Xcorr_circ(ishift) = sum(mapcorr.*circshift(mapbasecorr,xshift))/sqrt(sum(mapcorr.^2)*sum(mapbasecorr.^2));
% end