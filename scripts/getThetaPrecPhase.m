function [phsoffset, amp] = getThetaPrecPhase(vec,phs0)
vec = vec(:)';
nphsbins = numel(vec);
Phi = (0.5:(nphsbins-0.5))/nphsbins*2*pi;
Phi_num = sum((vec-mean(vec)).*sin(Phi));
Phi_den = sum((vec-mean(vec)).*cos(Phi));
phsoffset = 180 - 360*atan2(Phi_den,Phi_num)/(2*pi);
if nargin < 2
    phs0 = phsoffset;
end
amp = sum((vec-mean(vec)).*sin(Phi + (180-phs0)/360*2*pi))/sum(sin(Phi).^2);
end