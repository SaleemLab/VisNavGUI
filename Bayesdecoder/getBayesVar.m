function Xout = getBayesVar(X,latcorrection,sampleRate,Tsmth_win,smthtype,Fcircular,maxX)
X = circshift(X,[round(latcorrection/(1000/sampleRate)) 0]);
Xout = zeros(size(X));
if Fcircular
    for icell = 1:size(X,2)
        Xtemp = unwrap(X(:,icell)/maxX*2*pi)*maxX/(2*pi);
        Xtemp = smthInTime(Xtemp, sampleRate, Tsmth_win, 'same', [], smthtype);
        Xout(:,icell) = mod(Xtemp,maxX);
    end
else
    for icell = 1:size(X,2)
        Xout(:,icell) = smthInTime(X(:,icell), sampleRate, Tsmth_win, 'same', [], smthtype);
    end
end
end