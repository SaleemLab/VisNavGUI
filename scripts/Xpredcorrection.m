function [predave_out,x] = Xpredcorrection(predave,predave_ref,Xrange,dx)
predave_ref(predave_ref - (1:Xrange)' > floor(Xrange/2)) = predave_ref(predave_ref - (1:Xrange)' > floor(Xrange/2)) - Xrange;
predave_ref(predave_ref - (1:Xrange)' < -floor(Xrange/2)) = predave_ref(predave_ref - (1:Xrange)' < -floor(Xrange/2)) + Xrange;
predave(predave - (1:Xrange)' > floor(Xrange/2)) = predave(predave - (1:Xrange)' > floor(Xrange/2)) - Xrange;
predave(predave - (1:Xrange)' < -floor(Xrange/2)) = predave(predave - (1:Xrange)' < -floor(Xrange/2)) + Xrange;

predave_interp = interp1(linspace(0,numel(predave),numel(predave)+1), [predave(1);predave], 0:dx:(Xrange-dx));predave_interp(isnan(predave_interp)) = 0;
predaveref_interp = interp1(linspace(0,numel(predave_ref),numel(predave_ref)+1), [predave_ref(1);predave_ref], 0:dx:(Xrange-dx));predaveref_interp(isnan(predaveref_interp)) = 0;
predave_out = zeros(size(predave_interp));
x = 0:dx:(Xrange-dx);
predaveref_interp = [predaveref_interp(1:end-1)-Xrange predaveref_interp predaveref_interp(2:end)+Xrange];
xrep = [x(1:end-1)-Xrange x x(2:end)+Xrange];
for i = 1:numel(x)
    idxmatch = min(find(abs(predaveref_interp-predave_interp(i)) <= min(abs(predaveref_interp-predave_interp(i)))));
    idxmatch = idxmatch(abs(idxmatch-(numel(x)-1) - i) == min(abs(idxmatch-(numel(x)-1) - i)));
    predave_out(i) = xrep(round(idxmatch)) - x(i);
end
end