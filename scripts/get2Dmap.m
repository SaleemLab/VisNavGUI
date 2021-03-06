function [model, spred] = get2Dmap(zua, Xvariable, Yvariable, scaling, nGridX, binsX, nGridY, binsY, CVO, iter, sampFreq, winX, winY, FoptiSmooth, FcircularX, FcircularY,FuniformY,Xvariable_occ, Yvariable_occ)
% function to fit and test a 1D model, find the optimal smoothing
% function too
if nargin < 14
    FoptiSmooth = true;
end
if nargin < 15
    FcircularX = true;
    FcircularY = true;
end
if nargin < 17
    FuniformY = false;
end
if nargin < 18
    Xvariable_occ = [];
    Yvariable_occ = [];
end
X_train = Xvariable(CVO.train{iter}>0);
X_test  = Xvariable(CVO.test{iter}>0);
% X_cv    = Xvariable(CVO.cv{iter}>0);

Y_train = Yvariable(CVO.train{iter}>0);
Y_test  = Yvariable(CVO.test{iter}>0);
% Y_cv    = Yvariable(CVO.cv{iter}>0);

if ~isempty(Xvariable_occ)
    X_train_occ = Xvariable_occ(CVO.train{iter}>0);
    X_test_occ  = Xvariable_occ(CVO.test{iter}>0);
    Y_train_occ = Yvariable_occ(CVO.train{iter}>0);
    Y_test_occ  = Yvariable_occ(CVO.test{iter}>0);
else
    X_train_occ = X_train;
    X_test_occ  = X_test;
    Y_train_occ = Y_train;
    Y_test_occ  = Y_test;    
end

% for icell = 1:numCells
if nargin<11
    sampFreq = 200;
end
if nargin < 12
    winX = 10;
    smthX = 0;
    winY = 10;
    smthY = 0;
else 
    smthX = 1;
    smthY = 1;
end

% if smth
%     zua_train = smthInTime(zua(CVO.train{iter}>0), sampFreq, win);
%     zua_test  = smthInTime(zua(CVO.test{iter}>0), sampFreq, win);
%     zua_cv    = smthInTime(zua(CVO.cv{iter}>0), sampFreq, win);
% else
    zua_train = zua(CVO.train{iter}>0).*scaling(CVO.train{iter}>0);
    zua_test  = zua(CVO.test{iter}>0).*scaling(CVO.test{iter}>0);
%     zua_cv    = zua(CVO.cv{iter}>0).*scaling(CVO.cv{iter}>0);
% end
test = zua(CVO.test{iter});

% get spike count map
scMap = full(sparse(Y_train, X_train, (zua_train)', nGridY, nGridX));
% get occupancy map
occMap = full(sparse(Y_train_occ, X_train_occ, 1, nGridY, nGridX));
if FuniformY
    n = size(occMap,1);
    occMap = repmat(mean(occMap,1),[n 1]);
end

if FcircularX
    scMap = repmat(scMap,[1 3]);
    occMap = repmat(occMap,[1 3]);
else
    scMap = [zeros(size(scMap)) scMap zeros(size(scMap))];
    occMap = [zeros(size(occMap)) occMap zeros(size(occMap))];
end
if FcircularY
    scMap = repmat(scMap,[3 1]);
    occMap = repmat(occMap,[3 1]);
else
    scMap = [zeros(size(scMap));scMap;zeros(size(scMap))];
    occMap = [zeros(size(occMap));occMap;zeros(size(occMap))];
end

n1X = nGridX;
gridsX = 1:n1X;
n1Y = nGridY;
gridsY = 1:n1Y;

idxX = floor(numel(gridsX)/winX); 
idxY = floor(numel(gridsY)/winY);

if winX > 0 && ~isnan(winX)
    model.swinX = gridsX(min(idxX,end));
else
    model.swinX = NaN;
end
model.binsX = binsX;
if winY > 0 && ~isnan(winY)
    model.swinY = gridsY(min(idxY,end));
else
    model.swinY = NaN;
end
model.binsY = binsY;

model.tuning = special_smooth_2d(scMap, [1./model.swinY 1./model.swinX], 0, 0, [n1Y n1X])...
             ./special_smooth_2d(occMap, [1./model.swinY 1./model.swinX], 0, 0, [n1Y n1X]);
         
model.tuning = model.tuning(:,nGridX+1:2*nGridX);
model.tuning = model.tuning(nGridY+1:2*nGridY,:);
 
tuningvec = model.tuning(:);
pred  = tuningvec((X_test-1).*n1Y + Y_test);

spred = pred;

[model.EV, model.corr, model.L, model.Q, model.train_mean] = calCrossValExpVar(zua_train, zua_test, spred, test, pred);


P_i = 1./numel(model.tuning); %occMap./sum(occMap(:));% 
R_mean = nanmean(model.tuning(:)); % nanmean(X_train);
R_i = model.tuning; % scMap;

model.skaggs = sum(sum( P_i .* (R_i/R_mean) .* log2(R_i/R_mean),1),2);

end