function [map,x,mapstd,mapsem] = fast1Dmap(X,Y,dx,scaling,nbXbinsmth,Fcircular,nbbins)
if ~isempty(X)
    if nargin < 3
        dx = 2.5;%(max(X)/nbbins);
        scaling = 1;
    end
    if nargin < 4
        scaling = 1;
    end
    if nargin < 6
        Fcircular = true;
    end
    if nargin < 7
        nbbins = [];
    end
    if isnan(dx)
        if isempty(nbbins)
            nbbins = max(X);
        end
        dx = 1;
        x = 1:nbbins;
    else
        if isempty(nbbins)
            nbbins = round(max(floor(X/dx)+1));%round(100/dx);%
        end
        X = floor(X/dx)+1;
        x = 0:dx:dx*(nbbins-1);
    end
    
    if nargin < 5
        nbXbinsmth = [];
    end
    if isempty(nbXbinsmth)
        nbXbinsmth = nbbins;
    end
      
    Y = Y.*scaling;
    scMap = full(sparse(X, 1, Y', nbbins, 1));
    occMap = full(sparse(X, 1, 1, nbbins, 1));
    if nargout > 2
        scMapsqr = full(sparse(X, 1, (Y.^2)', nbbins, 1));
    end
    if Fcircular
        scMap = [scMap;scMap;scMap];
        occMap = [occMap;occMap;occMap];
        map = special_smooth_1d(scMap,1/(3*nbXbinsmth),0,(3*nbbins))./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins));
        if nargout > 2
            scMapsqr = [scMapsqr;scMapsqr;scMapsqr];
            mapsqr =  special_smooth_1d(scMapsqr,1/(3*nbXbinsmth),0,(3*nbbins))./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins));
        end
    else
        map = special_smooth_1d(scMap,1/(nbXbinsmth),0,(nbbins))./special_smooth_1d(occMap,1/(nbXbinsmth),0,nbbins);
        if nargout > 2
            mapsqr =  special_smooth_1d(scMapsqr,1/(nbXbinsmth),0,(nbbins))./special_smooth_1d(occMap,1/(nbXbinsmth),0,(nbbins));
        end
    end
    if nargout > 2
        mapstd = (mapsqr - map.^2).^0.5;%std
        mapsem = ((mapsqr - map.^2).^0.5)./special_smooth_1d(occMap,1/(3*nbXbinsmth),0,(3*nbbins)).^0.5;%sem
    end
   
    if Fcircular
        map = map(nbbins+1:2*nbbins);
        if nargout > 2
            mapstd = mapstd(nbbins+1:2*nbbins);
            mapsem = mapsem(nbbins+1:2*nbbins);
        end
    end
else
    map = [];
    x = [];
    mapstd = [];
    mapsem = [];
end

end
