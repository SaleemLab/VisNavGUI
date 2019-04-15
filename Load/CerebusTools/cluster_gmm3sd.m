function [f1h,gmm,d,sidx,sidx_tograph,C,S,U,V,Xrealigned] = cluster_gmm3sd( ...
          f1h,gmm,d,sidx,sidx_tograph,S,U,V,Xrealigned,                   ...
          iXraw,idx_tograph,xAxisInMsec,SamplingRateInKHZ,col,threshold,prev   )
% CLUSTER_GMM3SD.M FUNCTION

   %AZ 20081217
   set(f1h.fig.f1,'Pointer','watch');
   set(f1h.txt.status,'String','Auto-clustering...');
   drawnow;

   k  = get(f1h.pop.nclusters,'Value');

   f1h = sortnev2_controls_update(false,f1h);
    
   % AZ 2009-03-17: if loaded unit, xAxisInMsec will be 48 samples.  want 28
   % RESIZE
   switch size(xAxisInMsec,2)
      case 48
         % take only samples (7:34) (t = -0.1:1/30:0.8)
         Xrealigned = Xrealigned(7:34,:);
         xAxisInMsec = ((7:34)-10)/SamplingRateInKHZ; % center about thresholding point (10)
         [U,S,V] = svds(Xrealigned,3);
      case 28
         % do nothing
      otherwise
         warning('sortnev:cluster_gmm3sd:BadPCASize','Weird size of loaded PCA''s');
   end

	if isempty(gmm)
      % AZ20090303 %
      % K-means, then 3-D Gaussian Fit each
      C = cell(1,25); sidx = cell(1,5); gmm = []; d = zeros(1,25);
      for j = 1:25 % try 25x
         C{j}(iXraw,:) = kmeans(V(iXraw,:),k,'start','cluster');  %% ...  'replicates', 4

         for i = 1:4
            sidx{i} = find(C{j}==i);
         end

         if k == 1
            d = zeros(1,4);
         else
            d(j) = dprime(V,sidx);
         end

         set(f1h.txt.status,'String',['Auto-clustering... ',num2str(j*4),'%']);
         drawnow;
      end
      % Take result with best d'
      C = C{find(d==max(d),1)};
      for i = 1:4
         sidx{i}         = find(C==i);
         sidx_tograph{i} = sidx{i}(ismember(sidx{i},idx_tograph));
      end

      d = zeros(1,k);
      % if using kmeans, create GMM obj. also calculate d'
      for i = 1:k
   %       [gmm.obj.mu(i,:),gmm.obj.Sigma(i,:)] = normfit(V(find(C==i),:));
         % Create Gaussian Mixture Model Parameters
         gmm.mu(i,:)      = mean(V(logical(C==i),:),1);
         gmm.Sigma(:,:,i) =  cov(V(logical(C==i),:)  );
         d(i) = dprime(V,sidx{i},vertcat(sidx{setdiff(1:5,i)}));
      end
	end
      
   %AZ20090424
%   probThreshold = mvnpdf([0 0 0],[0 0 0],eye(3)); % 1x should = (2*pi)^(-3/2)
   probThreshold = zeros(1,5); P = zeros(size(V,1),k);
   for i = 1:k
      probThreshold(i) = mvnpdf(  gmm.mu(i,:),gmm.mu(i,:),gmm.Sigma(:,:,i));
      P(:,i)           = mvnpdf(V,gmm.mu(i,:),            gmm.Sigma(:,:,i));
   end
   % Threshold at 3 SD's (of largest distribution): ~99.7%
   probThreshold = (1-erf(3/sqrt(2)))*max(probThreshold);
   % C = index of maximum probability
   [a,C] = max(P,[],2); clear a;
   %indices where all probs for all classes are below threshold
   sidx{5} = find(sum(P<probThreshold,2)==k);
      
%       gmm.obj = gmdistribution(gmm.mu,gmm.Sigma,gmm.PComponents);
%       gmm = rmfield(gmm, {'mu'; 'Sigma'; 'PComponents'});

   %    C = cluster(gmm.obj,V);
   for i = 1:4
      sidx{i}         = find(C==i);
      sidx{i}         = setdiff(sidx{i},sidx{5});
      sidx_tograph{i} = sidx{i}(ismember(sidx{i},idx_tograph));
   end

   for i = 1:4
      d(i) = dprime(V,sidx{i},vertcat(sidx{setdiff(1:5,i)}));
   end

%    %   DEBUG   %
%    figure;hold on;
%    for i = [1:k 5]
%       scatter(V(sidx{i},1),V(sidx{i},2))
%    end
%    % END DEBUG %
   
   % AZ 20090316: RANK
   if k > 1
      [a,rank]               = sort(d,'descend');
      if sum(isnan(a))
         rank(isnan(a)) = [];
      end
   else
      % Ranks and sort by the proportion of elements in the class passed threshold
      rank = zeros(1,4);
      for i = 1:4
         rank(i)         = sum(ismember(sidx{i},find(iXraw)))/size(sidx{i},1);
      end
      rank(isnan(rank)) = [];
      [a,rank]               = sort(rank,'descend');
   end
   clear a;
   sidx(           sort(rank)  ) = sidx(               rank  );
   sidx_tograph(   sort(rank)  ) = sidx_tograph(       rank  );
   gmm.mu(         sort(rank),:) = gmm.mu(             rank,:);
   gmm.Sigma(:,:,  sort(rank)  ) = gmm.Sigma(    :,:,  rank  );
%    gmm.mu(         sort(rank),:) = gmm.obj.mu(         rank,:);
%    gmm.Sigma(:,:,  sort(rank)  ) = gmm.obj.Sigma(:,:,  rank  );
%    gmm.PComponents(sort(rank)  ) = gmm.obj.PComponents(rank  );
   d(              sort(rank)  ) = d(                  rank  );
   % END RANKING

%    gmm.obj = gmdistribution(gmm.mu,gmm.Sigma,gmm.PComponents);
%    gmm = rmfield(gmm, {'mu'; 'Sigma'; 'PComponents'});
%    gmm.obj
   
	[Xrealigned,idx_tograph,idx_tograph_only,sidx_tograph,sidx,f1h] = sortnev2_gui_update(...
    Xrealigned,V,iXraw,sidx_tograph,gmm,prev,xAxisInMsec,col,sidx,d,threshold,f1h,...
    'Clustering done!');
end
% END CLUSTER_GMM3SD FUNCTION