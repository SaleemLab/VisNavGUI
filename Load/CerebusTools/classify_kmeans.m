function [class_idx_temp,d_temp,best_i,Velec] = classify_kmeans(Xelec,iXrawelec,...
          class_idx_temp,k_temp)
% CLASSIFY_KMEANS.M: Classify (using k-means) function
% nClusters = # of classes

    % CLASSIFY: kmeans, then PCA->d'
    if k_temp > 1 && nnz(iXrawelec) > k_temp
        
%         C = kmeans(Xelec',k_temp,'start','cluster');  %% ...  'replicates', 4
        C = kmeans(Xelec(:,iXrawelec)',k_temp);  % 'start','cluster' uses 10% subset as first pass

        for a = 1:k_temp
            class_idx_temp{a} = find(C==a);
        end

        [U,S,Velec] = svds(Xelec(:,iXrawelec),3); % SVD
        % Apply svd from suprathreshold subset to create pca all traces
        Velec = (S\U'*Xelec)';
        
        % Calculate d' (for all permutations)
        [d_temp,best_i] = dprime(Velec(iXrawelec,:),class_idx_temp);

    else%if k_temp == 1 % i.e., if NOT clustering
        class_idx_temp{1} = (1:length(Xelec(:,iXrawelec)))';
        % make mean(d) = 0.5
        d_temp = 0.5;
        best_i = 1;
    end

end
% END Classify (using k-means) function