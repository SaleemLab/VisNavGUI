function [class_idx,d,fh,calculated,k,V] = cluster_electrodes(d,X,nClusters,...
      class_idx,A,iXraw,fh,calculated,V,elecs)
% CLUSTER_ELECTRODES.M: Cluster all electrodes
%     nClusters = get(fh.expnev1.cluster.popup,'Value') + 1;
    % Initialize
    if nargin < 10
       elecs = A';
       elecs = elecs(~isnan(elecs))';
    end

    if isempty(nClusters) || isnan(nClusters)
        KsToTry = 2:4;
    else
        KsToTry = nClusters;
    end
    
    classWasUpdated = false;

    for elec = elecs
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j);
            if ~isnan(elec) && ~isempty(X{elec}) && nnz(iXraw{elec}) > max(KsToTry)
                
                class_idx_temp    = cell(1,4);
                d_temp    = []; %k_temp    = [];

                for k_temp = KsToTry
                    [class_idx_temp,d_temp,best_i,V{elec}] = classify_kmeans(...
                               X{elec},iXraw{elec},class_idx_temp,k_temp);
                    if d_temp >= d(elec)
                        d(elec)         =         d_temp;
                        % BRING BEST_I CLUSTER TO END (so make swap in class_idx...
                        % + any others?)
                        if best_i ~= k_temp
                            for a = 1:k_temp
                                if     a == best_i
                                    class_idx(elec,a) = class_idx_temp(k_temp);
                                elseif a == k_temp
                                    class_idx(elec,a) = class_idx_temp(best_i);
                                else
                                    class_idx(elec,a) = class_idx_temp(a);
                                end
                            end
                        else
                            class_idx(elec,:) = class_idx_temp;
                        end
                        classWasUpdated = true;
                    end
                end
               
               if ~classWasUpdated
                  class_idx(elec,:) = class_idx_temp;
               end
                
               % AZ 20090812: adjust class_idx for iXraw indexing
               b = find(iXraw{elec});
               for a = 1:size(class_idx,2)
                  if ~isempty(class_idx{elec,a})
                     if all(ismember(class_idx{elec,a},b))
%                      try
                        class_idx{elec,a} = b(class_idx{elec,a});
%                      catch
                     else
                        [junk,class_idx{elec,a}] = intersect(b,class_idx{elec,a});
%                         %% was a warning, does it still need to be one?
%                         %% why does this happen?
%                         fprintf('Indexing for elec #%g is messed up.\n',elec);
%                      end
                     end
                  end
               end
            else
               d(elec) = 0;
            end
            
%             set(get(findobj(gcf,'Tag','status')),'String',['Clustering Electrode data... ',num2str((i-1)*10+j),'%']);
            set(get(findobj(gcf,'Tag','status')),'String',['Clustering Electrode data... ',...
               num2str(round(100*find(elec==elecs)/size(elecs,2))),'%']);
            drawnow;
%         end
    end
    
%     %Debug
%     for i = 1:96
%         for j = 1:4
%             a(i,j) = size(class_idx{i,j},1);
%         end
%     end
%     if ~isempty(find(sum(a,2)~=nwaves))
%         disp(['electrodes with less than ',num2str(nwaves),' points classified:']);
%         disp(find(sum(a,2)~=nwaves));
%     end
%     % END Debug
    
    k = howManyClasses(class_idx);
    
    if nClusters == 1 % i.e., if NOT clustering
        % make dMax-dMin nonzero (mean is 0.5)
        d(1)   = 0;
        d(end) = 1;
    end
    
    % Keep track of what was calculated during a given clustering session...
    calculated(1).D      = 0;
    calculated(1).Num    = 0;
    calculated(1).Traces = 0;
    calculated(1).Axes   = 0;  % make sure this is always 0
    calculated(2).D      = 0;
    calculated(2).Num    = 0;
    calculated(2).PCA    = 0;
    calculated(2).Axes   = 0;  % make sure this is always 0
end
% END Cluster electrodes function