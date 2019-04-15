function calculated = elecSubplotDestroyChildren(fignum,tagname,A,nClusters,calculated,elecs)
% ELECSUBPLOTDESTROYCHILDREN.M: Create Subplot for electrode function
    if nargin < 6
       elecs = A';
       elecs = elecs(~isnan(elecs))';
    end

    if isempty(nClusters) || isnan(nClusters)
        KsToTry = 2:4;
    else
        KsToTry = nClusters;
    end

    for elec = elecs
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
            if ~isnan(elec)
                h = allchild(findobj('Tag',['elec',num2str(elec),tagname{fignum}]));
                
                delete(h);
            end
%         end
    end
    
    calculated(fignum).D      = 0;
    calculated(fignum).Num    = 0;
    calculated(fignum).Traces = 0;
    calculated(fignum).Axes   = 0;  % make sure this is always 0
end
% END Create Subplot for electrode function