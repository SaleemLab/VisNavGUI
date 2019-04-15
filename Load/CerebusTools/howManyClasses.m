function k = howManyClasses(class_idx)
% HOWMANYCLASSES.M: Determine how many classes are in (a single, or all rows
% of) class_idx
    if size(class_idx,1) > 1
        for i = 1:size(class_idx,1)
            for j = 1:size(class_idx,2)
                a(i,j) = ~isempty(class_idx{i,j});
            end
        end
    else
        for i = 1:size(class_idx,2)
            a(i) = ~isempty(class_idx{i});
        end
    end
    k = sum(a,2);
end
% END howManyClasses