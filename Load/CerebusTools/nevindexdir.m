function nevindexdir(nevdir)

list = dir([nevdir '\*.nev']);

for(i=1:length(list))
    fn = [nevdir '\' list(i).name];
    if(list(i).name(1)=='u')
       % try
            nevindex(fn)
            disp(sprintf('%s has been indexed',fn))
        %catch
            disp(sprintf('Warning: Could not index %s',fn))
        %end
    end
end
