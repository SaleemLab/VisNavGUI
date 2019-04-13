function es = combineTwoexpts(es, es2)
if isempty(es)
    es = es2;
else
    fields = fieldnames(es);
    exceptionfields = {'sampleTimes', 'spikeTimes', 'trialID'};
    for i = 1:numel(fields)
        if sum(strcmp(fields{i}, exceptionfields)) == 0
            [~, catdim] = max(size(es.(fields{i})));
            es.(fields{i}) = cat(catdim, es.(fields{i}), es2.(fields{i}));
        end
    end
    %concatanating exceptions corresponding to cumulative fields
    [~, catdim] = max(size(es.sampleTimes));
    es.sampleTimes = cat(catdim, es.sampleTimes, es2.sampleTimes + es.sampleTimes(end));
    
    if isfield(es, 'spikeTimes')
        [~, catdim] = max(size(es.spikeTimes));
        es.spikeTimes = cat(catdim, es.spikeTimes, es2.spikeTimes + es.sampleTimes(end));
    end
    
    if isfield(es, 'trialID')
        [~, catdim] = max(size(es.trialID));
        es.trialID = cat(catdim, es.spikeTimes, es2.trialID + es.trialID(end));
    end
end