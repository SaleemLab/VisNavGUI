function S = updateFields(S,Supdate)
fields = fieldnames(Supdate);
for f = 1:numel(fields)
    S.(fields{f}) = Supdate.(fields{f});
end
end