% FIND_SIDX_NONEMPTY.M

sidx_nonempty = [];
for i = 1:5
   if ~isempty(sidx_tograph{i})
      sidx_nonempty = [sidx_nonempty i]; %#ok<AGROW>
   end
end