function nevclose

% Created DLR who knows when
% 2010-03-16 AZ Updated to include ns files

global pepNEV;

if ~isempty(pepNEV)
   if ~isfield(pepNEV,'ns')
      fclose(pepNEV.index.fid);
   end
   pepNEV = [];
end

return;