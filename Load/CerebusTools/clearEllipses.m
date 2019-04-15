% CLEARELLIPSES.M
% AZ20090501: Clear existing highlights
if isfield(f1h,'plots') && isfield(f1h.plots(1),'ellipse') && ...
      ~isempty(f1h.plots(1).ellipse(1,1)) && ishandle(f1h.plots(1).ellipse(1,1))
   for i = 1:5
      if ~isempty(f1h.plots(i).ellipse)
         for j = 1:3
            if isfield(f1h.plots(i),'ellipse') && ~isempty(f1h.plots(i).ellipse) && ...
              ishandle(f1h.plots(i).ellipse(j,1))
               delete( f1h.plots(i).ellipse(j,1));
               delete( f1h.plots(i).ellipse(j,2));
            end
         end
      end
   end

   refresh;
end