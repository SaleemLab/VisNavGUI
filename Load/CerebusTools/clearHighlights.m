% CLEARHIGHLIGHTS.M

% Clear existing highlights
% if isfield(f1h,'highlights')
   if isfield(f1h.highlights,'unc') && ~isempty(f1h.highlights.unc(1)) && ...
                                       ishandle(f1h.highlights.unc(1))
       delete(f1h.highlights.unc(1));
       delete(f1h.highlights.unc(2));

       refresh;
   end
   if isfield(f1h.highlights,'pc') && ~isempty(f1h.highlights.pc{1}(2)) && ...
                                      ishandle(f1h.highlights.pc{1}(2))
       delete(f1h.highlights.pc{1}(2));
       delete(f1h.highlights.pc{1}(1));
       delete(f1h.highlights.pc{2}(2));
       delete(f1h.highlights.pc{2}(1));
       delete(f1h.highlights.pc{3}(2));
       delete(f1h.highlights.pc{3}(1));

       refresh;
   end
% end