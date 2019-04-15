function showMeanTrace(hObject,eventdata)
% SHOWMEANTRACE.M
% Toggle mean trace visibility for a given class
%% TODO: un-globalize
global f1h

unit_class = get(hObject,'Tag');
old = strcmp(unit_class(7:end-3),'oldMeanWaves');
unit_class = str2double(unit_class(end-1));

if old
   if get(hObject,'Value')
      set(f1h.means.old.bord{unit_class},'Visible','on' );
      set(f1h.means.old.line{unit_class},'Visible','on' );
   else
      set(f1h.means.old.bord{unit_class},'Visible','off');
      set(f1h.means.old.line{unit_class},'Visible','off');
   end   
else
   if get(hObject,'Value')
      set(f1h.means.new.bord(unit_class),'Visible','on' );
      set(f1h.means.new.line(unit_class),'Visible','on' );
   else
      set(f1h.means.new.bord(unit_class),'Visible','off');
      set(f1h.means.new.line(unit_class),'Visible','off');
   end
end

% Clear existing highlights
if isfield(f1h.highlights,'unc') && ~isempty(f1h.highlights.unc(1)) && ishandle(f1h.highlights.unc(1))
    delete(f1h.highlights.unc(1));
    delete(f1h.highlights.unc(2));

    refresh;
end
if isfield(f1h.highlights,'pc') && ~isempty(f1h.highlights.pc{1}(2)) && ishandle(f1h.highlights.pc{1}(2))
    delete(f1h.highlights.pc{1}(2));
    delete(f1h.highlights.pc{1}(1));
    delete(f1h.highlights.pc{2}(2));
    delete(f1h.highlights.pc{2}(1));
    delete(f1h.highlights.pc{3}(2));
    delete(f1h.highlights.pc{3}(1));
    
    refresh;
end

% % Change Radio Button Selection
% set(f1h.radio.bringtotop(unit_class),'Value',1);

end %end function