function f1h = sortnev2_controls_update(up_or_down,f1h,d,col,sidx_tograph)

if up_or_down
   sidx_nonempty = [];
   for i = 1:5
      if ~isempty(sidx_tograph{i})
         sidx_nonempty = [sidx_nonempty i]; %#ok<AGROW>
      end
   end
   
   for i = sidx_nonempty
      set(f1h.chkbx.unitsave(    i),'Enable','On');
      set(f1h.chkbx.showtraces(  i),'Enable','On');
      set(f1h.chkbx.showtraces(  i),'Value' ,   1);
      set(f1h.chkbx.meanwave(    i),'Enable','On');
      if isfield(f1h.means,'old') && i <= max(size(f1h.means.old.line,2))
         set(f1h.chkbx.oldMeanWaves(i),'Enable','On');
      end
      set(f1h.radio.bringtotop(  i),'Enable','On');
      set(f1h.txtbox.unitnum(    i),'Enable','On');
      if i ~= 5
      set(f1h.txtbox.unitnum(    i),'String',num2str(i));
      set(f1h.txt.unitdprimeVal( i),'String',num2str(d(i),'%0.3g'));
      else
      set(f1h.txtbox.unitnum(    i),'String','unc');
      end
      set(f1h.txtbox.unitnum(    i),'BackgroundColor',col{i});
   end
   set(f1h.radio.bringtotop(min(sidx_nonempty)),'Value',1);
   if     exist('sidx_tograph','var') && isempty(sidx_tograph{4})
      set(f1h.pshbtn.addClearEllipses,'String','Add Ellipse')
   elseif exist('sidx_tograph','var') && ~isempty(sidx_tograph{4})
      set(f1h.pshbtn.addClearEllipses,'String','Clear Ellipses')
   end
else   %AZ20090501: disable 'save' etc buttons
   for i = 1:5
      set(f1h.chkbx.unitsave(    i),'Enable','Off');
      set(f1h.chkbx.showtraces(  i),'Enable','Off','Value' ,    0);
      set(f1h.chkbx.meanwave(    i),'Enable','Off','Value' ,    0);
      set(f1h.chkbx.oldMeanWaves(i),'Enable','Off','Value' ,    0);
      set(f1h.radio.bringtotop(  i),'Enable','Off');
      set(f1h.txtbox.unitnum(    i),'BackgroundColor',[1 1 1]);
      set(f1h.txtbox.unitnum(    i),'Enable','Off','String',   '');
      set(f1h.txt.unitdprimeVal( i),'String',   '');
   end
end