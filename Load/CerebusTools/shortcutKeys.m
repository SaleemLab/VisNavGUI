function shortcutKeys(hObject,eventdata)
% SHORTCUTKEYS.M: for sortnev2
   %% TODO: un-globalize
   global f1h

   switch eventdata.Key
%       case 'return' % If no file is loaded, pressing enter when figure (only) has 
%                     %  focus will "click" load
%          if hObject == f1h.fig.f1 && isempty(nevopen_outcome)
%             % Pretend we clicked 'Load' in f1
%             loadnev(findobj('Tag','pshbtn.loadnev'),eventdata);
%          %              set(f1h.fig.f1,'KeyPressFcn',[]);
%          end
      case 'leftarrow'
         if strcmp(eventdata.Modifier,'control')
            set(f1h.f1g.f1,'CurrentObject',f1h.pshbtn.prevelec);
            sortnev2(elecrement_Callback);
         else
            uicontrol(f1h.pshbtn.prevelec);
         end
      case 'rightarrow'
         if strcmp(eventdata.Modifier,'control')
            set(f1h.f1g.f1,'CurrentObject',f1h.pshbtn.nextelec);
            sortnev2(elecrement_Callback);
         else
            uicontrol(f1h.pshbtn.nextelec);
         end
      case 'a'
         uicontrol(f1h.pshbtn.clusterAuto);
      case 'c'
         uicontrol(f1h.pshbtn.cluster);
      case 'e'
         uicontrol(f1h.pshbtn.nextelec);
      case 'f'
         uicontrol(f1h.chkbx.fixedsc);
      case 'k'
         uicontrol(f1h.pop.nclusters);
%       case 'l'
%          uicontrol(f1h.pshbtn.loadnev);
      case 'm'
         uicontrol(f1h.pshbtn.clusterManual);
      case 'n'
         uicontrol(f1h.txtbox.elecnum);
      case 'r'
         uicontrol(f1h.pshbtn.reload);
      case 's'
         if strcmp(eventdata.Modifier,'control')
            sortnev2(save_Callback);
         else
            uicontrol(f1h.pshbtn.save);
         end
      case 't'
         uicontrol(f1h.txtbox.thresh);
      case 'u'
         uicontrol(f1h.txtbox.username);
      otherwise
         %do nothing
   end
end