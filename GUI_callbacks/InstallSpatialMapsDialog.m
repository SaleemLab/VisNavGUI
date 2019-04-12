function InstallSpatialMapsDialog(GUI, EXP, P, SpatialMapsDialog)
if numel(EXP.Nav.contrastValues) > 1
    Fcontrast = true;
else
    Fcontrast = false;
end
if numel(EXP.Nav.gainValues) > 1
    Fgain = true;
else
    Fgain = false;
end
if numel(EXP.Nav.roomlengthValues) > 1
    FroomLength = true;
else
    FroomLength = false;
end
if numel(EXP.Nav.outcomeValues) > 1
    Foutcome = true;
else
    Foutcome = false;
end
% cbase = find(EXP.SubsetVal.contrast == mode(EXP.data.es.contrast));
% gbase = find(EXP.SubsetVal.gain == mode(EXP.data.es.gain));
% rbase = find(EXP.SubsetVal.roomlength == mode(EXP.data.es.roomLength));
% obase = find(EXP.SubsetVal.outcome == 2);
% P.PlotParamsMaps.ChosenContrast = cbase;
% dialog2.UpdateUIcontrol('Gain','String', strsplit(num2str(EXP.SubsetVal.gain)), 'max', numel(EXP.SubsetVal.gain), 'min', 0, 'Val', gbase); 
% P.PlotParamsMaps.ChosenGain = gbase;
% dialog2.UpdateUIcontrol('Roomlength','String', strsplit(num2str(EXP.SubsetVal.roomlength)), 'max', numel(EXP.SubsetVal.roomlength), 'min', 0, 'Val', rbase);
% P.PlotParamsMaps.ChosenRoomlength = rbase;
% dialog2.UpdateUIcontrol('Outcome','String', strsplit(num2str(EXP.SubsetVal.outcome)), 'max', numel(EXP.SubsetVal.outcome), 'min', 0, 'Val', obase);
% P.PlotParamsMaps.ChosenOutcome = obase;

SpatialMapsDialog.reset;
SpatialMapsDialog.getlistbox('cell#',EXP.Spk.CellListString,[0.5 6],'Callback',@(hObject, eventdata)SelectCell_callback(hObject, eventdata, GUI, EXP, P),'Min', 1, 'Max', length(EXP.Spk.CellListString));
SpatialMapsDialog.getcheckbox('Pool cells',P.PlotParamsMaps.FpoolCell,0.3,'Callback',@(hObject, eventdata)FPoolCell_callback(hObject, eventdata, GUI, EXP, P));
SpatialMapsDialog.getlistbox('X variable',P.PlotParamsMaps.PlotVarList,[0.3 2],'Callback',@(hObject, eventdata)SelectVarX_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',1);
SpatialMapsDialog.getlistbox('Y variable',P.PlotParamsMaps.PlotVarList,[0.3 2],'Callback',@(hObject, eventdata)SelectVarY_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',1);
SpatialMapsDialog.getlistbox('Plots',P.PlotParamsMaps.PlotObjList,[0.3 2],'Callback',@(hObject, eventdata)SelectPlotObj_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',numel(P.PlotParamsMaps.PlotObjList));
SpatialMapsDialog.getcheckbox('Disp. Mat.',P.PlotParamsMaps.Fdispmat,0.3,'Callback',@(hObject, eventdata)Fdispmat_callback(hObject, eventdata, GUI, EXP, P));
SpatialMapsDialog.getvariable('Speed Thresh.',P.PlotParamsMaps.speed_th,[0.3 0.3],'Callback',@(hObject, eventdata)Getspeedth_callback(hObject, eventdata, GUI, EXP, P));
SpatialMapsDialog.getvariable('Xbin size',P.PlotParamsMaps.Xbinsize,[0.3 0.3],'Callback',@(hObject, eventdata)GetXbinsize_callback(hObject, eventdata, GUI, EXP, P));
SpatialMapsDialog.getcheckbox('Combine conditions',P.PlotParamsMaps.FXcond,0.3,'Callback',@(hObject, eventdata)FCombinedCond_callback(hObject, eventdata, GUI, EXP, P));
if Fcontrast
    SpatialMapsDialog.getlistbox('Contrast',strsplit(num2str(EXP.Nav.contrastValues)),[0.3 2],'Callback',@(hObject, eventdata)SelectContrast_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',numel(EXP.Nav.contrastValues));
    SpatialMapsDialog.getcheckbox('Pool Contrasts',P.PlotParamsMaps.FpoolContrast,0.5,'Callback',@(hObject, eventdata)FPoolContrast_callback(hObject, eventdata, GUI, EXP, P));
end
if Fgain
    SpatialMapsDialog.getlistbox('Gain',strsplit(num2str(EXP.Nav.gainValues)),[0.3 2],'Callback',@(hObject, eventdata)SelectGain_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',numel(EXP.Nav.gainValues));
    SpatialMapsDialog.getcheckbox('Pool Gains',P.PlotParamsMaps.FpoolGain,0.5,'Callback',@(hObject, eventdata)FPoolGain_callback(hObject, eventdata, GUI, EXP, P));
end
if FroomLength
    SpatialMapsDialog.getlistbox('Roomlength',strsplit(num2str(EXP.Nav.roomlengthValues)),[0.3 2],'Callback',@(hObject, eventdata)SelectRoomlength_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',numel(EXP.Nav.roomlengthValues));
    SpatialMapsDialog.getcheckbox('Pool Roomlengths',P.PlotParamsMaps.FpoolRoomlength,0.5,'Callback',@(hObject, eventdata)FPoolRoomlength_callback(hObject, eventdata, GUI, EXP, P));
end
if Foutcome
    SpatialMapsDialog.getlistbox('Outcome',strsplit(num2str(EXP.Nav.outcomeValues)),[0.3 2],'Callback',@(hObject, eventdata)SelectOutcome_callback(hObject, eventdata, GUI, EXP, P),'Min',1,'Max',numel(EXP.Nav.outcomeValues));
    SpatialMapsDialog.getcheckbox('Pool Outcomes',P.PlotParamsMaps.FpoolOutcome,0.5,'Callback',@(hObject, eventdata)FPoolOutcome_callback(hObject, eventdata, GUI, EXP, P));
end


%Callback functions
    function SelectCell_callback(hObject, eventdata, GUI, EXP, P)
        sel = get(hObject,'Value');
        P.PlotParamsMaps.ChosenCell = sel;
        if P.PlotParamsMaps.FpoolCell
            P.PlotParamsMaps.ChosenCell = P.PlotParamsMaps.ChosenCell(:);
        else
            P.PlotParamsMaps.ChosenCell = P.PlotParamsMaps.ChosenCell(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
        uicontrol(hObject);
    end

    function FPoolCell_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.FpoolCell = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolCell
            P.PlotParamsMaps.ChosenCell = P.PlotParamsMaps.ChosenCell(:);
        else
            P.PlotParamsMaps.ChosenCell = P.PlotParamsMaps.ChosenCell(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectVarX_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsMaps.ChosenVarX = str(val);
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectVarY_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsMaps.ChosenVarY = str(val);
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectPlotObj_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsMaps.ChosenObj = str(val);
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function Fdispmat_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.Fdispmat = get(hObject,'Value');
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function Getspeedth_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsMaps.speed_th = val;
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function GetXbinsize_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsMaps.Xbinsize = val;
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function FCombinedCond_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.FXcond = get(hObject,'Value');
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectContrast_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.ChosenContrast = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolContrast
            P.PlotParamsMaps.ChosenContrast = P.PlotParamsMaps.ChosenContrast(:);
        else
            P.PlotParamsMaps.ChosenContrast = P.PlotParamsMaps.ChosenContrast(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function FPoolContrast_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.FpoolContrast = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolContrast
            P.PlotParamsMaps.ChosenContrast = P.PlotParamsMaps.ChosenContrast(:);
        else
            P.PlotParamsMaps.ChosenContrast = P.PlotParamsMaps.ChosenContrast(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectGain_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.ChosenGain = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolGain
            P.PlotParamsMaps.ChosenGain = P.PlotParamsMaps.ChosenGain(:);
        else
            P.PlotParamsMaps.ChosenGain = P.PlotParamsMaps.ChosenGain(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function FPoolGain_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.FpoolGain = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolGain
            P.PlotParamsMaps.ChosenGain = P.PlotParamsMaps.ChosenGain(:);
        else
            P.PlotParamsMaps.ChosenGain = P.PlotParamsMaps.ChosenGain(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectRoomlength_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.ChosenRoomlength = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolRoomlength
            P.PlotParamsMaps.ChosenRoomlength = P.PlotParamsMaps.ChosenRoomlength(:);
        else
            P.PlotParamsMaps.ChosenRoomlength = P.PlotParamsMaps.ChosenRoomlength(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function FPoolRoomlength_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.FpoolRoomlength = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolRoomlength
            P.PlotParamsMaps.ChosenRoomlength = P.PlotParamsMaps.ChosenRoomlength(:);
        else
            P.PlotParamsMaps.ChosenRoomlength = P.PlotParamsMaps.ChosenRoomlength(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function SelectOutcome_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.ChosenOutcome = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolOutcome
            P.PlotParamsMaps.ChosenOutcome = P.PlotParamsMaps.ChosenOutcome(:);
        else
            P.PlotParamsMaps.ChosenOutcome = P.PlotParamsMaps.ChosenOutcome(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end

    function FPoolOutcome_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsMaps.FpoolOutcome = get(hObject,'Value');
        if P.PlotParamsMaps.FpoolOutcome
            P.PlotParamsMaps.ChosenOutcome = P.PlotParamsMaps.ChosenOutcome(:);
        else
            P.PlotParamsMaps.ChosenOutcome = P.PlotParamsMaps.ChosenOutcome(:)';
        end
        UpdateplotBehavSpks(P.PlotParamsMaps,EXP,GUI);
    end
end