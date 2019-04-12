function InstallBehaviorDialog(GUI, EXP, P, BehaviorDialog)

if numel(EXP.SubsetVal.contrast) > 1
    Fcontrast = true;
else
    Fcontrast = false;
end
if numel(EXP.SubsetVal.gain) > 1
    Fgain = true;
else
    Fgain = false;
end
if numel(EXP.SubsetVal.roomlength) > 1
    FroomLength = true;
else
    FroomLength = false;
end
if numel(EXP.SubsetVal.outcome) > 1
    Foutcome = true;
else
    Foutcome = false;
end

BehaviorDialog.reset;
BehaviorDialog.getlistbox('Plots',P.PlotParamsBehavior.PlotObjList, [0.3 5],'Min',1,'Max',numel(P.PlotParamsBehavior.PlotObjList), 'Callback', @(hObject, eventdata)SelectBehPlotObj_callback(hObject, eventdata, GUI, EXP, P));
BehaviorDialog.getvariable('Speed Thresh.', P.PlotParamsBehavior.speed_th, [0.3 0.3], 'Callback', @(hObject, eventdata)GetBehspeedth_callback(hObject, eventdata, GUI, EXP, P));
if Fcontrast
    BehaviorDialog.getlistbox('Contrast',strsplit(num2str(EXP.SubsetVal.contrast)), [0.3 3],'Min',1,'Max', numel(EXP.SubsetVal.contrast), 'Callback', @(hObject, eventdata)SelectBehContrast_callback(hObject, eventdata, GUI, EXP, P));
    BehaviorDialog.getcheckbox('Pool Contrasts', P.PlotParamsBehavior.FpoolContrast, 0.5, 'Callback', @(hObject, eventdata)FPoolBehContrast_callback(hObject, eventdata, GUI, EXP, P));
end
if Fgain
    BehaviorDialog.getlistbox('Gain',strsplit(num2str(EXP.SubsetVal.gain)), [0.3 2],'Min',1,'Max',numel(EXP.SubsetVal.gain), 'Callback', @(hObject, eventdata)SelectBehGain_callback(hObject, eventdata, GUI, EXP, P));
    BehaviorDialog.getcheckbox('Pool Gain', P.PlotParamsBehavior.FpoolGain, 0.5, 'Callback', @(hObject, eventdata)FPoolBehGain_callback(hObject, eventdata, GUI, EXP, P));
end
if FroomLength
    BehaviorDialog.getlistbox('Roomlength',strsplit(num2str(EXP.SubsetVal.roomlength)), [0.3 2],'Min', 1, 'Max', numel(EXP.SubsetVal.roomlength), 'Callback', @(hObject, eventdata)SelectBehRoomlength_callback(hObject, eventdata, GUI, EXP, P));
    BehaviorDialog.getcheckbox('Pool Roomlength', P.PlotParamsBehavior.FpoolRoomlength, 0.5, 'Callback', @(hObject, eventdata)FPoolBehRoomlength_callback(hObject, eventdata, GUI, EXP, P));
end
if Foutcome
    BehaviorDialog.getlistbox('Outcome',strsplit(num2str(EXP.SubsetVal.outcome)), [0.3 2],'Min', 1, 'Max', numel(EXP.SubsetVal.outcome), 'Callback', @(hObject, eventdata)SelectBehOutcome_callback(hObject, eventdata, GUI, EXP, P));
    BehaviorDialog.getcheckbox('Pool Outcome', P.PlotParamsBehavior.FpoolOutcome, 0.5, 'Callback', @(hObject, eventdata)FPoolBehOutcome_callback(hObject, eventdata, GUI, EXP, P));
end


%Callback functions
    function SelectBehPlotObj_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsBehavior.ChosenObj = str(val);
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function GetBehspeedth_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsBehavior.speed_th = val;
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function SelectBehContrast_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.ChosenContrast = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolContrast
            P.PlotParamsBehavior.ChosenContrast = P.PlotParamsBehavior.ChosenContrast(:);
        else
            P.PlotParamsBehavior.ChosenContrast = P.PlotParamsBehavior.ChosenContrast(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function SelectBehGain_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.ChosenGain = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolGain
            P.PlotParamsBehavior.ChosenGain = P.PlotParamsBehavior.ChosenGain(:);
        else
            P.PlotParamsBehavior.ChosenGain = P.PlotParamsBehavior.ChosenGain(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function SelectBehRoomlength_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.ChosenRoomlength = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolRoomlength
            P.PlotParamsBehavior.ChosenRoomlength = P.PlotParamsBehavior.ChosenRoomlength(:);
        else
            P.PlotParamsBehavior.ChosenRoomlength = P.PlotParamsBehavior.ChosenRoomlength(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function SelectBehOutcome_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.ChosenOutcome = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolOutcome
            P.PlotParamsBehavior.ChosenOutcome = P.PlotParamsBehavior.ChosenOutcome(:);
        else
            P.PlotParamsBehavior.ChosenOutcome = P.PlotParamsBehavior.ChosenOutcome(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function FPoolBehContrast_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.FpoolContrast = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolContrast
            P.PlotParamsBehavior.ChosenContrast = P.PlotParamsBehavior.ChosenContrast(:);
        else
            P.PlotParamsBehavior.ChosenContrast = P.PlotParamsBehavior.ChosenContrast(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function FPoolBehGain_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.FpoolGain = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolGain
            P.PlotParamsBehavior.ChosenGain = P.PlotParamsBehavior.ChosenGain(:);
        else
            P.PlotParamsBehavior.ChosenGain = P.PlotParamsBehavior.ChosenGain(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function FPoolBehRoomlength_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.FpoolRoomlength = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolRoomlength
            P.PlotParamsBehavior.ChosenRoomlength = P.PlotParamsBehavior.ChosenRoomlength(:);
        else
            P.PlotParamsBehavior.ChosenRoomlength = P.PlotParamsBehavior.ChosenRoomlength(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

    function FPoolBehOutcome_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsBehavior.FpoolOutcome = get(hObject,'Value');
        if P.PlotParamsBehavior.FpoolOutcome
            P.PlotParamsBehavior.ChosenOutcome = P.PlotParamsBehavior.ChosenOutcome(:);
        else
            P.PlotParamsBehavior.ChosenOutcome = P.PlotParamsBehavior.ChosenOutcome(:)';
        end
        P.PlotParamsBehavior = UpdateplotBehavior(P.PlotParamsBehavior,EXP,GUI);
    end

end