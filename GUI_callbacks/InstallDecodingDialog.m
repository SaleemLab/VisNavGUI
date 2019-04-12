function InstallDecodingDialog(GUI, EXP, P, DecodingDialog)

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

DecodingDialog.reset;
% DecodingDialog.getlistbox('Probe',{'CA1','V1'}, [0.5 1.5],'Min',1,'Max',2,'Callback',@(hObject, eventdata)GetChosenProbe_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getvariable('Tolerance for Max.', P.PlotParamsDecoder.maxtolerance, [0.5 0.5],'Callback', @(hObject, eventdata)Getmaxtolerance_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('use Xdec distri', P.PlotParamsDecoder.FdecXdistribution, 0.5, 'Callback', @(hObject, eventdata)FdecXdistribution_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('only good time bins', P.PlotParamsDecoder.Fgoodtimebins, 0.5, 'Callback', @(hObject, eventdata)Fgoodtimebins_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('Smooth spatiallly', P.PlotParamsDecoder.Fspatialsmooth, 0.5, 'Callback', @(hObject, eventdata)Fspatialsmooth_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('Use Theta Post.', P.PlotParamsDecoder.FthetaPost, 0.5, 'Callback', @(hObject, eventdata)FthetaPost_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getvariable('Theta Post Ref.', P.PlotParamsDecoder.thetaDecphase, [0.5 0.5], 'Callback', @(hObject, eventdata)GetthetaDecphase_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.addslider([0 360],[1 30],P.PlotParamsDecoder.thetaDecphase, 0.5);
DecodingDialog.getcheckbox('Use Theta bins', P.PlotParamsDecoder.Fthetabins, 0.5, 'Callback', @(hObject, eventdata)Fthetabins_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getvariable('Theta phase', P.PlotParamsDecoder.thetaphase, [0.5 0.5], 'Callback', @(hObject, eventdata)Getdecthetaphase_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.addslider([0 360],[1 30],P.PlotParamsDecoder.thetaphase, 0.5);
DecodingDialog.getvariable('# Theta phase bins', P.PlotParamsDecoder.nthetaphsbins, [0.5 0.5], 'Callback', @(hObject, eventdata)Getdecnthetaphsbins_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.addslider([1 24],[1 2],P.PlotParamsDecoder.nthetaphsbins, 0.5);
DecodingDialog.getvariable('Theta Channel', P.PlotParamsDecoder.thetaChannel, [0.5 0.5],'Callback', @(hObject, eventdata)GetdecthetaChannel_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.addslider([1 37],[1 4],P.PlotParamsDecoder.thetaChannel, 0.5);
DecodingDialog.getvariable('# Theta X bins', P.PlotParamsDecoder.nthetaXbins, [0.5 0.5], 'Callback', @(hObject, eventdata)GetdecnthetaXbins_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.addslider([1 100],[1 2],P.PlotParamsDecoder.nthetaXbins, 0.5);
DecodingDialog.getcheckbox('Use Speed bins', P.PlotParamsDecoder.Fspdbins, 0.5, 'Callback', @(hObject, eventdata)Fspdbins_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getvariable('Speed bin', P.PlotParamsDecoder.spdbin, [0.5 0.5], 'Callback', @(hObject, eventdata)Getspdbin_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.addslider([1 10],[1 2],P.PlotParamsDecoder.spdbin, 0.5);
DecodingDialog.getlistbox('Plots',P.PlotParamsDecoder.PlotObjList, [0.5 3], 'Min',1,'Max',numel(P.PlotParamsDecoder.PlotObjList), 'Callback', @(hObject, eventdata)SelectDecPlotObj_callback(hObject, eventdata, GUI, EXP, P));
if Fcontrast
    DecodingDialog.getlistbox('Contrast',strsplit(num2str(EXP.SubsetVal.contrast)), [0.5 2],'Min', 1,'Max', numel(EXP.SubsetVal.contrast), 'Callback', @(hObject, eventdata)SelectDecContrast_callback(hObject, eventdata, GUI, EXP, P));
    DecodingDialog.getcheckbox('Pool Contrasts', P.PlotParamsDecoder.FpoolContrast, 0.5, 'Callback', @(hObject, eventdata)FPoolDecContrast_callback(hObject, eventdata, GUI, EXP, P));
end
if Fgain
    DecodingDialog.getlistbox('Gain',strsplit(num2str(EXP.SubsetVal.gain)), [0.5 2], 'Min', 1, 'Max', numel(EXP.SubsetVal.gain), 'Callback', @(hObject, eventdata)SelectDecGain_callback(hObject, eventdata, GUI, EXP, P));
    DecodingDialog.getcheckbox('Pool Gain', P.PlotParamsDecoder.FpoolGain, 0.5, 'Callback', @(hObject, eventdata)FPoolDecGain_callback(hObject, eventdata, GUI, EXP, P));
end
if FroomLength
    DecodingDialog.getlistbox('Roomlength',strsplit(num2str(EXP.SubsetVal.roomlength)), [0.5 1], 'Min', 1, 'Max', numel(EXP.SubsetVal.roomlength), 'Callback', @(hObject, eventdata)SelectDecRoomlength_callback(hObject, eventdata, GUI, EXP, P));
    DecodingDialog.getcheckbox('Pool Roomlength', P.PlotParamsDecoder.FpoolRoomlength, 0.5, 'Callback', @(hObject, eventdata)FPoolDecRoomlength_callback(hObject, eventdata, GUI, EXP, P));
end
if Foutcome
    DecodingDialog.getlistbox('Outcome',strsplit(num2str(EXP.SubsetVal.outcome)),[0.5 2], 'Min', 1,'Max', numel(EXP.SubsetVal.outcome), 'Callback', @(hObject, eventdata)SelectDecOutcome_callback(hObject, eventdata, GUI, EXP, P));
    DecodingDialog.getcheckbox('Pool Outcome', P.PlotParamsDecoder.FpoolOutcome, 0.5, 'Callback', @(hObject, eventdata)FPoolDecOutcome_callback(hObject, eventdata, GUI, EXP, P));
end
DecodingDialog.getcheckbox('Display Log', P.PlotParamsDecoder.Fdisplaylog, 0.5, 'Callback', @(hObject, eventdata)Fdisplaylog_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getvariable('Climits', P.PlotParamsDecoder.Clim, [0.5 0.5], 'Callback', @(hObject, eventdata)GetClim_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('Normalize', P.PlotParamsDecoder.Fnormalize, 0.5, 'Callback', @(hObject, eventdata)Fnormalize_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getpopupmenu('Palette',{'RedWhiteBlue','parula','hot','hotcoldlin'}, 1, [0.5 0.7], 'Callback', @(hObject, eventdata)GetPalettename_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('overlap', P.PlotParamsDecoder.Foverlap, 0.5, 'Callback', @(hObject, eventdata)Foverlap_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('display PredMax', P.PlotParamsDecoder.FdisplayPredMax, 0.5, 'Callback', @(hObject, eventdata)FdisplayPredMax_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('display PredAve', P.PlotParamsDecoder.FdisplayPredAve, 0.5, 'Callback', @(hObject, eventdata)FdisplayPredAve_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('display Mat', P.PlotParamsDecoder.FdisplayMat, 0.5, 'Callback', @(hObject, eventdata)FdisplayMat_callback(hObject, eventdata, GUI, EXP, P));
DecodingDialog.getcheckbox('display LFP', P.PlotParamsDecoder.FdisplayLFP, 0.5, 'Callback', @(hObject, eventdata)FdisplayLFP_callback(hObject, eventdata, GUI, EXP, P));

%Callback functions
    function GetChosenProbe_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.ChosenProbe = get(hObject,'Val');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Getmaxtolerance_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.maxtolerance = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FdecXdistribution_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FdecXdistribution = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end


    function Fgoodtimebins_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Fgoodtimebins = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Fspatialsmooth_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Fspatialsmooth = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end


    function FthetaPost_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FthetaPost = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Fthetabins_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Fthetabins = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Fspdbins_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Fspdbins = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function GetdecthetaChannel_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.thetaChannel = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function GetthetaDecphase_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.thetaDecphase = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Getspdbin_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.spdbin = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Getdecthetaphase_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.thetaphase = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Getdecnthetaphsbins_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.nthetaphsbins = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function GetdecnthetaXbins_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.nthetaXbins = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecX_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsDecoder.ChosenVarX = str(val);
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecY_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsDecoder.ChosenVarY = str(val);
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecPlotObj_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = get(hObject,'String');
        P.PlotParamsDecoder.ChosenObj = str(val);
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecContrast_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.ChosenContrast = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolContrast
            P.PlotParamsDecoder.ChosenContrast = P.PlotParamsDecoder.ChosenContrast(:);
        else
            P.PlotParamsDecoder.ChosenContrast = P.PlotParamsDecoder.ChosenContrast(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecGain_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.ChosenGain = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolGain
            P.PlotParamsDecoder.ChosenGain = P.PlotParamsDecoder.ChosenGain(:);
        else
            P.PlotParamsDecoder.ChosenGain = P.PlotParamsDecoder.ChosenGain(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecRoomlength_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.ChosenRoomlength = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolRoomlength
            P.PlotParamsDecoder.ChosenRoomlength = P.PlotParamsDecoder.ChosenRoomlength(:);
        else
            P.PlotParamsDecoder.ChosenRoomlength = P.PlotParamsDecoder.ChosenRoomlength(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function SelectDecOutcome_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.ChosenOutcome = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolOutcome
            P.PlotParamsDecoder.ChosenOutcome = P.PlotParamsDecoder.ChosenOutcome(:);
        else
            P.PlotParamsDecoder.ChosenOutcome = P.PlotParamsDecoder.ChosenOutcome(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end


    function FPoolDecContrast_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FpoolContrast = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolContrast
            P.PlotParamsDecoder.ChosenContrast = P.PlotParamsDecoder.ChosenContrast(:);
        else
            P.PlotParamsDecoder.ChosenContrast = P.PlotParamsDecoder.ChosenContrast(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FPoolDecGain_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FpoolGain = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolGain
            P.PlotParamsDecoder.ChosenGain = P.PlotParamsDecoder.ChosenGain(:);
        else
            P.PlotParamsDecoder.ChosenGain = P.PlotParamsDecoder.ChosenGain(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FPoolDecRoomlength_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FpoolRoomlength = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolRoomlength
            P.PlotParamsDecoder.ChosenRoomlength = P.PlotParamsDecoder.ChosenRoomlength(:);
        else
            P.PlotParamsDecoder.ChosenRoomlength = P.PlotParamsDecoder.ChosenRoomlength(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FPoolDecOutcome_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FpoolOutcome = get(hObject,'Value');
        if P.PlotParamsDecoder.FpoolOutcome
            P.PlotParamsDecoder.ChosenOutcome = P.PlotParamsDecoder.ChosenOutcome(:);
        else
            P.PlotParamsDecoder.ChosenOutcome = P.PlotParamsDecoder.ChosenOutcome(:)';
        end
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Fdisplaylog_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Fdisplaylog = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function GetClim_callback(hObject, eventdata, GUI, EXP, P)
        val = str2num(get(hObject,'String'));
        P.PlotParamsDecoder.Clim = val;
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Fnormalize_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Fnormalize = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function GetPalettename_callback(hObject, eventdata, GUI, EXP, P)
        val = get(hObject,'Value');
        str = (get(hObject,'String'));
        P.PlotParamsDecoder.Palettename = str{val};
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function Foverlap_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.Foverlap = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FdisplayPredMax_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FdisplayPredMax = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FdisplayPredAve_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FdisplayPredAve = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FdisplayMat_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FdisplayMat = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

    function FdisplayLFP_callback(hObject, eventdata, GUI, EXP, P)
        P.PlotParamsDecoder.FdisplayLFP = get(hObject,'Value');
        UpdateplotDecodSpks(P.PlotParamsDecoder,EXP,GUI);
    end

end