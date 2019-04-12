function obj = defineCellProp(obj,Nperm)

nbcont = numel(obj.SubsetVal.contrast) + 1;
nbgain = numel(obj.SubsetVal.gain);
nbroomlength = numel(obj.SubsetVal.roomlength);
nboutcome = numel(obj.SubsetVal.outcome);

cbase = numel(obj.SubsetVal.contrast) + 1; %find(obj.SubsetVal.contrast == mode(obj.data.es.contrast));
gbase = find(obj.SubsetVal.gain == mode(obj.data.es.gain));
rbase = find(obj.SubsetVal.roomlength == mode(obj.data.es.roomLength));
obase = find(obj.SubsetVal.outcome == 2);
if isfield(obj.maps1d,'trajPercent')
    numbins = obj.maps1d.trajPercent{cbase, gbase, rbase, obase}.numBins;
else
    numbins = 100;
end

obj.CellInfo.Finterneuron = [];
obj.CellInfo.min2maxSpkwf = [];
obj.CellInfo.minwf = [];
obj.CellInfo.max1wf = [];
obj.CellInfo.max2wf = [];
obj.CellInfo.asymmetrywf = [];
obj.CellInfo.spkautocorr = [];
obj.CellInfo.burstindex = [];
obj.CellInfo.rateautocorr = [];
obj.CellInfo.Spatialmodulation = [];
obj.CellInfo.field = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field_half1 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field_half2 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field_set1 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field_set2 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldZ = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldXcorr = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldCOM = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldCOMSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsize = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsizeSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.rate = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldMax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldMin = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsinOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsinOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldsinOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.fieldShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfCOM = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfCOMSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsinOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsinOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsinOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsize = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfsizeSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldCOM_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldPos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.reliabilityCorr = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.fieldShift = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShiftSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfShift = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldShfShiftiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.rategain = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.rategainSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.rateShfgain = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.rateShfgainiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

%
obj.CellInfo.unwrapfield = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfield_half1 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfield_half2 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfield_set1 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfield_set2 = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldZ = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldXcorr = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldsize = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldsizeSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldMax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldMin = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldhaldiff = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldhaldiffSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.unwrapfieldShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfsize = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfsizeSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfhaldiff = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldShfhaldiffiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.unwrapfield_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapfieldPos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.unwrapreliabilityCorr = cell(nbcont, nbgain, nbroomlength, nboutcome);
%

obj.CellInfo.phsfieldPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldCOM = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldCOMSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfield = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsmodulation = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldsinOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldsinOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldsinOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.phsfieldShfPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfCOM = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfCOMSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsmodulationShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfsinOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfsinOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldShfsinOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.phsfield_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldCOM_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.phsfieldPos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

%
obj.CellInfo.spkphsfieldPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldCOM = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldCOMSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfield = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsmodulation = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldsinOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldsinOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldsinOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkphsfieldShfPos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfPosSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfCOM = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfCOMSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsmodulationShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfsinOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfsinOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldShfsinOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkphsfield_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldCOM_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkphsfieldPos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
%

obj.CellInfo.field2dXspd = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXspdSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dX = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dSpd = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dSpdSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dSpdpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dSpdposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dSpdmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dSpdmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

%
obj.CellInfo.field2dXPhstheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXPhsthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dPhstheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dslopeXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dphi0XY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2drhoXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dslopeXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dphi0XYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2drhoXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dslopeXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dphi0XYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2drhoXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.Phase = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXcorrtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXRefcorrtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetamaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dPhscorrtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXPhsthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXPhsthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhsthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dShfslopeXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfphi0XY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfrhoXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfslopeXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfphi0XYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfrhoXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfslopeXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfphi0XYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dShfrhoXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXcorrthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetaShfmaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXRefcorrthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXRefcorrthetaShfmaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dPhscorrthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXcorrtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXcorrthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dPhscorrtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dPhscorrthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dXPhstheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dXthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.field2dYtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dYthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.field2dYthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);


obj.CellInfo.PhaseRho = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.PhasePval = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.PhaseMax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.PhaseRayleighPval = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.PhaseRayleighZ = cell(nbcont, nbgain, nbroomlength, nboutcome);
%

obj.CellInfo.spkfield2dXPhstheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXPhsthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dPhstheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dslopeXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dphi0XY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2drhoXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dslopeXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dphi0XYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2drhoXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dslopeXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dphi0XYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2drhoXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkPhase = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXcorrtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXRefcorrtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetamaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dPhscorrtheta = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetapos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetamax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetamaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXPhsthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXPhsthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhsthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dShfslopeXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfphi0XY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfrhoXY = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfslopeXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfphi0XYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfrhoXYSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfslopeXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfphi0XYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dShfrhoXYiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXcorrthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetaShfmaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXRefcorrthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmp = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmpSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmpiter = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffset = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffsetSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffsetiter = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dPhscorrthetaShf = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaShfSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaShfpos = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaShfposSE = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaShfmax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetaShfmaxSE = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXcorrtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXcorrthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dPhscorrtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dPhscorrthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dXPhstheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dXthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.spkfield2dYtheta_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dYthetapos_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkfield2dYthetamax_2fold = cell(nbcont, nbgain, nbroomlength, nboutcome);


obj.CellInfo.spkPhaseRho = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkPhasePval = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkPhaseMax = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkPhaseRayleighPval = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.spkPhaseRayleighZ = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.phsfieldMUA = cell(1,2);
obj.CellInfo.phsfieldMUASE = cell(1,2);
obj.CellInfo.phsfieldZMUA = cell(1,2);
obj.CellInfo.phsfieldPosMUA = cell(1,2);
obj.CellInfo.LFP2Spike_phscorrMUA = cell(1,2);

obj.CellInfo.phsfieldMUAnorm = cell(1,2);
obj.CellInfo.phsfieldPosMUAnorm = cell(1,2);
obj.CellInfo.LFP2Spike_phscorrMUAnorm = cell(1,2);

obj.CellInfo.phsfieldMean = cell(1,2);
obj.CellInfo.phsfieldPosMean = cell(1,2);
obj.CellInfo.LFP2Spike_phscorrMean = cell(1,2);

obj.CellInfo.SSI = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SSIRef = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SpatialInfo = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SpatialInfoPerSpike = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SpatialInfoRef = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SpatialInfoPerSpikeRef = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldRef = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldStdRef = cell(nbcont, nbgain, nbroomlength, nboutcome);

obj.CellInfo.SSI_dec = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SpatialInfo_dec = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.SpatialInfoPerSpike_dec = cell(nbcont, nbgain, nbroomlength, nboutcome);
obj.CellInfo.fieldRef_dec = cell(nbcont, nbgain, nbroomlength, nboutcome);

if ~isfield(obj.CellInfo,'trajdecCell_Xave')
    obj.CellInfo.trajdecCell_Xave = cell(1,obj.CellInfo.NumCells);
    obj.CellInfo.trajdecCell_Xmax = cell(1,obj.CellInfo.NumCells);
    obj.CellInfo.trajdecAll_Xave = cell(1,obj.CellInfo.NumCells);
    obj.CellInfo.trajdecAll_Xmax = cell(1,obj.CellInfo.NumCells);
end

roomlength = round(max(obj.data.es.trajPercent)/50)*50;

if isfield(obj.data.es,'bestchan')
    obj.CellInfo.bestchan = cell2mat(obj.data.es.bestchan);
else
    obj.CellInfo.bestchan = NaN(1,obj.CellInfo.NumCells);
end
if isfield(obj.data.es,'waveform')
    obj.CellInfo.waveform = cell2mat(obj.data.es.waveform(:));
else
    nwpts = 61;
    obj.CellInfo.waveform = NaN(obj.CellInfo.NumCells,nwpts);
end
obj.CellInfo.contrastVal = obj.SubsetVal.contrast;
obj.CellInfo.gainVal = obj.SubsetVal.gain;
obj.CellInfo.outcomeVal = obj.SubsetVal.outcome;

obj.CellInfo.Nperm = Nperm;

obj.CellInfo.RFXpos = [];
obj.CellInfo.Xpos = [];
obj.CellInfo.XposZmax = [];
obj.CellInfo.globalXposrep = [];
obj.CellInfo.RFYpos = [];
obj.CellInfo.Ypos = [];
obj.CellInfo.YposZmax = [];
obj.CellInfo.globalYposrep = [];
if isfield(obj.data,'VStuning')
    if ~isempty(obj.data.VStuning)
        for icell = 1:obj.CellInfo.NumCells
            for iexp = 1:numel(obj.data.VStuning)
                switch obj.data.VStuning(iexp).name
                    case '10 x1'
                        obj.CellInfo.Xpos(icell) = mean(obj.data.VStuning(iexp).Xpos(obj.data.VStuning(iexp).mean(icell,1:end-1) == max(obj.data.VStuning(iexp).mean(icell,1:end-1))));
                        obj.CellInfo.RFXpos(icell,:) = obj.data.VStuning(iexp).mean(icell,1:end-1);
                        obj.CellInfo.XposZmax(icell) = max(obj.data.VStuning(iexp).zscore(icell,1:end-1));
                        obj.CellInfo.globalXposrep(icell,:) = obj.data.VStuning(iexp).globalrep(icell,:);
                    case '9 y1'
                        obj.CellInfo.Ypos(icell) = mean(obj.data.VStuning(iexp).Ypos(obj.data.VStuning(iexp).mean(icell,1:end-1) == max(obj.data.VStuning(iexp).mean(icell,1:end-1))));
                        obj.CellInfo.RFYpos(icell,:) = obj.data.VStuning(iexp).mean(icell,1:end-1);
                        obj.CellInfo.YposZmax(icell) = max(obj.data.VStuning(iexp).zscore(icell,1:end-1));
                        obj.CellInfo.globalYposrep(icell,:) = obj.data.VStuning(iexp).globalrep(icell,:);
                end
            end
        end
    end
end

for c = 1:numel(obj.SubsetVal.contrast) +1
    for g = 1:numel(obj.SubsetVal.gain)
        for r = 1:numel(obj.SubsetVal.roomlength)
            for o = 1:numel(obj.SubsetVal.outcome)
                if c > numel(obj.SubsetVal.contrast)
                    contidx = find(obj.SubsetVal.contrast>0);
                else
                    contidx = c;
                end
                tidx = obj.getSubsets(contidx, g, r, o);
                samplerate = 1./obj.data.es.sampleSize;
                dx = 1;
                nbXbinsmth = (numbins/1)/4;
                if (c == cbase || c > numel(obj.SubsetVal.contrast)) && g == gbase ...
                        && r == rbase && o == obase
                    usedNperm = Nperm;
                else
                    usedNperm = 0;
                end
                for iperm = 0:usedNperm
                    idx = find(tidx);
                    for icell = 1:obj.CellInfo.NumCells
                        meanrate = (sum(obj.data.es.spikeTrain(tidx,icell))/sum(tidx))*mean(samplerate);
                        if iperm == 0
                            spktrain = obj.data.es.spikeTrain(idx,icell);
                        else
                            try
                            spkISI = [1;diff(find(obj.data.es.spikeTrain(idx,icell) > 0))];
                            spktrain = 0*obj.data.es.spikeTrain(idx,icell);
                            if ~isempty(spktrain) && numel(spkISI) > 1
                                spktrain(cumsum(spkISI(randperm(numel(spkISI))))) = 1;
                            end
                            catch
                                keyboard
                            end
                        end
                        
                        [SInfo, SInfoperSpk, map] = SpatialInformation(obj.data.es.trajPercent(tidx), spktrain(:,:), dx, samplerate(tidx),nbXbinsmth,obj.data.es.CircularMaze);
                        SSI = (max(map) - mean(map))/mean(map);
                        map = map(:)';
                        if ~isempty(map)
                            if iperm == 0
                                try
                                    obj.CellInfo.SSI{c, g, r, o}(icell) = SSI;
                                    obj.CellInfo.SpatialInfo{c, g, r, o}(icell) = SInfo;
                                    obj.CellInfo.SpatialInfoPerSpike{c, g, r, o}(icell) = SInfoperSpk;
                                    obj.CellInfo.fieldRef{c, g, r, o}(icell,:) = 0*map;
                                    obj.CellInfo.fieldStdRef{c, g, r, o}(icell,:) = 0*map;
                                catch
                                    keyboard
                                end
                            else
                                try
                                    obj.CellInfo.SSIRef{c, g, r, o}(icell,iperm) = SSI;
                                    obj.CellInfo.SpatialInfoRef{c, g, r, o}(icell,iperm) = SInfo;
                                    obj.CellInfo.SpatialInfoPerSpikeRef{c, g, r, o}(icell,iperm) = SInfoperSpk;
                                    obj.CellInfo.fieldRef{c, g, r, o}(icell,:) = obj.CellInfo.fieldRef{c, g, r, o}(icell,:) + map/Nperm;
                                    obj.CellInfo.fieldStdRef{c, g, r, o}(icell,:) = obj.CellInfo.fieldStdRef{c, g, r, o}(icell,:) + map.^2/Nperm;
                                catch
                                    keyboard
                                end
                            end
                        else
                            obj.CellInfo.SpatialInfo{c, g, r, o}(icell) = NaN;
                            obj.CellInfo.SpatialInfoPerSpike{c, g, r, o}(icell) = NaN;
                            obj.CellInfo.fieldRef{c, g, r, o}(icell,:) = 0;
                            obj.CellInfo.fieldStdRef{c, g, r, o}(icell,:) = 0;
                        end
                        try
                            if icell <= numel(obj.CellInfo.trajdecCell_Xave)
                                if ~isempty(obj.CellInfo.trajdecCell_Xave{icell})
                                    [SInfo_dec, SInfoperSpk_dec, map_dec] = SpatialInformation(obj.CellInfo.trajdecCell_Xave{icell}(tidx)-1, spktrain(:,:), dx, samplerate(tidx),nbXbinsmth,obj.data.es.CircularMaze);
                                    SSI_dec = (max(map_dec) - mean(map_dec))/mean(map_dec);
                                    map_dec = map_dec(:)';
                                    if ~isempty(map_dec)
                                        if iperm == 0
                                            try
                                                obj.CellInfo.SSI_dec{c, g, r, o}(icell) = SSI_dec;
                                                obj.CellInfo.SpatialInfo_dec{c, g, r, o}(icell) = SInfo_dec;
                                                obj.CellInfo.SpatialInfoPerSpike_dec{c, g, r, o}(icell) = SInfoperSpk_dec;
                                                obj.CellInfo.field_dec{c, g, r, o}(icell,:) = 0*map_dec;
                                            catch
                                                keyboard
                                            end
                                        end
                                    else
                                        obj.CellInfo.SSI_dec{c, g, r, o}(icell) = NaN;
                                        obj.CellInfo.SpatialInfo_dec{c, g, r, o}(icell) = NaN;
                                        obj.CellInfo.SpatialInfoPerSpike_dec{c, g, r, o}(icell) = NaN;
                                        obj.CellInfo.fieldRef_dec{c, g, r, o}(icell,:) = 0;
                                    end
                                end
                            end
                        catch
                            keyboard
                        end
                    end
                end
                if Nperm > 0
                    obj.CellInfo.fieldStdRef{c, g, r, o} = (Nperm/(Nperm-1)*(obj.CellInfo.fieldStdRef{c, g, r, o} - obj.CellInfo.fieldRef{c, g, r, o}.^2)).^0.5;
                else
                    obj.CellInfo.fieldStdRef{c, g, r, o} = NaN(size(obj.CellInfo.fieldStdRef{c, g, r, o}));
                end
                
                idxall = find(tidx);
                randidxset = randperm(numel(idxall));
                tidx_set1 = false(size(tidx));
                tidx_set2 = false(size(tidx));
                tidx_set1(idxall(randidxset(1:floor(numel(randidxset)/2)))) = true;
                tidx_set2(idxall(randidxset(floor(numel(randidxset)/2)+1:end))) = true;
                tidx_set1 = tidx_set1 & obj.data.es.completetrial;
                tidx_set2 = tidx_set2 & obj.data.es.completetrial;
                
                for icell = 1:obj.CellInfo.NumCells
                    map = [];
                    if isfield(obj.maps1d,'trajPercent')
                        if ~isempty(obj.maps1d.trajPercent{c, g, r, o}.model)
                            map = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).meanrespModel;
                            mapbase = obj.maps1d.trajPercent{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel;
                            mapSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SErespModel;
                            fieldCOM = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            fieldCOMSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            fieldmax = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            fieldmaxSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            fieldAmp = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            fieldAmpSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            fieldAmpiter = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            fieldsinAmp = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).meanrespModelXsinAmp;
                            fieldsinAmpSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SErespModelXsinAmp;
                            fieldsinAmpiter = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).respModelXsinAmp;
                            fieldsinOffset = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).meanrespModelXsinOffset;
                            fieldsinOffsetSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SErespModelXsinOffset;
                            fieldsinOffsetiter = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).respModelXsinOffset;
                            fieldsize = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).fieldsize;
                            fieldsizeSE = obj.maps1d.trajPercent{c, g, r, o}.model.tuning(icell).SEfieldsize;
                            x = obj.maps1d.trajPercent{c, g, r, o}.bins;
                        end
                    end
                    map_Shf = [];
                    if isfield(obj.maps1d,'trajPercent_Shf')
                        if ~isempty(obj.maps1d.trajPercent_Shf{c, g, r, o}.model)
                            map_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).meanrespModel;
                            mapSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SErespModel;
                            fieldCOM_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            fieldCOMSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            fieldmax_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            fieldmaxSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            fieldAmp_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            fieldAmpSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            fieldAmpiter_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            fieldsinAmp_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXsinAmp;
                            fieldsinAmpSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SErespModelXsinAmp;
                            fieldsinAmpiter_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).respModelXsinAmp;
                            fieldsinOffset_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXsinOffset;
                            fieldsinOffsetSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SErespModelXsinOffset;
                            fieldsinOffsetiter_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).respModelXsinOffset;
                            fieldsize_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).fieldsize;
                            fieldsizeSE_Shf = obj.maps1d.trajPercent_Shf{c, g, r, o}.model.tuning(icell).SEfieldsize;
                        end
                    end
                    map_2fold = [];
                    if isfield(obj.maps1d,'trajPercent_2fold')
                        if ~isempty(obj.maps1d.trajPercent_2fold{c, g, r, o}.model)
                            map_2fold = obj.maps1d.trajPercent_2fold{c, g, r, o}.model.tuning(icell).respModel;
                            fieldCOM_2fold = obj.maps1d.trajPercent_2fold{c, g, r, o}.model.tuning(icell).respModelXpos;
                            fieldmax_2fold = obj.maps1d.trajPercent_2fold{c, g, r, o}.model.tuning(icell).respModelXmax;
                        end
                    end
                    if c == cbase && g == gbase && r == rbase && o == obase
                        %normaveform = (waveform)./repmat(max(abs(waveform),[],1),[size(waveform,1) 1]);
                        Spkwfcenter = floor(size(obj.CellInfo.waveform,2)/2);
%                         obj.CellInfo.min2maxSpkwf(icell) = (Spkwfcenter+find(obj.CellInfo.waveform(icell,(Spkwfcenter+1):end)==max(obj.CellInfo.waveform(icell,(Spkwfcenter+1):end)),1,'first'))-find(obj.CellInfo.waveform(icell,:)==min(obj.CellInfo.waveform(icell,:)),1,'first');
                        
                        [pks,minidx] = findpeaks(-obj.CellInfo.waveform(icell,:));%find(S.CellInfo.waveform(icell,:) == min(S.CellInfo.waveform(icell,:)));
                        pks = pks(minidx >= 10 & minidx <= size(obj.CellInfo.waveform,2)-10);
                        minidx = minidx(minidx >= 10 & minidx <= size(obj.CellInfo.waveform,2)-10);
                        minidx = minidx(pks == max(pks));
                        if minidx >= 10 && minidx <= size(obj.CellInfo.waveform,2)-10
                            minidx = minidx(1);
                            mini =  min(obj.CellInfo.waveform(icell,minidx));
                            obj.CellInfo.minwf(icell) = mini(1);
                            maxi = max(obj.CellInfo.waveform(icell,1:(minidx-1)));
                            obj.CellInfo.max1wf(icell) = maxi(1);
                            [maxi,maxidx] = max(obj.CellInfo.waveform(icell,(minidx+1):end));
                            obj.CellInfo.max2wf(icell) = maxi(1);
                            obj.CellInfo.asymmetrywf(icell) = (obj.CellInfo.max2wf(icell) - obj.CellInfo.max1wf(icell))/(obj.CellInfo.max2wf(icell) + obj.CellInfo.max2wf(icell));
                            obj.CellInfo.min2maxSpkwf(icell) = maxidx;
                        else
                            obj.CellInfo.minwf(icell) = NaN;
                            obj.CellInfo.max1wf(icell) = NaN;
                            obj.CellInfo.max2wf(icell) = NaN;
                            obj.CellInfo.asymmetrywf(icell) = NaN;
                            obj.CellInfo.min2maxSpkwf(icell) = NaN;
                        end
                        
                        spkautocorr = 60*xcorr(obj.data.es.spikeTrain(:,icell),60,'coeff');
                        spkautocorr = spkautocorr(62:end);
                        obj.CellInfo.spkautocorr(icell,:) = spkautocorr;
                        obj.CellInfo.burstindex(icell) = max(spkautocorr(1:3))/max(spkautocorr(6:end));
                        obj.CellInfo.rateautocorr(icell) = nanmean(spkautocorr(end-5:end));
                        
                        [minFRmapbase,~] = fast1Dmap(obj.data.es.trajPercent(tidx), obj.data.es.spikeTrain(idx,icell), dx, samplerate(tidx), nbXbinsmth, obj.data.es.CircularMaze);
                        minFR = min(minFRmapbase);
                        obj.CellInfo.minFR(icell) = minFR;
                        if ~isempty(minFRmapbase)
                            obj.CellInfo.Spatialmodulation(icell) = (max(minFRmapbase) - mean(minFRmapbase))/mean(minFRmapbase);
                        else
                            obj.CellInfo.Spatialmodulation(icell) = NaN;
                        end
                        if obj.CellInfo.min2maxSpkwf(icell) < 18 && obj.CellInfo.Spatialmodulation(icell) < 2 ... %(((obj.CellInfo.min2maxSpkwf(icell) < 18 || obj.CellInfo.burstindex(icell) < 1) && (obj.CellInfo.minFR(icell) > 4 && obj.CellInfo.Spatialmodulation(icell) < 0.5)) || (obj.CellInfo.minFR(icell) > 15)) ... %obj.CellInfo.min2maxSpkwf(icell) <= 15 ... %
                           && obj.CellInfo.Probe(icell) == 1
                            obj.CellInfo.Finterneuron(icell) = true;
                            if ~contains(obj.CellInfo.CellListString{icell},'IN')
                                obj.CellInfo.CellListString{icell} = [obj.CellInfo.CellListString{icell} '_IN'];
                            end
                        else
                            stridx = strfind(obj.CellInfo.CellListString{icell},'IN');
                            if ~isempty(stridx)
                                obj.CellInfo.CellListString{icell} = obj.CellInfo.CellListString{icell}(1:stridx(1)-2);
                            end
                            obj.CellInfo.Finterneuron(icell) = false;
                        end
                    end
                    
                    if ~isempty(map)
                        obj.CellInfo.fieldPos{c, g, r, o}(icell) = fieldmax;
                        obj.CellInfo.fieldPosSE{c, g, r, o}(icell) = fieldmaxSE;
                        obj.CellInfo.fieldCOM{c, g, r, o}(icell) = fieldCOM;
                        obj.CellInfo.fieldCOMSE{c, g, r, o}(icell) = fieldCOMSE;
                        obj.CellInfo.field{c, g, r, o}(icell,:) = map;
                        obj.CellInfo.fieldSE{c, g, r, o}(icell,:) = mapSE;
                        obj.CellInfo.fieldAmp{c, g, r, o}(icell) = fieldAmp;
                        obj.CellInfo.fieldAmpSE{c, g, r, o}(icell) = fieldAmpSE;
                        obj.CellInfo.fieldAmpiter{c, g, r, o}(icell,:) = fieldAmpiter;
                        obj.CellInfo.fieldsinAmp{c, g, r, o}(icell) = fieldsinAmp;
                        obj.CellInfo.fieldsinAmpSE{c, g, r, o}(icell) = fieldsinAmpSE;
                        obj.CellInfo.fieldsinAmpiter{c, g, r, o}(icell,:) = fieldsinAmpiter;
                        obj.CellInfo.fieldsinOffset{c, g, r, o}(icell) = fieldsinOffset;
                        obj.CellInfo.fieldsinOffsetSE{c, g, r, o}(icell) = fieldsinOffsetSE;
                        obj.CellInfo.fieldsinOffsetiter{c, g, r, o}(icell,:) = fieldsinOffsetiter;
                        
                        obj.CellInfo.fieldsize{c, g, r, o}(icell) = fieldsize;
                        obj.CellInfo.fieldsizeSE{c, g, r, o}(icell) = fieldsizeSE;
                        mapcorr = map - nanmean(map);
                        mapcorr = mapcorr./sqrt(sum(mapcorr.^2));
                        mapbasecorr = mapbase - nanmean(mapbase);
                        mapbasecorr = mapbasecorr./sqrt(sum(mapbasecorr.^2));
                        fieldXcorr = zeros(1,numel(mapbase));
                        xshiftlim = floor(numel(mapbase)/2);
                        ishift = 0;
                        for xshift = -xshiftlim:xshiftlim-1
                            ishift = ishift + 1;
                            fieldXcorr(ishift) = sum(mapcorr.*circshift(mapbasecorr,xshift))/sqrt(sum(mapcorr.^2)*sum(mapbasecorr.^2));
                        end
                        obj.CellInfo.fieldXcorr{c, g, r, o}(icell,:) = fieldXcorr;
                        [~,imax] = max(map);
                        obj.CellInfo.fieldZ{c, g, r, o}(icell) = (map(imax) - mean(map))./mapSE(imax);
                        obj.CellInfo.rate{c, g, r, o}(icell) = (sum(obj.data.es.spikeTrain(tidx,icell))/sum(tidx))*mean(samplerate(tidx));
                        obj.CellInfo.fieldMax{c, g, r, o}(icell) = max(map);
                        obj.CellInfo.fieldMin{c, g, r, o}(icell) = min(map);
                    else
                        obj.CellInfo.fieldPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldCOM{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldCOMSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldsize{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldsizeSE{c, g, r, o}(icell) = NaN;
                        
                        obj.CellInfo.fieldAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldsinOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldsinOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldsinOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.fieldXcorr{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldZ{c, g, r, o}(icell) = 0;
                        obj.CellInfo.rate{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldMax{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldMin{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(map_Shf)
                        obj.CellInfo.fieldShfPos{c, g, r, o}(icell) = fieldmax_Shf;
                        obj.CellInfo.fieldShfPosSE{c, g, r, o}(icell) = fieldmaxSE_Shf;
                        obj.CellInfo.fieldShfCOM{c, g, r, o}(icell) = fieldCOM_Shf;
                        obj.CellInfo.fieldShfCOMSE{c, g, r, o}(icell) = fieldCOMSE_Shf;
                        obj.CellInfo.fieldShf{c, g, r, o}(icell,:) = map_Shf;
                        obj.CellInfo.fieldShfSE{c, g, r, o}(icell,:) = mapSE_Shf;
                        obj.CellInfo.fieldShfsize{c, g, r, o}(icell) = fieldsize_Shf;
                        obj.CellInfo.fieldShfsizeSE{c, g, r, o}(icell) = fieldsizeSE_Shf;
                        
                        obj.CellInfo.fieldShfAmp{c, g, r, o}(icell) = fieldAmp_Shf;
                        obj.CellInfo.fieldShfAmpSE{c, g, r, o}(icell) = fieldAmpSE_Shf;
                        obj.CellInfo.fieldShfAmpiter{c, g, r, o}(icell,:) = fieldAmpiter_Shf;
                        obj.CellInfo.fieldShfsinAmp{c, g, r, o}(icell) = fieldsinAmp_Shf;
                        obj.CellInfo.fieldShfsinAmpSE{c, g, r, o}(icell) = fieldsinAmpSE_Shf;
                        obj.CellInfo.fieldShfsinAmpiter{c, g, r, o}(icell,:) = fieldsinAmpiter_Shf;
                        obj.CellInfo.fieldShfsinOffset{c, g, r, o}(icell) = fieldsinOffset_Shf;
                        obj.CellInfo.fieldShfsinOffsetSE{c, g, r, o}(icell) = fieldsinOffsetSE_Shf;
                        obj.CellInfo.fieldShfsinOffsetiter{c, g, r, o}(icell,:) = fieldsinOffsetiter_Shf;
                    else
                        obj.CellInfo.fieldShfPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShfPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShfCOM{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShfCOMSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldShfsize{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldShfsizeSE{c, g, r, o}(icell) = NaN;
                        
                        obj.CellInfo.fieldShfAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldShfAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShfAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldShfsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldShfsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShfsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.fieldShfsinOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldShfsinOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.fieldShfsinOffsetiter{c, g, r, o}(icell,:) = 0;
                    end
                    if ~isempty(map_2fold)
                        obj.CellInfo.field_2fold{c, g, r, o}(icell,:,:) = map_2fold;
                        obj.CellInfo.fieldCOM_2fold{c, g, r, o}(icell,:) = fieldCOM_2fold;
                        obj.CellInfo.fieldPos_2fold{c, g, r, o}(icell,:) = fieldmax_2fold;
                        fieldset1 = squeeze(map_2fold(1,:));
                        fieldset2 = squeeze(map_2fold(2,:));
                        reliabilityCorr = sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).*...
                            (fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])),2)./...
                            sqrt(sum((fieldset1-repmat(nanmean(fieldset1,2),[1 size(fieldset1,2)])).^2,2).*...
                            sum((fieldset2-repmat(nanmean(fieldset2,2),[1 size(fieldset2,2)])).^2,2));
                        obj.CellInfo.reliabilityCorr{c, g, r, o}(icell) = reliabilityCorr;
                        
                    else
                        obj.CellInfo.field_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.fieldCOM_2fold{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.fieldPos_2fold{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.reliabilityCorr{c, g, r, o}(icell) = NaN;
                    end
                    
                    map1dgain = [];
                    if isfield(obj.maps1d,'gain')
                        if ~isempty(obj.maps1d.gain{c, 2, r, o}.model)
                            map1dgain = obj.maps1d.gain{c, 2, r, o}.model.tuning(icell).meanrespModel(g);
                            map1dgainSE = obj.maps1d.gain{c, 2, r, o}.model.tuning(icell).SErespModel(g);
                        end
                    end
                    map1dgain_Shf = [];
                    if isfield(obj.maps1d,'gain_Shf')
                        if ~isempty(obj.maps1d.gain_Shf{c, 2, r, o}.model)
                            map1dgain_Shf = obj.maps1d.gain_Shf{c, 2, r, o}.model.tuning(icell).meanrespModel(g);
                            map1dgainiter_Shf = obj.maps1d.gain_Shf{c, 2, r, o}.model.tuning(icell).respModel(:,g);
                        end
                    end
                    
                    if ~isempty(map1dgain)
                        obj.CellInfo.rategain{c, g, r, o}(icell) = map1dgain;
                        obj.CellInfo.rategainSE{c, g, r, o}(icell) = map1dgainSE;
                    else
                        obj.CellInfo.rategain{c, g, r, o}(icell) = 0;
                        obj.CellInfo.rategainSE{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(map1dgain_Shf)
                        obj.CellInfo.rateShfgain{c, g, r, o}(icell) = map1dgain_Shf;
                        obj.CellInfo.rateShfgainiter{c, g, r, o}(icell,:) = map1dgainiter_Shf;
                    else
                        obj.CellInfo.rateShfgain{c, g, r, o}(icell) = 0;
                        obj.CellInfo.rateShfgainiter{c, g, r, o}(icell,:) = 0;
                    end
                    
                    map2dgain_Xshift = [];
                    if isfield(obj.maps2d,'trajPercent_gain')
                        if ~isempty(obj.maps2d.trajPercent_gain{c, 2, r, o}.model)
                            map2dgain_Xshift = obj.maps2d.trajPercent_gain{c, 2, r, o}.model.tuning(icell).meancorrModelXmax(g);
                            map2dgainSE_Xshift = obj.maps2d.trajPercent_gain{c, 2, r, o}.model.tuning(icell).SEcorrModelXmax(g);
                        end
                    end
                    map2dgain_Xshift_Shf = [];
                    if isfield(obj.maps2d,'trajPercent_gain_Shf')
                        if ~isempty(obj.maps2d.trajPercent_gain_Shf{c, 2, r, o}.model)
                            map2dgain_Xshift_Shf = obj.maps2d.trajPercent_gain_Shf{c, 2, r, o}.model.tuning(icell).meancorrModelXmax(g);
                            map2dgainiter_Xshift_Shf = obj.maps2d.trajPercent_gain_Shf{c, 2, r, o}.model.tuning(icell).corrModelXmax(:,g);
                        end
                    end
                    
                    if ~isempty(map2dgain_Xshift)
                        obj.CellInfo.fieldShift{c, g, r, o}(icell) = map2dgain_Xshift;
                        obj.CellInfo.fieldShiftSE{c, g, r, o}(icell) = map2dgainSE_Xshift;
                    else
                        obj.CellInfo.fieldShift{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldShiftSE{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(map2dgain_Xshift_Shf)
                        obj.CellInfo.fieldShfShift{c, g, r, o}(icell) = map2dgain_Xshift_Shf;
                        obj.CellInfo.fieldShfShiftiter{c, g, r, o}(icell,:) = map2dgainiter_Xshift_Shf;
                    else
                        obj.CellInfo.fieldShfShift{c, g, r, o}(icell) = 0;
                        obj.CellInfo.fieldShfShiftiter{c, g, r, o}(icell,:) = 0;
                    end
                    
                    
                    unwrapmap = [];
                    if isfield(obj.maps1d,'trajPercentunwrapped')
                        if ~isempty(obj.maps1d.trajPercentunwrapped{c, g, r, o}.model)
                            unwrapmap = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).meanrespModel;
                            unwrapmapbase = obj.maps1d.trajPercentunwrapped{cbase, gbase, rbase, obase}.model.tuning(icell).meanrespModel;
                            unwrapmapSE = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).SErespModel;
                            unwrapfieldmax = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            unwrapfieldmaxSE = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            unwrapfieldAmp = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            unwrapfieldAmpSE = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            unwrapfieldAmpiter = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            unwrapfieldsize = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).fieldsize;
                            unwrapfieldsizeSE = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).SEfieldsize;
                            unwrapfieldhaldiff = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).meanrespModelXhalfdiff;
                            unwrapfieldhaldiffSE = obj.maps1d.trajPercentunwrapped{c, g, r, o}.model.tuning(icell).SErespModelXhalfdiff;
                            x = obj.maps1d.trajPercentunwrapped{c, g, r, o}.bins;
                        end
                    end
                    unwrapmap_Shf = [];
                    if isfield(obj.maps1d,'trajPercentunwrapped_Shf')
                        if ~isempty(obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model)
                            unwrapmap_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).meanrespModel;
                            unwrapmapSE_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).SErespModel;
                            unwrapfieldmax_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            unwrapfieldmaxSE_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            unwrapfieldAmp_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            unwrapfieldAmpSE_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            unwrapfieldAmpiter_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            
                            %remove unnecessary thgs from
                            %trajPercentunwrapped but most importantly save
                            % meanrespModelXhalfdiff and co.
                            
                            unwrapfieldsize_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).fieldsize;
                            unwrapfieldsizeSE_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).SEfieldsize;
                            unwrapfieldhaldiff_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXhalfdiff;
                            unwrapfieldhaldiffiter_Shf = obj.maps1d.trajPercentunwrapped_Shf{c, g, r, o}.model.tuning(icell).respModelXhalfdiff;
                        end
                    end
                    unwrapmap_2fold = [];
                    if isfield(obj.maps1d,'trajPercentunwrapped_2fold')
                        if ~isempty(obj.maps1d.trajPercentunwrapped_2fold{c, g, r, o}.model)
                            unwrapmap_2fold = obj.maps1d.trajPercentunwrapped_2fold{c, g, r, o}.model.tuning(icell).respModel;
                            unwrapfieldmax_2fold = obj.maps1d.trajPercentunwrapped_2fold{c, g, r, o}.model.tuning(icell).respModelXmax;
                        end
                    end
                    
                    if ~isempty(unwrapmap)
                        obj.CellInfo.unwrapfieldPos{c, g, r, o}(icell) = unwrapfieldmax;
                        obj.CellInfo.unwrapfieldPosSE{c, g, r, o}(icell) = unwrapfieldmaxSE;
                        obj.CellInfo.unwrapfield{c, g, r, o}(icell,:) = unwrapmap;
                        obj.CellInfo.unwrapfieldSE{c, g, r, o}(icell,:) = unwrapmapSE;
                        obj.CellInfo.unwrapfieldAmp{c, g, r, o}(icell) = unwrapfieldAmp;
                        obj.CellInfo.unwrapfieldAmpSE{c, g, r, o}(icell) = unwrapfieldAmpSE;
                        obj.CellInfo.unwrapfieldAmpiter{c, g, r, o}(icell,:) = unwrapfieldAmpiter;
                        
                        obj.CellInfo.unwrapfieldsize{c, g, r, o}(icell) = unwrapfieldsize;
                        obj.CellInfo.unwrapfieldsizeSE{c, g, r, o}(icell) = unwrapfieldsizeSE;
                        
                        obj.CellInfo.unwrapfieldhaldiff{c, g, r, o}(icell) = unwrapfieldhaldiff;
                        obj.CellInfo.unwrapfieldhaldiffSE{c, g, r, o}(icell) = unwrapfieldhaldiffSE;
                        
                        unwrapmapcorr = unwrapmap - nanmean(unwrapmap);
                        unwrapmapcorr = unwrapmapcorr./sqrt(sum(unwrapmapcorr.^2));
                        unwrapmapbasecorr = unwrapmapbase - nanmean(unwrapmapbase);
                        unwrapmapbasecorr = unwrapmapbasecorr./sqrt(sum(unwrapmapbasecorr.^2));
                        unwrapfieldXcorr = zeros(1,numel(unwrapmapbase));
                        xshiftlim = floor(numel(unwrapmapbase)/2);
                        ishift = 0;
                        for xshift = -xshiftlim:xshiftlim-1
                            ishift = ishift + 1;
                            unwrapfieldXcorr(ishift) = sum(unwrapmapcorr.*circshift(unwrapmapbasecorr,xshift))/sqrt(sum(unwrapmapcorr.^2)*sum(unwrapmapbasecorr.^2));
                        end
                        obj.CellInfo.unwrapfieldXcorr{c, g, r, o}(icell,:) = unwrapfieldXcorr;
                        [~,unwrapimax] = max(unwrapmap);
                        obj.CellInfo.unwrapfieldZ{c, g, r, o}(icell) = (unwrapmap(unwrapimax) - mean(unwrapmap))./unwrapmapSE(unwrapimax);
                        
                        obj.CellInfo.unwrapfieldMax{c, g, r, o}(icell) = max(unwrapmap);
                        obj.CellInfo.unwrapfieldMin{c, g, r, o}(icell) = min(unwrapmap);
                    else
                        obj.CellInfo.unwrapfieldPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfieldPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfield{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.unwrapfieldSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.unwrapfieldsize{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldsizeSE{c, g, r, o}(icell) = NaN;
                        
                        obj.CellInfo.unwrapfieldAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfieldAmpiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.unwrapfieldhaldiff{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldhaldiffSE{c, g, r, o}(icell) = NaN;
                        
                        obj.CellInfo.unwrapfieldXcorr{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.unwrapfieldZ{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldMax{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfieldMin{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(unwrapmap_Shf)
                        obj.CellInfo.unwrapfieldShfPos{c, g, r, o}(icell) = unwrapfieldmax_Shf;
                        obj.CellInfo.unwrapfieldShfPosSE{c, g, r, o}(icell) = unwrapfieldmaxSE_Shf;
                        obj.CellInfo.unwrapfieldShf{c, g, r, o}(icell,:) = unwrapmap_Shf;
                        obj.CellInfo.unwrapfieldShfSE{c, g, r, o}(icell,:) = unwrapmapSE_Shf;
                        obj.CellInfo.unwrapfieldShfsize{c, g, r, o}(icell) = unwrapfieldsize_Shf;
                        obj.CellInfo.unwrapfieldShfsizeSE{c, g, r, o}(icell) = unwrapfieldsizeSE_Shf;
                        
                        obj.CellInfo.unwrapfieldShfAmp{c, g, r, o}(icell) = unwrapfieldAmp_Shf;
                        obj.CellInfo.unwrapfieldShfAmpSE{c, g, r, o}(icell) = unwrapfieldAmpSE_Shf;
                        obj.CellInfo.unwrapfieldShfAmpiter{c, g, r, o}(icell,:) = unwrapfieldAmpiter_Shf;
                        
                        obj.CellInfo.unwrapfieldShfhaldiff{c, g, r, o}(icell) = unwrapfieldhaldiff_Shf;
                        obj.CellInfo.unwrapfieldShfhaldiffiter{c, g, r, o}(icell,:) = unwrapfieldhaldiffiter_Shf;
                    else
                        obj.CellInfo.unwrapfieldShfPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfieldShfPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfieldShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.unwrapfieldShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.unwrapfieldShfsize{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldShfsizeSE{c, g, r, o}(icell) = NaN;
                        
                        obj.CellInfo.unwrapfieldShfAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldShfAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.unwrapfieldShfAmpiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.unwrapfieldShfhaldiff{c, g, r, o}(icell) = 0;
                        obj.CellInfo.unwrapfieldShfhaldiffiter{c, g, r, o}(icell,:) = 0;
                    end
                    if ~isempty(unwrapmap_2fold)
                        obj.CellInfo.unwrapfield_2fold{c, g, r, o}(icell,:,:) = unwrapmap_2fold;
                        obj.CellInfo.unwrapfieldPos_2fold{c, g, r, o}(icell,:) = unwrapfieldmax_2fold;
                        unwrapfieldset1 = squeeze(unwrapmap_2fold(1,:));
                        unwrapfieldset2 = squeeze(unwrapmap_2fold(2,:));
                        unwrapreliabilityCorr = sum((unwrapfieldset1-repmat(nanmean(unwrapfieldset1,2),[1 size(unwrapfieldset1,2)])).*...
                            (unwrapfieldset2-repmat(nanmean(unwrapfieldset2,2),[1 size(unwrapfieldset2,2)])),2)./...
                            sqrt(sum((unwrapfieldset1-repmat(nanmean(unwrapfieldset1,2),[1 size(unwrapfieldset1,2)])).^2,2).*...
                            sum((unwrapfieldset2-repmat(nanmean(unwrapfieldset2,2),[1 size(unwrapfieldset2,2)])).^2,2));
                        obj.CellInfo.unwrapreliabilityCorr{c, g, r, o}(icell) = unwrapreliabilityCorr;
                    else
                        obj.CellInfo.unwrapfield_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.unwrapfieldPos_2fold{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.unwrapreliabilityCorr{c, g, r, o}(icell) = NaN;
                    end
                    
                    map2d = [];
                    if isfield(obj.maps2d,'trajPercent_smthBallSpd')
                        if ~isempty(obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model)
                            map2d = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModel;
                            map2dSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModel;
                            map2dX = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModelX;
                            map2dXSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModelX;
                            map2dXPos = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            map2dXPosSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            map2dXmax = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            map2dXmaxSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            
                            map2dY = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModelY;
                            map2dYSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModelY;
                            map2dYPos = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModelYpos;
                            map2dYPosSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModelYpos;
                            map2dYmax = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).meanrespModelYmax;
                            map2dYmaxSE = obj.maps2d.trajPercent_smthBallSpd{c, g, r, o}.model.tuning(icell).SErespModelYmax;
                        end
                    end
                    
                    if ~isempty(map2d)
                        obj.CellInfo.field2dXspd{c, g, r, o}(icell,:,:) = map2d;
                        obj.CellInfo.field2dXspdSE{c, g, r, o}(icell,:,:) = map2dSE;
                        obj.CellInfo.field2dX{c, g, r, o}(icell,:) = map2dX(:)';
                        obj.CellInfo.field2dXSE{c, g, r, o}(icell,:) = map2dXSE(:)';
                        obj.CellInfo.field2dXpos{c, g, r, o}(icell,:) = map2dXPos(:)';
                        obj.CellInfo.field2dXposSE{c, g, r, o}(icell,:) = map2dXPosSE(:)';
                        obj.CellInfo.field2dXmax{c, g, r, o}(icell,:) = map2dXmax(:)';
                        obj.CellInfo.field2dXmaxSE{c, g, r, o}(icell,:) = map2dXmaxSE(:)';
                        
                        obj.CellInfo.field2dSpd{c, g, r, o}(icell,:) = map2dY(:)';
                        obj.CellInfo.field2dSpdSE{c, g, r, o}(icell,:) = map2dYSE(:)';
                        obj.CellInfo.field2dSpdpos{c, g, r, o}(icell,:) = map2dYPos(:)';
                        obj.CellInfo.field2dSpdposSE{c, g, r, o}(icell,:) = map2dYPosSE(:)';
                        obj.CellInfo.field2dSpdmax{c, g, r, o}(icell,:) = map2dYmax(:)';
                        obj.CellInfo.field2dSpdmaxSE{c, g, r, o}(icell,:) = map2dYmaxSE(:)';
                    else
                        obj.CellInfo.field2dXspd{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXspdSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dX{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXmaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dSpd{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dSpdSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dSpdpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dSpdposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dSpdmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dSpdmaxSE{c, g, r, o}(icell,:) = 0;
                    end
                    
                    %
                    map2d_theta = [];
                    if isfield(obj.maps2d,'trajPercent_LFPphase2')
                        if ~isempty(obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model)
                            map2d_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModel;                           
                            
                            map2dSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModel;
                            map2dX_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelX;
                            map2dXSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelX;
                            map2dXPos_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            map2dXPosSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            map2dXmax_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            map2dXmaxSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            
                            map2dY_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelY;
                            map2dYSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelY;
                            map2dYPos_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelYpos;
                            map2dYPosSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelYpos;
                            map2dYmax_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelYmax;
                            map2dYmaxSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelYmax;
                            
                            map2dslopeXY = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelslopeXY;
                            map2dphi0XY = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelphi0XY;
                            map2drhoXY = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelrhoXY;
                            map2dslopeXYSE = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelslopeXY;
                            map2dphi0XYSE = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelphi0XY;
                            map2drhoXYSE = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelrhoXY;
                            map2dslopeXYiter = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).respModelslopeXY;
                            map2dphi0XYiter = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).respModelphi0XY;
                            map2drhoXYiter = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).respModelrhoXY;
                            
                            map2dXcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelX;
                            map2dXcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelX;
                            map2dXPoscorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXpos;
                            map2dXPoscorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXpos;
                            map2dXmaxcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXmax;
                            map2dXmaxcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXmax;
                            map2dXmaxAmpcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXmaxAmp;
                            map2dXmaxAmpcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxAmp;
                            map2dXmaxAmpcorriter_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).corrModelXmaxAmp;
                            map2dXmaxsinAmpcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXmaxsinAmp;
                            map2dXmaxsinAmpcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxsinAmp;
                            map2dXmaxsinAmpcorriter_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).corrModelXmaxsinAmp;
                            map2dXmaxOffsetcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXmaxOffset;
                            map2dXmaxOffsetcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxOffset;
                            map2dXmaxOffsetcorriter_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).corrModelXmaxOffset;
                            
                            map2dXRefcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXRef;
                            map2dXRefcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXRef;
                            map2dXRefPoscorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXRefpos;
                            map2dXRefPoscorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXRefpos;
                            map2dXRefmaxcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXRefmax;
                            map2dXRefmaxcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmax;
                            map2dXRefmaxAmpcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxAmp;
                            map2dXRefmaxAmpcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxAmp;
                            map2dXRefmaxAmpcorriter_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).corrModelXRefmaxAmp;
                            map2dXRefmaxsinAmpcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxsinAmp;
                            map2dXRefmaxsinAmpcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxsinAmp;
                            map2dXRefmaxsinAmpcorriter_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).corrModelXRefmaxsinAmp;
                            map2dXRefmaxOffsetcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxOffset;
                            map2dXRefmaxOffsetcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxOffset;
                            map2dXRefmaxOffsetcorriter_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).corrModelXRefmaxOffset;
                            
                            map2dYcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelY;
                            map2dYcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelY;
                            map2dYPoscorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelYpos;
                            map2dYPoscorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelYpos;
                            map2dYmaxcorr_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).meancorrModelYmax;
                            map2dYmaxcorrSE_theta = obj.maps2d.trajPercent_LFPphase2{c, g, r, o}.model.tuning(icell).SEcorrModelYmax;
                            
                            xmap = map2dX_theta;
                            [~,imax] = max(xmap);
                            Xmax = max(obj.data.es.trajPercent);
                            varXtemp = mod(obj.data.es.trajPercent - imax + Xmax/2,Xmax)-Xmax/2;
                            xmap = circshift(xmap,-imax+floor(numel(xmap)/2));
                            fieldstart = x(find(xmap(1:floor(numel(xmap)/2))<=mean(xmap),1,'last'))-floor(numel(xmap)/2);
                            fieldend = x(find(xmap(floor(numel(xmap)/2)+1:end)<=mean(xmap),1,'first'));
                            spktraintemp = obj.data.es.spikeTrain(:,icell);
                            spktraintemp(varXtemp<=fieldstart | varXtemp>=fieldend) = 0;
                            [rho,pval] = circ_corrcc(obj.data.es.spikePhase(tidx & spktraintemp>0.5,icell)/360*2*pi, varXtemp(tidx & spktraintemp>0.5)/Xmax*2*pi);
                            
                            phsmap = map2dY_theta;
                            obj.CellInfo.Phase{c, g, r, o}(icell,:) = phsmap(:)';
                            y = linspace(0,360,numel(phsmap));
                            phsmap = phsmap/sum(phsmap)*sum(obj.data.es.spikeTrain(tidx,icell),1);
                            [rtest_pval,rtest_z] = circ_rtest(y/360*2*pi,phsmap);
                            
                            [~, phsidx] = max(map2dY_theta);
                            phsmax = y(phsidx);
                        end
                    end
                    map2d_theta_Shf = [];
                    if isfield(obj.maps2d,'trajPercent_LFPphase2_Shf')
                        if ~isempty(obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model)
                            map2d_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModel;                           
                            
                            map2dSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModel;
                            map2dX_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelX;
                            map2dXSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelX;
                            map2dXPos_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            map2dXPosSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            map2dXmax_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            map2dXmaxSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            
                            map2dY_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelY;
                            map2dYSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelY;
                            map2dYPos_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelYpos;
                            map2dYPosSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelYpos;
                            map2dYmax_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelYmax;
                            map2dYmaxSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelYmax;
                            
                            map2dslopeXY_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelslopeXY;
                            map2dphi0XY_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelphi0XY;
                            map2drhoXY_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelrhoXY;
                            map2dslopeXYSE_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelslopeXY;
                            map2dphi0XYSE_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelphi0XY;
                            map2drhoXYSE_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelrhoXY;
                            map2dslopeXYiter_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).respModelslopeXY;
                            map2dphi0XYiter_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).respModelphi0XY;
                            map2drhoXYiter_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).respModelrhoXY;
                            
                            map2dXcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelX;
                            map2dXcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelX;
                            map2dXPoscorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXpos;
                            map2dXPoscorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXpos;
                            map2dXmaxcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmax;
                            map2dXmaxcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmax;
                            map2dXmaxAmpcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmaxAmp;
                            map2dXmaxAmpcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxAmp;
                            map2dXmaxAmpcorriter_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).corrModelXmaxAmp;
                            map2dXmaxsinAmpcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmaxsinAmp;
                            map2dXmaxsinAmpcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxsinAmp;
                            map2dXmaxsinAmpcorriter_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).corrModelXmaxsinAmp;
                            map2dXmaxOffsetcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmaxOffset;
                            map2dXmaxOffsetcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxOffset;
                            map2dXmaxOffsetcorriter_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).corrModelXmaxOffset;
                            
                            map2dXRefcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRef;
                            map2dXRefcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRef;
                            map2dXRefPoscorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefpos;
                            map2dXRefPoscorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefpos;
                            map2dXRefmaxcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmax;
                            map2dXRefmaxcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmax;
                            map2dXRefmaxAmpcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxAmp;
                            map2dXRefmaxAmpcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxAmp;
                            map2dXRefmaxAmpcorriter_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).corrModelXRefmaxAmp;
                            map2dXRefmaxsinAmpcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxsinAmp;
                            map2dXRefmaxsinAmpcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxsinAmp;
                            map2dXRefmaxsinAmpcorriter_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).corrModelXRefmaxsinAmp;
                            map2dXRefmaxOffsetcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxOffset;
                            map2dXRefmaxOffsetcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxOffset;
                            map2dXRefmaxOffsetcorriter_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).corrModelXRefmaxOffset;
                            
                            map2dYcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelY;
                            map2dYcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelY;
                            map2dYPoscorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelYpos;
                            map2dYPoscorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelYpos;
                            map2dYmaxcorr_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meancorrModelYmax;
                            map2dYmaxcorrSE_theta_Shf = obj.maps2d.trajPercent_LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelYmax;
                        end
                    end
                    map2d_theta_2fold = [];
                    if isfield(obj.maps2d,'trajPercent_LFPphase2_2fold')
                        if ~isempty(obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model)
                            map2dXcorr_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).corrModelX;
                            map2dXPoscorr_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).corrModelXpos;
                            map2dXmaxcorr_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).corrModelXmax;
                            
                            map2dYcorr_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).corrModelY;
                            map2dYPoscorr_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).corrModelYpos;
                            map2dYmaxcorr_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).corrModelYmax;
                            
                            map2d_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModel;
                            map2dX_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelX;
                            map2dXPos_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelXpos;
                            map2dXmax_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelXmax;
                            
                            map2dY_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelY;
                            map2dYPos_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelYpos;
                            map2dYmax_theta_2fold = obj.maps2d.trajPercent_LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelYmax;
                        end
                    end
                    if ~isempty(map2d_theta)
                        obj.CellInfo.field2dXPhstheta{c, g, r, o}(icell,:,:) = map2d_theta;
                        obj.CellInfo.field2dXPhsthetaSE{c, g, r, o}(icell,:,:) = map2dSE_theta;
                        obj.CellInfo.field2dXtheta{c, g, r, o}(icell,:) = map2dX_theta(:)';
                        obj.CellInfo.field2dXthetaSE{c, g, r, o}(icell,:) = map2dXSE_theta(:)';
                        obj.CellInfo.field2dXthetapos{c, g, r, o}(icell,:) = map2dXPos_theta(:)';
                        obj.CellInfo.field2dXthetaposSE{c, g, r, o}(icell,:) = map2dXPosSE_theta(:)';
                        obj.CellInfo.field2dXthetamax{c, g, r, o}(icell,:) = map2dXmax_theta(:)';
                        obj.CellInfo.field2dXthetamaxSE{c, g, r, o}(icell,:) = map2dXmaxSE_theta(:)';
                        
                        obj.CellInfo.field2dPhstheta{c, g, r, o}(icell,:) = map2dY_theta(:)';
                        obj.CellInfo.field2dPhsthetaSE{c, g, r, o}(icell,:) = map2dYSE_theta(:)';
                        obj.CellInfo.field2dPhsthetapos{c, g, r, o}(icell,:) = map2dYPos_theta(:)';
                        obj.CellInfo.field2dPhsthetaposSE{c, g, r, o}(icell,:) = map2dYPosSE_theta(:)';
                        obj.CellInfo.field2dPhsthetamax{c, g, r, o}(icell,:) = map2dYmax_theta(:)';
                        obj.CellInfo.field2dPhsthetamaxSE{c, g, r, o}(icell,:) = map2dYmaxSE_theta(:)';
                        
                        obj.CellInfo.field2dslopeXY{c, g, r, o}(icell) = map2dslopeXY;
                        obj.CellInfo.field2dphi0XY{c, g, r, o}(icell) = map2dphi0XY;
                        obj.CellInfo.field2drhoXY{c, g, r, o}(icell) = map2drhoXY;
                        obj.CellInfo.field2dslopeXYSE{c, g, r, o}(icell) = map2dslopeXYSE;
                        obj.CellInfo.field2dphi0XYSE{c, g, r, o}(icell) = map2dphi0XYSE;
                        obj.CellInfo.field2drhoXYSE{c, g, r, o}(icell) = map2drhoXYSE;
                        obj.CellInfo.field2dslopeXYiter{c, g, r, o}(icell,:) = map2dslopeXYiter;
                        obj.CellInfo.field2dphi0XYiter{c, g, r, o}(icell,:) = map2dphi0XYiter;
                        obj.CellInfo.field2drhoXYiter{c, g, r, o}(icell,:) = map2drhoXYiter;
                        
                        obj.CellInfo.field2dXcorrtheta{c, g, r, o}(icell,:,:) = map2dXcorr_theta;
                        obj.CellInfo.field2dXcorrthetaSE{c, g, r, o}(icell,:,:) = map2dXcorrSE_theta;
                        obj.CellInfo.field2dXcorrthetapos{c, g, r, o}(icell,:) = map2dXPoscorr_theta;
                        obj.CellInfo.field2dXcorrthetaposSE{c, g, r, o}(icell,:) = map2dXPoscorrSE_theta;
                        obj.CellInfo.field2dXcorrthetamax{c, g, r, o}(icell,:) = map2dXmaxcorr_theta;
                        obj.CellInfo.field2dXcorrthetamaxSE{c, g, r, o}(icell,:) = map2dXmaxcorrSE_theta;
                        obj.CellInfo.field2dXcorrthetamaxAmp{c, g, r, o}(icell) = map2dXmaxAmpcorr_theta;
                        obj.CellInfo.field2dXcorrthetamaxAmpSE{c, g, r, o}(icell) = map2dXmaxAmpcorrSE_theta;
                        obj.CellInfo.field2dXcorrthetamaxAmpiter{c, g, r, o}(icell,:) = map2dXmaxAmpcorriter_theta;
                        obj.CellInfo.field2dXcorrthetamaxsinAmp{c, g, r, o}(icell) = map2dXmaxsinAmpcorr_theta;
                        obj.CellInfo.field2dXcorrthetamaxsinAmpSE{c, g, r, o}(icell) = map2dXmaxsinAmpcorrSE_theta;
                        obj.CellInfo.field2dXcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = map2dXmaxsinAmpcorriter_theta;
                        obj.CellInfo.field2dXcorrthetamaxOffset{c, g, r, o}(icell) = map2dXmaxOffsetcorr_theta;
                        obj.CellInfo.field2dXcorrthetamaxOffsetSE{c, g, r, o}(icell) = map2dXmaxOffsetcorrSE_theta;
                        obj.CellInfo.field2dXcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = map2dXmaxOffsetcorriter_theta;
                        
                        obj.CellInfo.field2dXRefcorrtheta{c, g, r, o}(icell,:,:) = map2dXRefcorr_theta;
                        obj.CellInfo.field2dXRefcorrthetaSE{c, g, r, o}(icell,:,:) = map2dXRefcorrSE_theta;
                        obj.CellInfo.field2dXRefcorrthetapos{c, g, r, o}(icell,:) = map2dXRefPoscorr_theta;
                        obj.CellInfo.field2dXRefcorrthetaposSE{c, g, r, o}(icell,:) = map2dXRefPoscorrSE_theta;
                        obj.CellInfo.field2dXRefcorrthetamax{c, g, r, o}(icell,:) = map2dXRefmaxcorr_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxSE{c, g, r, o}(icell,:) = map2dXRefmaxcorrSE_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxAmp{c, g, r, o}(icell) = map2dXRefmaxAmpcorr_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxAmpSE{c, g, r, o}(icell) = map2dXRefmaxAmpcorrSE_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxAmpiter{c, g, r, o}(icell,:) = map2dXRefmaxAmpcorriter_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxsinAmp{c, g, r, o}(icell) = map2dXRefmaxsinAmpcorr_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxsinAmpSE{c, g, r, o}(icell) = map2dXRefmaxsinAmpcorrSE_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = map2dXRefmaxsinAmpcorriter_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxOffset{c, g, r, o}(icell) = map2dXRefmaxOffsetcorr_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxOffsetSE{c, g, r, o}(icell) = map2dXRefmaxOffsetcorrSE_theta;
                        obj.CellInfo.field2dXRefcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = map2dXRefmaxOffsetcorriter_theta;
                        
                        obj.CellInfo.field2dPhscorrtheta{c, g, r, o}(icell,:,:) = map2dYcorr_theta;
                        obj.CellInfo.field2dPhscorrthetaSE{c, g, r, o}(icell,:,:) = map2dYcorrSE_theta;
                        obj.CellInfo.field2dPhscorrthetapos{c, g, r, o}(icell,:) = map2dYPoscorr_theta;
                        obj.CellInfo.field2dPhscorrthetaposSE{c, g, r, o}(icell,:) = map2dYPoscorrSE_theta;
                        obj.CellInfo.field2dPhscorrthetamax{c, g, r, o}(icell,:) = map2dYmaxcorr_theta;
                        obj.CellInfo.field2dPhscorrthetamaxSE{c, g, r, o}(icell,:) = map2dYmaxcorrSE_theta;
                        
                        obj.CellInfo.PhaseRho{c, g, r, o}(icell) = rho;
                        obj.CellInfo.PhasePval{c, g, r, o}(icell) = pval;
                        obj.CellInfo.PhaseRayleighPval{c, g, r, o}(icell) = rtest_pval;
                        obj.CellInfo.PhaseRayleighZ{c, g, r, o}(icell) = rtest_z;
                        obj.CellInfo.PhaseMax{c, g, r, o}(icell) = phsmax;
                    else
                        obj.CellInfo.field2dXPhstheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXPhsthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXtheta{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetamaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dPhstheta{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetamaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dslopeXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dphi0XY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2drhoXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dslopeXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dphi0XYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2drhoXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dslopeXYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dphi0XYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2drhoXYiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dXcorrtheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXcorrthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXcorrthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetamaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetamaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXcorrthetamaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXcorrthetamaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetamaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXcorrthetamaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetamaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXcorrthetamaxOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dXRefcorrtheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXRefcorrthetamaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXRefcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxOffsetSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dPhscorrtheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dPhscorrthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhscorrthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhscorrthetamaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.PhaseRho{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.PhasePval{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.PhaseRayleighPval{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.PhaseRayleighZ{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.PhaseMax{c, g, r, o}(icell) = NaN;
                    end
                    
                    if ~isempty(map2d_theta_Shf)
                        obj.CellInfo.field2dXPhsthetaShf{c, g, r, o}(icell,:,:) = map2d_theta_Shf;
                        obj.CellInfo.field2dXPhsthetaShfSE{c, g, r, o}(icell,:,:) = map2dSE_theta_Shf;
                        obj.CellInfo.field2dXthetaShf{c, g, r, o}(icell,:) = map2dX_theta_Shf(:)';
                        obj.CellInfo.field2dXthetaShfSE{c, g, r, o}(icell,:) = map2dXSE_theta_Shf(:)';
                        obj.CellInfo.field2dXthetaShfpos{c, g, r, o}(icell,:) = map2dXPos_theta_Shf(:)';
                        obj.CellInfo.field2dXthetaShfposSE{c, g, r, o}(icell,:) = map2dXPosSE_theta_Shf(:)';
                        obj.CellInfo.field2dXthetaShfmax{c, g, r, o}(icell,:) = map2dXmax_theta_Shf(:)';
                        obj.CellInfo.field2dXthetaShfmaxSE{c, g, r, o}(icell,:) = map2dXmaxSE_theta_Shf(:)';
                        
                        obj.CellInfo.field2dPhsthetaShf{c, g, r, o}(icell,:) = map2dY_theta_Shf(:)';
                        obj.CellInfo.field2dPhsthetaShfSE{c, g, r, o}(icell,:) = map2dYSE_theta_Shf(:)';
                        obj.CellInfo.field2dPhsthetaShfpos{c, g, r, o}(icell,:) = map2dYPos_theta_Shf(:)';
                        obj.CellInfo.field2dPhsthetaShfposSE{c, g, r, o}(icell,:) = map2dYPosSE_theta_Shf(:)';
                        obj.CellInfo.field2dPhsthetaShfmax{c, g, r, o}(icell,:) = map2dYmax_theta_Shf(:)';
                        obj.CellInfo.field2dPhsthetaShfmaxSE{c, g, r, o}(icell,:) = map2dYmaxSE_theta_Shf(:)';
                        
                        obj.CellInfo.field2dShfslopeXY{c, g, r, o}(icell) = map2dslopeXY_Shf;
                        obj.CellInfo.field2dShfphi0XY{c, g, r, o}(icell) = map2dphi0XY_Shf;
                        obj.CellInfo.field2dShfrhoXY{c, g, r, o}(icell) = map2drhoXY_Shf;
                        obj.CellInfo.field2dShfslopeXYSE{c, g, r, o}(icell) = map2dslopeXYSE_Shf;
                        obj.CellInfo.field2dShfphi0XYSE{c, g, r, o}(icell) = map2dphi0XYSE_Shf;
                        obj.CellInfo.field2dShfrhoXYSE{c, g, r, o}(icell) = map2drhoXYSE_Shf;
                        obj.CellInfo.field2dShfslopeXYiter{c, g, r, o}(icell,:) = map2dslopeXYiter_Shf;
                        obj.CellInfo.field2dShfphi0XYiter{c, g, r, o}(icell,:) = map2dphi0XYiter_Shf;
                        obj.CellInfo.field2dShfrhoXYiter{c, g, r, o}(icell,:) = map2drhoXYiter_Shf;
                        
                        obj.CellInfo.field2dXcorrthetaShf{c, g, r, o}(icell,:,:) = map2dXcorr_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfSE{c, g, r, o}(icell,:,:) = map2dXcorrSE_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfpos{c, g, r, o}(icell,:) = map2dXPoscorr_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfposSE{c, g, r, o}(icell,:) = map2dXPoscorrSE_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmax{c, g, r, o}(icell,:) = map2dXmaxcorr_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxSE{c, g, r, o}(icell,:) = map2dXmaxcorrSE_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxAmp{c, g, r, o}(icell) = map2dXmaxAmpcorr_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = map2dXmaxAmpcorrSE_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = map2dXmaxAmpcorriter_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = map2dXmaxsinAmpcorr_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = map2dXmaxsinAmpcorrSE_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = map2dXmaxsinAmpcorriter_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxOffset{c, g, r, o}(icell) = map2dXmaxOffsetcorr_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = map2dXmaxOffsetcorrSE_theta_Shf;
                        obj.CellInfo.field2dXcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = map2dXmaxOffsetcorriter_theta_Shf;
                        
                        obj.CellInfo.field2dXRefcorrthetaShf{c, g, r, o}(icell,:,:) = map2dXRefcorr_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfSE{c, g, r, o}(icell,:,:) = map2dXRefcorrSE_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfpos{c, g, r, o}(icell,:) = map2dXRefPoscorr_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfposSE{c, g, r, o}(icell,:) = map2dXRefPoscorrSE_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmax{c, g, r, o}(icell,:) = map2dXRefmaxcorr_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxSE{c, g, r, o}(icell,:) = map2dXRefmaxcorrSE_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxAmp{c, g, r, o}(icell) = map2dXRefmaxAmpcorr_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = map2dXRefmaxAmpcorrSE_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = map2dXRefmaxAmpcorriter_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = map2dXRefmaxsinAmpcorr_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = map2dXRefmaxsinAmpcorrSE_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = map2dXRefmaxsinAmpcorriter_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxOffset{c, g, r, o}(icell) = map2dXRefmaxOffsetcorr_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = map2dXRefmaxOffsetcorrSE_theta_Shf;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = map2dXRefmaxOffsetcorriter_theta_Shf;
                        
                        obj.CellInfo.field2dPhscorrthetaShf{c, g, r, o}(icell,:,:) = map2dYcorr_theta_Shf;
                        obj.CellInfo.field2dPhscorrthetaShfSE{c, g, r, o}(icell,:,:) = map2dYcorrSE_theta_Shf;
                        obj.CellInfo.field2dPhscorrthetaShfpos{c, g, r, o}(icell,:) = map2dYPoscorr_theta_Shf;
                        obj.CellInfo.field2dPhscorrthetaShfposSE{c, g, r, o}(icell,:) = map2dYPoscorrSE_theta_Shf;
                        obj.CellInfo.field2dPhscorrthetaShfmax{c, g, r, o}(icell,:) = map2dYmaxcorr_theta_Shf;
                        obj.CellInfo.field2dPhscorrthetaShfmaxSE{c, g, r, o}(icell,:) = map2dYmaxcorrSE_theta_Shf;                     
                    else
                        obj.CellInfo.field2dXPhsthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXPhsthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXthetaShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dPhsthetaShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhsthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dShfslopeXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dShfphi0XY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dShfrhoXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dShfslopeXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dShfphi0XYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dShfrhoXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dShfslopeXYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dShfphi0XYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dShfrhoXYiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dXcorrthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dXRefcorrthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.field2dXRefcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.field2dPhscorrthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.field2dPhscorrthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                    end
                    
                    if ~isempty(map2d_theta_2fold)
                        obj.CellInfo.field2dXcorrtheta_2fold{c, g, r, o}(icell,:,:,:) = map2dXcorr_theta_2fold;
                        obj.CellInfo.field2dXcorrthetapos_2fold{c, g, r, o}(icell,:,:) = map2dXPoscorr_theta_2fold;
                        obj.CellInfo.field2dXcorrthetamax_2fold{c, g, r, o}(icell,:,:) = map2dXmaxcorr_theta_2fold;
                        
                        obj.CellInfo.field2dPhscorrtheta_2fold{c, g, r, o}(icell,:,:,:) = map2dYcorr_theta_2fold;
                        obj.CellInfo.field2dPhscorrthetapos_2fold{c, g, r, o}(icell,:,:) = map2dYPoscorr_theta_2fold;
                        obj.CellInfo.field2dPhscorrthetamax_2fold{c, g, r, o}(icell,:,:) = map2dYmaxcorr_theta_2fold;
                        
                        obj.CellInfo.field2dXPhstheta_2fold{c, g, r, o}(icell,:,:,:) = map2d_theta_2fold;
                        obj.CellInfo.field2dXtheta_2fold{c, g, r, o}(icell,:,:) = map2dX_theta_2fold;
                        obj.CellInfo.field2dXthetapos_2fold{c, g, r, o}(icell,:,:) = map2dXPos_theta_2fold;
                        obj.CellInfo.field2dXthetamax_2fold{c, g, r, o}(icell,:,:) = map2dXmax_theta_2fold;
                        
                        obj.CellInfo.field2dYtheta_2fold{c, g, r, o}(icell,:,:) = map2dY_theta_2fold;
                        obj.CellInfo.field2dYthetapos_2fold{c, g, r, o}(icell,:,:) = map2dYPos_theta_2fold;
                        obj.CellInfo.field2dYthetamax_2fold{c, g, r, o}(icell,:,:) = map2dYmax_theta_2fold;
                    else
                        obj.CellInfo.field2dXcorrtheta_2fold{c, g, r, o}(icell,:,:,:) = 0;
                        obj.CellInfo.field2dXcorrthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXcorrthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                        
                        obj.CellInfo.field2dPhscorrtheta_2fold{c, g, r, o}(icell,:,:,:) = 0;
                        obj.CellInfo.field2dPhscorrthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dPhscorrthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                        
                        obj.CellInfo.field2dXPhstheta_2fold{c, g, r, o}(icell,:,:,:) = 0;
                        obj.CellInfo.field2dXtheta_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dXthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                        
                        obj.CellInfo.field2dYtheta_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dYthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.field2dYthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                    end
                    %
                    
                    spkmap2d_theta = [];
                    if isfield(obj.maps2d,'spiketrajPercent_spikePhase')
                        if ~isempty(obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model)
                            spkmap2d_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModel;                           
                            
                            spkmap2dSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModel;
                            spkmap2dX_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelX;
                            spkmap2dXSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelX;
                            spkmap2dXPos_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            spkmap2dXPosSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            spkmap2dXmax_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            spkmap2dXmaxSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            
                            spkmap2dY_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelY;
                            spkmap2dYSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelY;
                            spkmap2dYPos_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelYpos;
                            spkmap2dYPosSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelYpos;
                            spkmap2dYmax_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelYmax;
                            spkmap2dYmaxSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelYmax;
                            
                            spkmap2dslopeXY = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelslopeXY;
                            spkmap2dphi0XY = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelphi0XY;
                            spkmap2drhoXY = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelrhoXY;
                            spkmap2dslopeXYSE = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelslopeXY;
                            spkmap2dphi0XYSE = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelphi0XY;
                            spkmap2drhoXYSE = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SErespModelrhoXY;
                            spkmap2dslopeXYiter = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).respModelslopeXY;
                            spkmap2dphi0XYiter = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).respModelphi0XY;
                            spkmap2drhoXYiter = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).respModelrhoXY;
                            
                            spkmap2dXcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelX;
                            spkmap2dXcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelX;
                            spkmap2dXPoscorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXpos;
                            spkmap2dXPoscorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXpos;
                            spkmap2dXmaxcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXmax;
                            spkmap2dXmaxcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXmax;
                            spkmap2dXmaxAmpcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXmaxAmp;
                            spkmap2dXmaxAmpcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxAmp;
                            spkmap2dXmaxAmpcorriter_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).corrModelXmaxAmp;
                            spkmap2dXmaxsinAmpcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXmaxsinAmp;
                            spkmap2dXmaxsinAmpcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxsinAmp;
                            spkmap2dXmaxsinAmpcorriter_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).corrModelXmaxsinAmp;
                            spkmap2dXmaxOffsetcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXmaxOffset;
                            spkmap2dXmaxOffsetcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxOffset;
                            spkmap2dXmaxOffsetcorriter_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).corrModelXmaxOffset;
                            
                            spkmap2dXRefcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXRef;
                            spkmap2dXRefcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXRef;
                            spkmap2dXRefPoscorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXRefpos;
                            spkmap2dXRefPoscorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXRefpos;
                            spkmap2dXRefmaxcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXRefmax;
                            spkmap2dXRefmaxcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmax;
                            spkmap2dXRefmaxAmpcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxAmp;
                            spkmap2dXRefmaxAmpcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxAmp;
                            spkmap2dXRefmaxAmpcorriter_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).corrModelXRefmaxAmp;
                            spkmap2dXRefmaxsinAmpcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxsinAmp;
                            spkmap2dXRefmaxsinAmpcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxsinAmp;
                            spkmap2dXRefmaxsinAmpcorriter_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).corrModelXRefmaxsinAmp;
                            spkmap2dXRefmaxOffsetcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxOffset;
                            spkmap2dXRefmaxOffsetcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxOffset;
                            spkmap2dXRefmaxOffsetcorriter_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).corrModelXRefmaxOffset;
                            
                            spkmap2dYcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelY;
                            spkmap2dYcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelY;
                            spkmap2dYPoscorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelYpos;
                            spkmap2dYPoscorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelYpos;
                            spkmap2dYmaxcorr_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).meancorrModelYmax;
                            spkmap2dYmaxcorrSE_theta = obj.maps2d.spiketrajPercent_spikePhase{c, g, r, o}.model.tuning(icell).SEcorrModelYmax;
                            
                            xmap = spkmap2dX_theta;
                            [~,imax] = max(xmap);
                            Xmax = max(obj.data.es.spiketrajPercent(:,icell));
                            varXtemp = mod(obj.data.es.spiketrajPercent(:,icell) - imax + Xmax/2,Xmax)-Xmax/2;
                            xmap = circshift(xmap,-imax+floor(numel(xmap)/2));
                            fieldstart = x(find(xmap(1:floor(numel(xmap)/2))<=mean(xmap),1,'last'))-floor(numel(xmap)/2);
                            fieldend = x(find(xmap(floor(numel(xmap)/2)+1:end)<=mean(xmap),1,'first'));
                            spktraintemp = obj.data.es.spikeTrain(:,icell);
                            spktraintemp(varXtemp<=fieldstart | varXtemp>=fieldend) = 0;
                            [spkrho,spkpval] = circ_corrcc(obj.data.es.spikePhase(tidx & spktraintemp>0.5,icell)/360*2*pi, varXtemp(tidx & spktraintemp>0.5)/Xmax*2*pi);
                            
                            spkphsmap = spkmap2dY_theta;
                            obj.CellInfo.spkPhase{c, g, r, o}(icell,:) = spkphsmap(:)';
                            y = linspace(0,360,numel(spkphsmap));
                            spkphsmap = spkphsmap/sum(spkphsmap)*sum(obj.data.es.spikeTrain(tidx,icell),1);
                            [spkrtest_pval,spkrtest_z] = circ_rtest(y/360*2*pi,spkphsmap);
                            
                            [~, spkphsidx] = max(spkmap2dY_theta);
                            spkphsmax = y(spkphsidx);
                        end
                    end
                    spkmap2d_theta_Shf = [];
                    if isfield(obj.maps2d,'spiketrajPercent_spikePhase_Shf')
                        if ~isempty(obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model)
                            spkmap2d_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModel;                           
                            
                            spkmap2dSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModel;
                            spkmap2dX_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelX;
                            spkmap2dXSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelX;
                            spkmap2dXPos_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            spkmap2dXPosSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            spkmap2dXmax_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            spkmap2dXmaxSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            
                            spkmap2dY_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelY;
                            spkmap2dYSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelY;
                            spkmap2dYPos_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelYpos;
                            spkmap2dYPosSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelYpos;
                            spkmap2dYmax_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelYmax;
                            spkmap2dYmaxSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelYmax;
                            
                            spkmap2dslopeXY_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelslopeXY;
                            spkmap2dphi0XY_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelphi0XY;
                            spkmap2drhoXY_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelrhoXY;
                            spkmap2dslopeXYSE_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelslopeXY;
                            spkmap2dphi0XYSE_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelphi0XY;
                            spkmap2drhoXYSE_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelrhoXY;
                            spkmap2dslopeXYiter_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).respModelslopeXY;
                            spkmap2dphi0XYiter_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).respModelphi0XY;
                            spkmap2drhoXYiter_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).respModelrhoXY;
                            
                            spkmap2dXcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelX;
                            spkmap2dXcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelX;
                            spkmap2dXPoscorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXpos;
                            spkmap2dXPoscorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXpos;
                            spkmap2dXmaxcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmax;
                            spkmap2dXmaxcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmax;
                            spkmap2dXmaxAmpcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmaxAmp;
                            spkmap2dXmaxAmpcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxAmp;
                            spkmap2dXmaxAmpcorriter_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).corrModelXmaxAmp;
                            spkmap2dXmaxsinAmpcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmaxsinAmp;
                            spkmap2dXmaxsinAmpcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxsinAmp;
                            spkmap2dXmaxsinAmpcorriter_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).corrModelXmaxsinAmp;
                            spkmap2dXmaxOffsetcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXmaxOffset;
                            spkmap2dXmaxOffsetcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXmaxOffset;
                            spkmap2dXmaxOffsetcorriter_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).corrModelXmaxOffset;
                            
                            spkmap2dXRefcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRef;
                            spkmap2dXRefcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRef;
                            spkmap2dXRefPoscorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefpos;
                            spkmap2dXRefPoscorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefpos;
                            spkmap2dXRefmaxcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmax;
                            spkmap2dXRefmaxcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmax;
                            spkmap2dXRefmaxAmpcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxAmp;
                            spkmap2dXRefmaxAmpcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxAmp;
                            spkmap2dXRefmaxAmpcorriter_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).corrModelXRefmaxAmp;
                            spkmap2dXRefmaxsinAmpcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxsinAmp;
                            spkmap2dXRefmaxsinAmpcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxsinAmp;
                            spkmap2dXRefmaxsinAmpcorriter_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).corrModelXRefmaxsinAmp;
                            spkmap2dXRefmaxOffsetcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelXRefmaxOffset;
                            spkmap2dXRefmaxOffsetcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelXRefmaxOffset;
                            spkmap2dXRefmaxOffsetcorriter_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).corrModelXRefmaxOffset;
                            
                            spkmap2dYcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelY;
                            spkmap2dYcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelY;
                            spkmap2dYPoscorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelYpos;
                            spkmap2dYPoscorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelYpos;
                            spkmap2dYmaxcorr_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).meancorrModelYmax;
                            spkmap2dYmaxcorrSE_theta_Shf = obj.maps2d.spiketrajPercent_spikePhase_Shf{c, g, r, o}.model.tuning(icell).SEcorrModelYmax;
                        end
                    end
                    spkmap2d_theta_2fold = [];
                    if isfield(obj.maps2d,'spiketrajPercent_spikePhase_2fold')
                        if ~isempty(obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model)
                            spkmap2dXcorr_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).corrModelX;
                            spkmap2dXPoscorr_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).corrModelXpos;
                            spkmap2dXmaxcorr_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).corrModelXmax;
                            
                            spkmap2dYcorr_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).corrModelY;
                            spkmap2dYPoscorr_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).corrModelYpos;
                            spkmap2dYmaxcorr_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).corrModelYmax;
                            
                            spkmap2d_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModel;
                            spkmap2dX_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelX;
                            spkmap2dXPos_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelXpos;
                            spkmap2dXmax_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelXmax;
                            
                            spkmap2dY_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelY;
                            spkmap2dYPos_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelYpos;
                            spkmap2dYmax_theta_2fold = obj.maps2d.spiketrajPercent_spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelYmax;
                        end
                    end
                    if ~isempty(spkmap2d_theta)
                        obj.CellInfo.spkfield2dXPhstheta{c, g, r, o}(icell,:,:) = spkmap2d_theta;
                        obj.CellInfo.spkfield2dXPhsthetaSE{c, g, r, o}(icell,:,:) = spkmap2dSE_theta;
                        obj.CellInfo.spkfield2dXtheta{c, g, r, o}(icell,:) = spkmap2dX_theta(:)';
                        obj.CellInfo.spkfield2dXthetaSE{c, g, r, o}(icell,:) = spkmap2dXSE_theta(:)';
                        obj.CellInfo.spkfield2dXthetapos{c, g, r, o}(icell,:) = spkmap2dXPos_theta(:)';
                        obj.CellInfo.spkfield2dXthetaposSE{c, g, r, o}(icell,:) = spkmap2dXPosSE_theta(:)';
                        obj.CellInfo.spkfield2dXthetamax{c, g, r, o}(icell,:) = spkmap2dXmax_theta(:)';
                        obj.CellInfo.spkfield2dXthetamaxSE{c, g, r, o}(icell,:) = spkmap2dXmaxSE_theta(:)';
                        
                        obj.CellInfo.spkfield2dPhstheta{c, g, r, o}(icell,:) = spkmap2dY_theta(:)';
                        obj.CellInfo.spkfield2dPhsthetaSE{c, g, r, o}(icell,:) = spkmap2dYSE_theta(:)';
                        obj.CellInfo.spkfield2dPhsthetapos{c, g, r, o}(icell,:) = spkmap2dYPos_theta(:)';
                        obj.CellInfo.spkfield2dPhsthetaposSE{c, g, r, o}(icell,:) = spkmap2dYPosSE_theta(:)';
                        obj.CellInfo.spkfield2dPhsthetamax{c, g, r, o}(icell,:) = spkmap2dYmax_theta(:)';
                        obj.CellInfo.spkfield2dPhsthetamaxSE{c, g, r, o}(icell,:) = spkmap2dYmaxSE_theta(:)';
                        
                        obj.CellInfo.spkfield2dslopeXY{c, g, r, o}(icell) = spkmap2dslopeXY;
                        obj.CellInfo.spkfield2dphi0XY{c, g, r, o}(icell) = spkmap2dphi0XY;
                        obj.CellInfo.spkfield2drhoXY{c, g, r, o}(icell) = spkmap2drhoXY;
                        obj.CellInfo.spkfield2dslopeXYSE{c, g, r, o}(icell) = spkmap2dslopeXYSE;
                        obj.CellInfo.spkfield2dphi0XYSE{c, g, r, o}(icell) = spkmap2dphi0XYSE;
                        obj.CellInfo.spkfield2drhoXYSE{c, g, r, o}(icell) = spkmap2drhoXYSE;
                        obj.CellInfo.spkfield2dslopeXYiter{c, g, r, o}(icell,:) = spkmap2dslopeXYiter;
                        obj.CellInfo.spkfield2dphi0XYiter{c, g, r, o}(icell,:) = spkmap2dphi0XYiter;
                        obj.CellInfo.spkfield2drhoXYiter{c, g, r, o}(icell,:) = spkmap2drhoXYiter;
                        
                        obj.CellInfo.spkfield2dXcorrtheta{c, g, r, o}(icell,:,:) = spkmap2dXcorr_theta;
                        obj.CellInfo.spkfield2dXcorrthetaSE{c, g, r, o}(icell,:,:) = spkmap2dXcorrSE_theta;
                        obj.CellInfo.spkfield2dXcorrthetapos{c, g, r, o}(icell,:) = spkmap2dXPoscorr_theta;
                        obj.CellInfo.spkfield2dXcorrthetaposSE{c, g, r, o}(icell,:) = spkmap2dXPoscorrSE_theta;
                        obj.CellInfo.spkfield2dXcorrthetamax{c, g, r, o}(icell,:) = spkmap2dXmaxcorr_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxSE{c, g, r, o}(icell,:) = spkmap2dXmaxcorrSE_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxAmp{c, g, r, o}(icell) = spkmap2dXmaxAmpcorr_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxAmpSE{c, g, r, o}(icell) = spkmap2dXmaxAmpcorrSE_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxAmpiter{c, g, r, o}(icell,:) = spkmap2dXmaxAmpcorriter_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxsinAmp{c, g, r, o}(icell) = spkmap2dXmaxsinAmpcorr_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxsinAmpSE{c, g, r, o}(icell) = spkmap2dXmaxsinAmpcorrSE_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = spkmap2dXmaxsinAmpcorriter_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxOffset{c, g, r, o}(icell) = spkmap2dXmaxOffsetcorr_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxOffsetSE{c, g, r, o}(icell) = spkmap2dXmaxOffsetcorrSE_theta;
                        obj.CellInfo.spkfield2dXcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = spkmap2dXmaxOffsetcorriter_theta;
                        
                        obj.CellInfo.spkfield2dXRefcorrtheta{c, g, r, o}(icell,:,:) = spkmap2dXRefcorr_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetaSE{c, g, r, o}(icell,:,:) = spkmap2dXRefcorrSE_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetapos{c, g, r, o}(icell,:) = spkmap2dXRefPoscorr_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetaposSE{c, g, r, o}(icell,:) = spkmap2dXRefPoscorrSE_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamax{c, g, r, o}(icell,:) = spkmap2dXRefmaxcorr_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxSE{c, g, r, o}(icell,:) = spkmap2dXRefmaxcorrSE_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxAmp{c, g, r, o}(icell) = spkmap2dXRefmaxAmpcorr_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxAmpSE{c, g, r, o}(icell) = spkmap2dXRefmaxAmpcorrSE_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxAmpiter{c, g, r, o}(icell,:) = spkmap2dXRefmaxAmpcorriter_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmp{c, g, r, o}(icell) = spkmap2dXRefmaxsinAmpcorr_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmpSE{c, g, r, o}(icell) = spkmap2dXRefmaxsinAmpcorrSE_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = spkmap2dXRefmaxsinAmpcorriter_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxOffset{c, g, r, o}(icell) = spkmap2dXRefmaxOffsetcorr_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxOffsetSE{c, g, r, o}(icell) = spkmap2dXRefmaxOffsetcorrSE_theta;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = spkmap2dXRefmaxOffsetcorriter_theta;
                        
                        obj.CellInfo.spkfield2dPhscorrtheta{c, g, r, o}(icell,:,:) = spkmap2dYcorr_theta;
                        obj.CellInfo.spkfield2dPhscorrthetaSE{c, g, r, o}(icell,:,:) = spkmap2dYcorrSE_theta;
                        obj.CellInfo.spkfield2dPhscorrthetapos{c, g, r, o}(icell,:) = spkmap2dYPoscorr_theta;
                        obj.CellInfo.spkfield2dPhscorrthetaposSE{c, g, r, o}(icell,:) = spkmap2dYPoscorrSE_theta;
                        obj.CellInfo.spkfield2dPhscorrthetamax{c, g, r, o}(icell,:) = spkmap2dYmaxcorr_theta;
                        obj.CellInfo.spkfield2dPhscorrthetamaxSE{c, g, r, o}(icell,:) = spkmap2dYmaxcorrSE_theta;
                        
                        obj.CellInfo.spkPhaseRho{c, g, r, o}(icell) = spkrho;
                        obj.CellInfo.spkPhasePval{c, g, r, o}(icell) = spkpval;
                        obj.CellInfo.spkPhaseRayleighPval{c, g, r, o}(icell) = spkrtest_pval;
                        obj.CellInfo.spkPhaseRayleighZ{c, g, r, o}(icell) = spkrtest_z;
                        obj.CellInfo.spkPhaseMax{c, g, r, o}(icell) = spkphsmax;
                    else
                        obj.CellInfo.spkfield2dXPhstheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXPhsthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXtheta{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetamaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dPhstheta{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetamaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dslopeXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dphi0XY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2drhoXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dslopeXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dphi0XYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2drhoXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dslopeXYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dphi0XYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2drhoXYiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dXcorrtheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXcorrthetamaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamaxOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dXRefcorrtheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxOffsetSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetamaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dPhscorrtheta{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetapos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetamax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetamaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkPhaseRho{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkPhasePval{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkPhaseRayleighPval{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkPhaseRayleighZ{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkPhaseMax{c, g, r, o}(icell) = NaN;
                    end
                    
                    if ~isempty(spkmap2d_theta_Shf)
                        obj.CellInfo.spkfield2dXPhsthetaShf{c, g, r, o}(icell,:,:) = spkmap2d_theta_Shf;
                        obj.CellInfo.spkfield2dXPhsthetaShfSE{c, g, r, o}(icell,:,:) = spkmap2dSE_theta_Shf;
                        obj.CellInfo.spkfield2dXthetaShf{c, g, r, o}(icell,:) = spkmap2dX_theta_Shf(:)';
                        obj.CellInfo.spkfield2dXthetaShfSE{c, g, r, o}(icell,:) = spkmap2dXSE_theta_Shf(:)';
                        obj.CellInfo.spkfield2dXthetaShfpos{c, g, r, o}(icell,:) = spkmap2dXPos_theta_Shf(:)';
                        obj.CellInfo.spkfield2dXthetaShfposSE{c, g, r, o}(icell,:) = spkmap2dXPosSE_theta_Shf(:)';
                        obj.CellInfo.spkfield2dXthetaShfmax{c, g, r, o}(icell,:) = spkmap2dXmax_theta_Shf(:)';
                        obj.CellInfo.spkfield2dXthetaShfmaxSE{c, g, r, o}(icell,:) = spkmap2dXmaxSE_theta_Shf(:)';
                        
                        obj.CellInfo.spkfield2dPhsthetaShf{c, g, r, o}(icell,:) = spkmap2dY_theta_Shf(:)';
                        obj.CellInfo.spkfield2dPhsthetaShfSE{c, g, r, o}(icell,:) = spkmap2dYSE_theta_Shf(:)';
                        obj.CellInfo.spkfield2dPhsthetaShfpos{c, g, r, o}(icell,:) = spkmap2dYPos_theta_Shf(:)';
                        obj.CellInfo.spkfield2dPhsthetaShfposSE{c, g, r, o}(icell,:) = spkmap2dYPosSE_theta_Shf(:)';
                        obj.CellInfo.spkfield2dPhsthetaShfmax{c, g, r, o}(icell,:) = spkmap2dYmax_theta_Shf(:)';
                        obj.CellInfo.spkfield2dPhsthetaShfmaxSE{c, g, r, o}(icell,:) = spkmap2dYmaxSE_theta_Shf(:)';
                        
                        obj.CellInfo.spkfield2dShfslopeXY{c, g, r, o}(icell) = spkmap2dslopeXY_Shf;
                        obj.CellInfo.spkfield2dShfphi0XY{c, g, r, o}(icell) = spkmap2dphi0XY_Shf;
                        obj.CellInfo.spkfield2dShfrhoXY{c, g, r, o}(icell) = spkmap2drhoXY_Shf;
                        obj.CellInfo.spkfield2dShfslopeXYSE{c, g, r, o}(icell) = spkmap2dslopeXYSE_Shf;
                        obj.CellInfo.spkfield2dShfphi0XYSE{c, g, r, o}(icell) = spkmap2dphi0XYSE_Shf;
                        obj.CellInfo.spkfield2dShfrhoXYSE{c, g, r, o}(icell) = spkmap2drhoXYSE_Shf;
                        obj.CellInfo.spkfield2dShfslopeXYiter{c, g, r, o}(icell,:) = spkmap2dslopeXYiter_Shf;
                        obj.CellInfo.spkfield2dShfphi0XYiter{c, g, r, o}(icell,:) = spkmap2dphi0XYiter_Shf;
                        obj.CellInfo.spkfield2dShfrhoXYiter{c, g, r, o}(icell,:) = spkmap2drhoXYiter_Shf;
                        
                        obj.CellInfo.spkfield2dXcorrthetaShf{c, g, r, o}(icell,:,:) = spkmap2dXcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfSE{c, g, r, o}(icell,:,:) = spkmap2dXcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfpos{c, g, r, o}(icell,:) = spkmap2dXPoscorr_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfposSE{c, g, r, o}(icell,:) = spkmap2dXPoscorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmax{c, g, r, o}(icell,:) = spkmap2dXmaxcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxSE{c, g, r, o}(icell,:) = spkmap2dXmaxcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxAmp{c, g, r, o}(icell) = spkmap2dXmaxAmpcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = spkmap2dXmaxAmpcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = spkmap2dXmaxAmpcorriter_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = spkmap2dXmaxsinAmpcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = spkmap2dXmaxsinAmpcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = spkmap2dXmaxsinAmpcorriter_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxOffset{c, g, r, o}(icell) = spkmap2dXmaxOffsetcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = spkmap2dXmaxOffsetcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = spkmap2dXmaxOffsetcorriter_theta_Shf;
                        
                        obj.CellInfo.spkfield2dXRefcorrthetaShf{c, g, r, o}(icell,:,:) = spkmap2dXRefcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfSE{c, g, r, o}(icell,:,:) = spkmap2dXRefcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfpos{c, g, r, o}(icell,:) = spkmap2dXRefPoscorr_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfposSE{c, g, r, o}(icell,:) = spkmap2dXRefPoscorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmax{c, g, r, o}(icell,:) = spkmap2dXRefmaxcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxSE{c, g, r, o}(icell,:) = spkmap2dXRefmaxcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmp{c, g, r, o}(icell) = spkmap2dXRefmaxAmpcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = spkmap2dXRefmaxAmpcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = spkmap2dXRefmaxAmpcorriter_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = spkmap2dXRefmaxsinAmpcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = spkmap2dXRefmaxsinAmpcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = spkmap2dXRefmaxsinAmpcorriter_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffset{c, g, r, o}(icell) = spkmap2dXRefmaxOffsetcorr_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = spkmap2dXRefmaxOffsetcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = spkmap2dXRefmaxOffsetcorriter_theta_Shf;
                        
                        obj.CellInfo.spkfield2dPhscorrthetaShf{c, g, r, o}(icell,:,:) = spkmap2dYcorr_theta_Shf;
                        obj.CellInfo.spkfield2dPhscorrthetaShfSE{c, g, r, o}(icell,:,:) = spkmap2dYcorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dPhscorrthetaShfpos{c, g, r, o}(icell,:) = spkmap2dYPoscorr_theta_Shf;
                        obj.CellInfo.spkfield2dPhscorrthetaShfposSE{c, g, r, o}(icell,:) = spkmap2dYPoscorrSE_theta_Shf;
                        obj.CellInfo.spkfield2dPhscorrthetaShfmax{c, g, r, o}(icell,:) = spkmap2dYmaxcorr_theta_Shf;
                        obj.CellInfo.spkfield2dPhscorrthetaShfmaxSE{c, g, r, o}(icell,:) = spkmap2dYmaxcorrSE_theta_Shf;                   
                    else
                        obj.CellInfo.spkfield2dXPhsthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXPhsthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXthetaShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dPhsthetaShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhsthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dShfslopeXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dShfphi0XY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dShfrhoXY{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dShfslopeXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dShfphi0XYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dShfrhoXYSE{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dShfslopeXYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dShfphi0XYiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dShfrhoXYiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dXcorrthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dXRefcorrthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmp{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxsinAmpiter{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffset{c, g, r, o}(icell) = 0;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkfield2dXRefcorrthetaShfmaxOffsetiter{c, g, r, o}(icell,:) = 0;
                        
                        obj.CellInfo.spkfield2dPhscorrthetaShf{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaShfSE{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaShfpos{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaShfposSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaShfmax{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetaShfmaxSE{c, g, r, o}(icell,:) = 0;
                    end
                    
                    if ~isempty(spkmap2d_theta_2fold)
                        obj.CellInfo.spkfield2dXcorrtheta_2fold{c, g, r, o}(icell,:,:,:) = spkmap2dXcorr_theta_2fold;
                        obj.CellInfo.spkfield2dXcorrthetapos_2fold{c, g, r, o}(icell,:,:) = spkmap2dXPoscorr_theta_2fold;
                        obj.CellInfo.spkfield2dXcorrthetamax_2fold{c, g, r, o}(icell,:,:) = spkmap2dXmaxcorr_theta_2fold;
                        
                        obj.CellInfo.spkfield2dPhscorrtheta_2fold{c, g, r, o}(icell,:,:,:) = spkmap2dYcorr_theta_2fold;
                        obj.CellInfo.spkfield2dPhscorrthetapos_2fold{c, g, r, o}(icell,:,:) = spkmap2dYPoscorr_theta_2fold;
                        obj.CellInfo.spkfield2dPhscorrthetamax_2fold{c, g, r, o}(icell,:,:) = spkmap2dYmaxcorr_theta_2fold;
                        
                        obj.CellInfo.spkfield2dXPhstheta_2fold{c, g, r, o}(icell,:,:,:) = spkmap2d_theta_2fold;
                        obj.CellInfo.spkfield2dXtheta_2fold{c, g, r, o}(icell,:,:) = spkmap2dX_theta_2fold;
                        obj.CellInfo.spkfield2dXthetapos_2fold{c, g, r, o}(icell,:,:) = spkmap2dXPos_theta_2fold;
                        obj.CellInfo.spkfield2dXthetamax_2fold{c, g, r, o}(icell,:,:) = spkmap2dXmax_theta_2fold;
                        
                        obj.CellInfo.spkfield2dYtheta_2fold{c, g, r, o}(icell,:,:) = spkmap2dY_theta_2fold;
                        obj.CellInfo.spkfield2dYthetapos_2fold{c, g, r, o}(icell,:,:) = spkmap2dYPos_theta_2fold;
                        obj.CellInfo.spkfield2dYthetamax_2fold{c, g, r, o}(icell,:,:) = spkmap2dYmax_theta_2fold;
                    else
                        obj.CellInfo.spkfield2dXcorrtheta_2fold{c, g, r, o}(icell,:,:,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXcorrthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                        
                        obj.CellInfo.spkfield2dPhscorrtheta_2fold{c, g, r, o}(icell,:,:,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dPhscorrthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                        
                        obj.CellInfo.spkfield2dXPhstheta_2fold{c, g, r, o}(icell,:,:,:) = 0;
                        obj.CellInfo.spkfield2dXtheta_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dXthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                        
                        obj.CellInfo.spkfield2dYtheta_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dYthetapos_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkfield2dYthetamax_2fold{c, g, r, o}(icell,:,:) = 0;
                    end
                    
                    phsmap = [];
                    if isfield(obj.maps1d,'LFPphase2')
                        if ~isempty(obj.maps1d.LFPphase2{c, g, r, o}.model)
                            phsmap = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModel;
                            phsmapSE = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).SErespModel;
                            phsfieldCOM = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            phsfieldCOMSE = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            phsfieldmax = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            phsfieldmaxSE = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            phsfieldAmp = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            phsfieldAmpSE = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            phsfieldAmpiter = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            phsfieldsinAmp = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXsinAmp;
                            phsfieldsinAmpSE = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXsinAmp;
                            phsfieldsinAmpiter = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).respModelXsinAmp;
                            phsfieldsinOffset = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).meanrespModelXsinOffset;
                            phsfieldsinOffsetSE = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).SErespModelXsinOffset;
                            phsfieldsinOffsetiter = obj.maps1d.LFPphase2{c, g, r, o}.model.tuning(icell).respModelXsinOffset;
                        end
                    end
                    phsmap_Shf = [];
                    if isfield(obj.maps1d,'LFPphase2_Shf')
                        if ~isempty(obj.maps1d.LFPphase2_Shf{c, g, r, o}.model)
                            phsmap_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModel;
                            phsmapSE_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModel;
                            phsfieldCOM_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            phsfieldCOMSE_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            phsfieldmax_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            phsfieldmaxSE_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            phsfieldAmp_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            phsfieldAmpSE_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            phsfieldAmpiter_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            phsfieldsinAmp_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXsinAmp;
                            phsfieldsinAmpSE_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXsinAmp;
                            phsfieldsinAmpiter_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).respModelXsinAmp;
                            phsfieldsinOffset_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXsinOffset;
                            phsfieldsinOffsetSE_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).SErespModelXsinOffset;
                            phsfieldsinOffsetiter_Shf = obj.maps1d.LFPphase2_Shf{c, g, r, o}.model.tuning(icell).respModelXsinOffset;
                        end
                    end
                    phsmap_2fold = [];
                    if isfield(obj.maps1d,'LFPphase2_2fold')
                        if ~isempty(obj.maps1d.LFPphase2_2fold{c, g, r, o}.model)
                            phsmap_2fold = obj.maps1d.LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModel;
                            phsfieldCOM_2fold = obj.maps1d.LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelXpos;
                            phsfieldmax_2fold = obj.maps1d.LFPphase2_2fold{c, g, r, o}.model.tuning(icell).respModelXmax;
                        end
                    end
                    if ~isempty(phsmap)
                        obj.CellInfo.phsfieldPos{c, g, r, o}(icell) = phsfieldmax;
                        obj.CellInfo.phsfieldPosSE{c, g, r, o}(icell) = phsfieldmaxSE;
                        obj.CellInfo.phsfieldCOM{c, g, r, o}(icell) = phsfieldCOM;
                        obj.CellInfo.phsfieldCOMSE{c, g, r, o}(icell) = phsfieldCOMSE;
                        obj.CellInfo.phsfield{c, g, r, o}(icell,:) = phsmap;
                        obj.CellInfo.phsfieldSE{c, g, r, o}(icell,:) = phsmapSE;
                        obj.CellInfo.phsfieldAmp{c, g, r, o}(icell) = phsfieldAmp;
                        obj.CellInfo.phsfieldAmpSE{c, g, r, o}(icell) = phsfieldAmpSE;
                        obj.CellInfo.phsfieldAmpiter{c, g, r, o}(icell,:) = phsfieldAmpiter;
                        obj.CellInfo.phsfieldsinAmp{c, g, r, o}(icell) = phsfieldsinAmp;
                        obj.CellInfo.phsfieldsinAmpSE{c, g, r, o}(icell) = phsfieldsinAmpSE;
                        obj.CellInfo.phsfieldsinAmpiter{c, g, r, o}(icell,:) = phsfieldsinAmpiter;
                        obj.CellInfo.phsfieldsinOffset{c, g, r, o}(icell) = phsfieldsinOffset;
                        obj.CellInfo.phsfieldsinOffsetSE{c, g, r, o}(icell) = phsfieldsinOffsetSE;
                        obj.CellInfo.phsfieldsinOffsetiter{c, g, r, o}(icell,:) = phsfieldsinOffsetiter;
                        obj.CellInfo.phsmodulation{c, g, r, o}(icell) = (max(phsmap) - mean(phsmap))/mean(phsmap);
                    else
                        obj.CellInfo.phsfieldPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldCOM{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldCOMSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfieldsinAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldsinAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfieldsinOffset{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldsinOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldsinOffsetiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfield{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.phsfieldSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.phsmodulation{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(phsmap_Shf)
                        obj.CellInfo.phsfieldShfPos{c, g, r, o}(icell) = phsfieldmax_Shf;
                        obj.CellInfo.phsfieldShfPosSE{c, g, r, o}(icell) = phsfieldmaxSE_Shf;
                        obj.CellInfo.phsfieldShfCOM{c, g, r, o}(icell) = phsfieldCOM_Shf;
                        obj.CellInfo.phsfieldShfCOMSE{c, g, r, o}(icell) = phsfieldCOMSE_Shf;
                        obj.CellInfo.phsfieldShf{c, g, r, o}(icell,:) = phsmap_Shf;
                        obj.CellInfo.phsfieldShfSE{c, g, r, o}(icell,:) = phsmapSE_Shf;
                        obj.CellInfo.phsfieldShfAmp{c, g, r, o}(icell) = phsfieldAmp_Shf;
                        obj.CellInfo.phsfieldShfAmpSE{c, g, r, o}(icell) = phsfieldAmpSE_Shf;
                        obj.CellInfo.phsfieldShfAmpiter{c, g, r, o}(icell,:) = phsfieldAmpiter_Shf;
                        obj.CellInfo.phsfieldShfsinAmp{c, g, r, o}(icell) = phsfieldsinAmp_Shf;
                        obj.CellInfo.phsfieldShfsinAmpSE{c, g, r, o}(icell) = phsfieldsinAmpSE_Shf;
                        obj.CellInfo.phsfieldShfsinAmpiter{c, g, r, o}(icell,:) = phsfieldsinAmpiter_Shf;
                        obj.CellInfo.phsfieldShfsinOffset{c, g, r, o}(icell) = phsfieldsinOffset_Shf;
                        obj.CellInfo.phsfieldShfsinOffsetSE{c, g, r, o}(icell) = phsfieldsinOffsetSE_Shf;
                        obj.CellInfo.phsfieldShfsinOffsetiter{c, g, r, o}(icell,:) = phsfieldsinOffsetiter_Shf;
                        obj.CellInfo.phsmodulationShf{c, g, r, o}(icell) = (max(phsmap_Shf) - mean(phsmap_Shf))/mean(phsmap_Shf);
                    else
                        obj.CellInfo.phsfieldShfPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfCOM{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfCOMSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfieldShfsinAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfsinAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfieldShfsinOffset{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfsinOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.phsfieldShfsinOffsetiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfieldShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.phsfieldShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.phsmodulationShf{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(phsmap_2fold)
                        obj.CellInfo.phsfield_2fold{c, g, r, o}(icell,:,:) = phsmap_2fold;
                        obj.CellInfo.phsfieldCOM_2fold{c, g, r, o}(icell,:) = phsfieldCOM_2fold;
                        obj.CellInfo.phsfieldPos_2fold{c, g, r, o}(icell,:) = phsfieldmax_2fold;
                    else
                        obj.CellInfo.phsfield_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.phsfieldCOM_2fold{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.phsfieldPos_2fold{c, g, r, o}(icell,:) = NaN;
                    end
                    
                    %
                    spkphsmap = [];
                    if isfield(obj.maps1d,'spikePhase')
                        if ~isempty(obj.maps1d.spikePhase{c, g, r, o}.model)
                            spkphsmap = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).meanrespModel;
                            spkphsmapSE = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).SErespModel;
                            spkphsfieldCOM = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            spkphsfieldCOMSE = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            spkphsfieldmax = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            spkphsfieldmaxSE = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            spkphsfieldAmp = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            spkphsfieldAmpSE = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            spkphsfieldAmpiter = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            spkphsfieldsinAmp = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXsinAmp;
                            spkphsfieldsinAmpSE = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXsinAmp;
                            spkphsfieldsinAmpiter = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).respModelXsinAmp;
                            spkphsfieldsinOffset = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).meanrespModelXsinOffset;
                            spkphsfieldsinOffsetSE = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).SErespModelXsinOffset;
                            spkphsfieldsinOffsetiter = obj.maps1d.spikePhase{c, g, r, o}.model.tuning(icell).respModelXsinOffset;
                        end
                    end
                    spkphsmap_Shf = [];
                    if isfield(obj.maps1d,'spikePhase_Shf')
                        if ~isempty(obj.maps1d.spikePhase_Shf{c, g, r, o}.model)
                            spkphsmap_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModel;
                            spkphsmapSE_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModel;
                            spkphsfieldCOM_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXpos;
                            spkphsfieldCOMSE_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXpos;
                            spkphsfieldmax_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXmax;
                            spkphsfieldmaxSE_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXmax;
                            spkphsfieldAmp_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXAmp;
                            spkphsfieldAmpSE_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXAmp;
                            spkphsfieldAmpiter_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).respModelXAmp;
                            spkphsfieldsinAmp_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXsinAmp;
                            spkphsfieldsinAmpSE_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXsinAmp;
                            spkphsfieldsinAmpiter_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).respModelXsinAmp;
                            spkphsfieldsinOffset_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).meanrespModelXsinOffset;
                            spkphsfieldsinOffsetSE_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).SErespModelXsinOffset;
                            spkphsfieldsinOffsetiter_Shf = obj.maps1d.spikePhase_Shf{c, g, r, o}.model.tuning(icell).respModelXsinOffset;
                        end
                    end
                    spkphsmap_2fold = [];
                    if isfield(obj.maps1d,'spikePhase_2fold')
                        if ~isempty(obj.maps1d.spikePhase_2fold{c, g, r, o}.model)
                            spkphsmap_2fold = obj.maps1d.spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModel;
                            spkphsfieldCOM_2fold = obj.maps1d.spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelXpos;
                            spkphsfieldmax_2fold = obj.maps1d.spikePhase_2fold{c, g, r, o}.model.tuning(icell).respModelXmax;
                        end
                    end
                    if ~isempty(spkphsmap)
                        obj.CellInfo.spkphsfieldPos{c, g, r, o}(icell) = spkphsfieldmax;
                        obj.CellInfo.spkphsfieldPosSE{c, g, r, o}(icell) = spkphsfieldmaxSE;
                        obj.CellInfo.spkphsfieldCOM{c, g, r, o}(icell) = spkphsfieldCOM;
                        obj.CellInfo.spkphsfieldCOMSE{c, g, r, o}(icell) = spkphsfieldCOMSE;
                        obj.CellInfo.spkphsfield{c, g, r, o}(icell,:) = spkphsmap;
                        obj.CellInfo.spkphsfieldSE{c, g, r, o}(icell,:) = spkphsmapSE;
                        obj.CellInfo.spkphsfieldAmp{c, g, r, o}(icell) = spkphsfieldAmp;
                        obj.CellInfo.spkphsfieldAmpSE{c, g, r, o}(icell) = spkphsfieldAmpSE;
                        obj.CellInfo.spkphsfieldAmpiter{c, g, r, o}(icell,:) = spkphsfieldAmpiter;
                        obj.CellInfo.spkphsfieldsinAmp{c, g, r, o}(icell) = spkphsfieldsinAmp;
                        obj.CellInfo.spkphsfieldsinAmpSE{c, g, r, o}(icell) = spkphsfieldsinAmpSE;
                        obj.CellInfo.spkphsfieldsinAmpiter{c, g, r, o}(icell,:) = spkphsfieldsinAmpiter;
                        obj.CellInfo.spkphsfieldsinOffset{c, g, r, o}(icell) = spkphsfieldsinOffset;
                        obj.CellInfo.spkphsfieldsinOffsetSE{c, g, r, o}(icell) = spkphsfieldsinOffsetSE;
                        obj.CellInfo.spkphsfieldsinOffsetiter{c, g, r, o}(icell,:) = spkphsfieldsinOffsetiter;
                        obj.CellInfo.spkphsmodulation{c, g, r, o}(icell) = (max(spkphsmap) - mean(spkphsmap))/mean(spkphsmap);
                    else
                        obj.CellInfo.spkphsfieldPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldCOM{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldCOMSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfieldsinAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldsinAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfieldsinOffset{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldsinOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldsinOffsetiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfield{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkphsfieldSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkphsmodulation{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(spkphsmap_Shf)
                        obj.CellInfo.spkphsfieldShfPos{c, g, r, o}(icell) = spkphsfieldmax_Shf;
                        obj.CellInfo.spkphsfieldShfPosSE{c, g, r, o}(icell) = spkphsfieldmaxSE_Shf;
                        obj.CellInfo.spkphsfieldShfCOM{c, g, r, o}(icell) = spkphsfieldCOM_Shf;
                        obj.CellInfo.spkphsfieldShfCOMSE{c, g, r, o}(icell) = spkphsfieldCOMSE_Shf;
                        obj.CellInfo.spkphsfieldShf{c, g, r, o}(icell,:) = spkphsmap_Shf;
                        obj.CellInfo.spkphsfieldShfSE{c, g, r, o}(icell,:) = spkphsmapSE_Shf;
                        obj.CellInfo.spkphsfieldShfAmp{c, g, r, o}(icell) = spkphsfieldAmp_Shf;
                        obj.CellInfo.spkphsfieldShfAmpSE{c, g, r, o}(icell) = spkphsfieldAmpSE_Shf;
                        obj.CellInfo.spkphsfieldShfAmpiter{c, g, r, o}(icell,:) = spkphsfieldAmpiter_Shf;
                        obj.CellInfo.spkphsfieldShfsinAmp{c, g, r, o}(icell) = spkphsfieldsinAmp_Shf;
                        obj.CellInfo.spkphsfieldShfsinAmpSE{c, g, r, o}(icell) = spkphsfieldsinAmpSE_Shf;
                        obj.CellInfo.spkphsfieldShfsinAmpiter{c, g, r, o}(icell,:) = spkphsfieldsinAmpiter_Shf;
                        obj.CellInfo.spkphsfieldShfsinOffset{c, g, r, o}(icell) = spkphsfieldsinOffset_Shf;
                        obj.CellInfo.spkphsfieldShfsinOffsetSE{c, g, r, o}(icell) = spkphsfieldsinOffsetSE_Shf;
                        obj.CellInfo.spkphsfieldShfsinOffsetiter{c, g, r, o}(icell,:) = spkphsfieldsinOffsetiter_Shf;
                        obj.CellInfo.spkphsmodulationShf{c, g, r, o}(icell) = (max(spkphsmap_Shf) - mean(spkphsmap_Shf))/mean(spkphsmap_Shf);
                    else
                        obj.CellInfo.spkphsfieldShfPos{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfPosSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfCOM{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfCOMSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfieldShfsinAmp{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfsinAmpSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfsinAmpiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfieldShfsinOffset{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfsinOffsetSE{c, g, r, o}(icell) = NaN;
                        obj.CellInfo.spkphsfieldShfsinOffsetiter{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfieldShf{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkphsfieldShfSE{c, g, r, o}(icell,:) = 0;
                        obj.CellInfo.spkphsmodulationShf{c, g, r, o}(icell) = NaN;
                    end
                    if ~isempty(spkphsmap_2fold)
                        obj.CellInfo.spkphsfield_2fold{c, g, r, o}(icell,:,:) = spkphsmap_2fold;
                        obj.CellInfo.spkphsfieldCOM_2fold{c, g, r, o}(icell,:) = spkphsfieldCOM_2fold;
                        obj.CellInfo.spkphsfieldPos_2fold{c, g, r, o}(icell,:) = spkphsfieldmax_2fold;
                    else
                        obj.CellInfo.spkphsfield_2fold{c, g, r, o}(icell,:,:) = 0;
                        obj.CellInfo.spkphsfieldCOM_2fold{c, g, r, o}(icell,:) = NaN;
                        obj.CellInfo.spkphsfieldPos_2fold{c, g, r, o}(icell,:) = NaN;
                    end
                end
                
                
                if isfield(obj.data.es,'LFPphase')
                    if c == cbase && g == gbase && r == rbase && o == obase
                        tidx_phs = obj.data.es.smthBallSpd > obj.SpeedThresh & ~isnan(obj.data.es.smthBallSpd) & ~isnan(obj.data.es.smthTrajSpd);%tidx;
                        zphsfield_th = 0;
                        wfphsfield_th = 15;
                        for iprobe = 1:2
                            try
                            phscellidx = obj.CellInfo.Goodcluster & obj.CellInfo.Probe==iprobe & ~obj.CellInfo.Finterneuron & obj.CellInfo.min2maxSpkwf > wfphsfield_th;%& obj.CellInfo.phsfieldZ{c, g, r, o}>=zphsfield_th 
                            catch
                                keyboard
                            end
                            
                            [phsmapall,~,~,phsmapSEall] = fast1Dmap(obj.data.es.LFPphase(tidx_phs), sum(obj.data.es.spikeTrain(tidx_phs,phscellidx),2), 20, samplerate(tidx_phs),(360/20)/2,true);
                            obj.CellInfo.phsfieldMUA{iprobe} = phsmapall;
                            obj.CellInfo.phsfieldMUASE{iprobe} = phsmapSEall;
                            [~,imax] = max(phsmapall);
                            obj.CellInfo.phsfieldZMUA{iprobe} = (phsmapall(imax) - mean(phsmapall))./phsmapSEall(imax);
                            obj.CellInfo.phsfieldPosMUA{iprobe} = getCircularAverage(obj.CellInfo.phsfieldMUA{iprobe}(:),0,1);%getCircularAverage(obj.CellInfo.phsfieldMUA{iprobe}(:),0,0.01,0.05);
                            obj.CellInfo.LFP2Spike_phscorrMUA{iprobe} = obj.CellInfo.phsfieldPosMUA{iprobe}*(360/numel(obj.CellInfo.phsfieldMUA{iprobe})) - 180;
                            
                            phsgoodcells = find(phscellidx);
                            phsmapcell = NaN(numel(phsgoodcells),numel(phsmapall));
                            for icell = 1:numel(phsgoodcells)
                                [phsmapcell(icell,:),~,~,~] = fast1Dmap(obj.data.es.LFPphase(tidx_phs), obj.data.es.spikeTrain(tidx_phs,phsgoodcells(icell)), 20, samplerate(tidx_phs),(360/20)/2,true);
                            end
                            obj.CellInfo.phsfieldMUAnorm{iprobe} = nanmean(phsmapcell./repmat(nanmean(phsmapcell,2),[1 size(phsmapcell,2)]),1);
                            obj.CellInfo.phsfieldPosMUAnorm{iprobe} = getCircularAverage(obj.CellInfo.phsfieldMUAnorm{iprobe}(:),0,1);%getCircularAverage(obj.CellInfo.phsfieldMUAnorm{iprobe}(:),0,0.01,0.05);
                            obj.CellInfo.LFP2Spike_phscorrMUAnorm{iprobe} = obj.CellInfo.phsfieldPosMUAnorm{iprobe}*(360/numel(obj.CellInfo.phsfieldMUAnorm{iprobe})) - 180;
                            
                            obj.CellInfo.phsfieldMean{iprobe} = nanmean(obj.CellInfo.phsfield{c, g, r, o}(phscellidx,:),1);
                            obj.CellInfo.phsfieldPosMean{iprobe} = getCircularAverage(obj.CellInfo.phsfieldMean{iprobe}(:),0,1);%getCircularAverage(obj.CellInfo.phsfieldMean{iprobe}(:),0,0.01,0.05);
                            obj.CellInfo.LFP2Spike_phscorrMean{iprobe} = obj.CellInfo.phsfieldPosMean{iprobe}*(360/numel(obj.CellInfo.phsfieldMean{iprobe})) - 180;
                        end
                    end
                else
                    if c == cbase && g == gbase && r == rbase && o == obase
                        for iprobe = 1:2
                            obj.CellInfo.phsfieldMUA{iprobe} = 0;
                            obj.CellInfo.phsfieldMUASE{iprobe} = 0;
                            obj.CellInfo.phsfieldZMUA{iprobe} = NaN;
                            obj.CellInfo.phsfieldPosMUA{iprobe} = NaN;
                            obj.CellInfo.LFP2Spike_phscorrMUA{iprobe} = NaN;
                            
                            obj.CellInfo.phsfieldMUAnorm{iprobe} = 0;
                            obj.CellInfo.phsfieldPosMUAnorm{iprobe} = NaN;
                            obj.CellInfo.LFP2Spike_phscorrMUAnorm{iprobe} = NaN;
                                                        
                            obj.CellInfo.phsfieldMean{iprobe} = 0;
                            obj.CellInfo.phsfieldPosMean{iprobe} = NaN;
                            obj.CellInfo.LFP2Spike_phscorrMean{iprobe} = NaN;
                        end
                    end
                end
            end
        end
    end
end
end