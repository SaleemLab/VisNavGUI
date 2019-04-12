function EXP = CreateEXPStructure
%create and define fields of the data structure used troughout the GUI

EXP = Tstructure('Data');
EXP.addprop('animal');
EXP.addprop('series');
EXP.addprop('exp');

EXP.addprop('Nav');
EXP.addprop('Vis');
EXP.addprop('Eye');
EXP.addprop('Spk');
EXP.addprop('Lfp');

end
