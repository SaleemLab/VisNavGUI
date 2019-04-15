function         [prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,sorted_series] = ...
    eliminateUnit(prev,animal_spikesdir,prevSavedUnits,already_sorted_to_load,sorted_series,unitIDX)
% ELIMINATEUNIT.M: Moves unlucky unit to archive, updates sortnev2 variables (prevUnitProcess)

global DIRS

UnitArchive(prev(unitIDX).unit,DIRS.spikes);

% Update Variables to reflect elimination of unit
animal_spikesdir(prevSavedUnits(unitIDX):end-1) = animal_spikesdir(prevSavedUnits(unitIDX)+1:end);
animal_spikesdir = animal_spikesdir(1:end-1);

prevSavedUnits(prevSavedUnits > prevSavedUnits(unitIDX)) = ...
    prevSavedUnits(prevSavedUnits > prevSavedUnits(unitIDX)) + 1;
already_sorted_to_load(already_sorted_to_load(:,5) > already_sorted_to_load(unitIDX,5),5) = ...
    already_sorted_to_load(already_sorted_to_load(:,5) > already_sorted_to_load(unitIDX,5),5) + 1;

% already_sorted_to_load(unitIDX,:) = [];
%         prevSavedUnits(unitIDX  ) = [];
%          sorted_series(unitIDX  ) = [];
% TODO: change size of prev
prev(unitIDX).unit   = [];
prev(unitIDX).loaded = [];