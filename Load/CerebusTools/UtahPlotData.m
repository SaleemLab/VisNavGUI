function fdata = UtahPlotData(animal, iseries, e, dataFieldName, enumFieldName, rotate, matchScale)
% fdata = UtahPlotData(animal, iseries, e, dataFieldName, enumFieldName) 
% 
% plots data according to the utah array layout
% plots electrode # according to utah array layout
%
% fdata = UtahPlotData(animal, iseries, e, dataFieldName, enumFieldName)
%   e is an array of structs, and has at least two fields containing
%       - electrode number
%       - data to plot
%   the names of the fields have to be passed as strings in dataFieldName and
%   enumFieldName
%
%   returns a handle to the data figure
%
% fdata = UtahPlotData(animal, iseries, e, dataFieldName, enumFieldName, rotate)
%   lets you rotate the array layout, options are 'ud' (up-down), 'lr'
%   (left-right) or a combination of both 'udlr' or 'lrud'
%
% fdata = UtahPlotData(animal, iseries, e, dataFieldName, enumFieldName, [], 'sameScale')
%   lets you specifty to plot all subplots with the same ylim (in case of
%   line plots) or clims (in case of images). Default: individual scales
%
% 2008-02-22 LB
% 2008-03-18 LB introduced rotate option
% 2008-08-25 LB corrected a bug concerning the figure with the enums that
%   was causing wrong numbering if one did not plot all the electrodes
% 2008-08-26 LB introduced the sameScale option

if nargin < 7
    matchScale = 'individualScale';
end
if nargin < 6
    rotate = [];
end

if ~isfield(e, dataFieldName) || ~isfield(e, enumFieldName)
    error('<UtahPlotData> Check format of e struct');
end

% get layout first
utahLayout = UtahGetLayout(animal, iseries);
if ~isempty(rotate)
    switch rotate
        case 'ud'
            utahLayout = flipud(utahLayout);
        case 'lr'
            utahLayout = fliplr(utahLayout);
        case {'udlr', 'lrup'}
            utahLayout = flipud(fliplr(utahLayout));
        otherwise
            fprintf('<UtahPlotData> Invalid string argument rotate %s', rotate);
    end
end

fdata = figure; hold on;
flabels = figure; hold on;

ax = nan(1, length(e));
for ie = 1 : length(e)
    if ie > (numel(utahLayout) - numel(find(isnan(utahLayout))))
        fprintf('Warning: data contains more electrodes than channels in the array\n');
        ax(ie) = [];
        continue;
    end
    [irow, icol] = find(utahLayout==e(ie).(enumFieldName));
    
    % figure containing the data
    figure(fdata);
    ax(ie) = gridplot(size(utahLayout,1), size(utahLayout,2), irow, icol);
    if find(size(e(1).(dataFieldName)) == 1) % array of values to plot
        plot(e(ie).(dataFieldName), 'k');
    else % images to plot
        imagesc(e(ie).(dataFieldName));
        axis off;
    end
    axis tight; axis square;
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    
    
    % figure containing the labels
    figure(flabels);
    gridplot(size(utahLayout,1), size(utahLayout,2), irow, icol);
    text(0.5,0.5, num2str(e(ie).(enumFieldName)));
    axis off;
    supertitle(sprintf('Utah array layout %s %d', animal, iseries), 0.98);
end

% set to the same scale
if any(size(e(1).(dataFieldName)) == 1) && strcmp(matchScale, 'sameScale')
    figure(fdata);
    matchy(ax);
end
if ~any(size(e(1).(dataFieldName)) == 1) && strcmp(matchScale, 'sameScale') && length(e) > 1
    figure(fdata);
    cl = get(ax, 'CLim');
    minC = min([cl{:}]);
    maxC = max([cl{:}]);
    set(ax, 'CLim', [minC maxC]);
end
