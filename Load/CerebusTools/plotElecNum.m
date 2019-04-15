function plotElecNum(fignum,tagname,A,d,elecs)
%% Plot electrode number function
% To use for figure 1, run plotElecNum(1,A,d,hObject), where (optionally)
% hObject contains the handles to all the subplots:
% hObject = findobj('-regexp','Tag','^elec[0-9]{1-2}plotNum');

% This requires a 96x1 vector d, which contains a d' value for each electrode
% (can remove this if you don't care about color of electrode number)

    if nargin < 5
       elecs = A';
       elecs = elecs(~isnan(elecs))';
    end
    
    dMax    = max(d);
    dMin    = min(d);
   for elec = elecs
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
            if ~isnan(elec)
                h = findobj('Tag',sprintf('elec%d%s',elec,tagname{fignum}));
%                 set(gcf,'CurrentAxes',h(1));
                subplot('Position',get(h,'Position'))

                dColor = [1 1-(d(elec)-dMin)/(dMax-dMin) 1-(d(elec)-dMin)/(dMax-dMin)];
                text('Tag',sprintf('elec%d%sNum',elec,tagname{fignum}),'Position',[0.5 0.5],...
                    'String',num2str(elec),'Units','normalized','Parent',h,...
                    'HorizontalAlignment','Center','VerticalAlignment','Middle',...
                    'FontWeight','Bold','FontUnits','normalized','BackgroundColor',...
                    'none','FontSize',1.25,'Color',dColor,'Visible','on');
            end
%         end
    end
end
% END Plot electrode number function