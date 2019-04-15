function plotElecD(fignum,tagname,A,d,elecs)
% Plot electrode d' function
    if nargin < 5
       elecs = A';
       elecs = elecs(~isnan(elecs))';
    end
    
   for elec = elecs
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
            if ~isnan(elec)
                h = findobj('Tag',['elec',num2str(elec),tagname{fignum}]);
                set(gcf,'CurrentAxes',h(1));
                
                text('Tag',['elec',num2str(elec),tagname{fignum},'D'],'Position',...
                    [0.99 0.03],'String',num2str(d(elec),'%2.2g'),'Units','normalized',...
                    'HorizontalAlignment','Right','VerticalAlignment','Bottom',...
                    'FontWeight','Bold','FontUnits','normalized','BackgroundColor',...
                    [1 0.7 0.7],'FontSize',0.15,'Color',[1 0 0],'Visible','on');%,...
%                     'TooltipString',sprintf('Elec: %d',elec));
            end
%         end
    end
end
% END Plot electrode d' function