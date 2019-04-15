function calculated = elecSubplotCreate(fignum,tagname,A,calculated,elecs)
% ELECSUBPLOTCREATE.M: Create Subplot for electrode function
    if nargin < 5
       elecs = A';
       elecs = elecs(~isnan(elecs))';
    end
    
   num_y = size(A,1);
	num_x = size(A,2);
   
   % Select the right figure window
   set(0,'CurrentFigure',findobj(get(0,'Children'),'Tag','expnev1.fig'));
            
   for elec = elecs
%     for i = 1:10
%         for j = 1:10
%             elec = A(i,j); %axes(findobj('Tag','Special scaling_params'))
            [y,x] = ind2sub(size(A),find(A==elec));
            if ~isnan(elec)
                subplot('Position',...%[Xoffset+(j-1)*xoffset 1-Yoffset-(i-1)*yoffset xsize ysize],... %                    [(x-1)*(0.8)/num_x-0.8/num_x/20 1-y*(1.0)/num_y-1.0/num_y/20 ... %                           (0.8)/num_x+0.8/num_x/10     (1.0)/num_y+1.0/num_y/10],...
                     [(x-1)*(0.8)/num_x 1-y*(1.0)/num_y  ...
                            (0.8)/num_x     (1.0)/num_y],...
                    'Tag',sprintf('elec%d%s',elec,tagname{fignum}),'DrawMode','fast',...
                    'XTick',[],'YTick',[],'XTickLabel',{},'YTickLabel',{},'Visible','off');
                hold on;
                
%                 set(gca,'OuterPosition',get(gca,'Position'));
                set(gca, 'Color', get(gcf, 'Color'));
%                 if     y == 10% && fignum == 1
%                     set(gca,'XTick',[],'XTickLabel',{});
%                 elseif x ==  1% && fignum == 1
%                     set(gca,'YTick',[],'YTickLabel',{});
%                 else
%                     set(gca,'XTick',[],'YTick',[],'XTickLabel',{},'YTickLabel',{});
%                 end
            end
%         end
    end
    
    calculated(fignum).D      = 0;
    calculated(fignum).Num    = 0;
    calculated(fignum).Traces = 0;
    calculated(fignum).Axes   = 0;  % make sure this is always 0
end
% END Create Subplot for electrode function