function [Xrealigned,iXraw] = changeThreshold(threshold,Xraw,f1h,nevopen_outcome)
% CHANGETHRESHOLD: Use Dynamic Multiphasic Filter method of thresholding
% (from Blanche et al. Manuscript)
%
% Outputs:
% Xrealigned is Xraw(:,~iXraw) combined with realigned Xraw(:,iXraw)
% Xrealigned(:,iXraw) are all the traces (realigned) above threshold
% Xraw(:,iXraw) are all the raw traces (not realigned) above threshold
% 
% 2009-01 AZ Created
% 2010-03 AZ Added alignToSampleNum variable, changed it from 12th to 10th sample
if nargin < 4
   nevopen_outcome = true;
end

alignToSampleNum = 10;

if isnan(threshold)
   if ~isempty(f1h)
      angry;
   end
   
   threshold = 32;
   set(f1h.txtbox.thresh,'String',num2str(threshold));
end

   if nevopen_outcome
      if ~isempty(f1h), set(gcf,'Pointer','watch'); drawnow; end
      iXraw = false(size(Xraw,2),1);
      if threshold
%          % BIPOLAR AMPLITUDE THRESHOLD Method
%          iXraw = zeros(size(Xraw,2),1);
%          % Index all waveforms with spikes crossing +/- threshold
%          % [s_i(t) = Xraw, f = threshold]
%          % Check if abs(s_i(t)) > f
%          iXraw( intersect( find(any(Xraw >  abs(threshold),1)) ,...
%                            find(any(Xraw < -abs(threshold),1)) ) ) = 1;
%          % End BIPOLAR AMPLITUDE THRESHOLD Method
         
         % DYNAMIC MULTIPHASIC FILTER Method
         deltat = 8;
         for i = 1:size(Xraw,2)
         % Noise Estimation
%             threshold(i) = 3.6*median(abs(Xraw(:,i)));
         % END Noise Estimation
         % If noise is always 10, dynamic thresh 'f' is always in 31-62 range
         % Uncomment the following line if not interested in noise estimation
%             threshold(i) = threshold(1);
            xr = find(Xraw(:,i) < -abs(threshold));
            yr = find(Xraw(:,i) >  abs(threshold));
            if ~isempty(xr) || ~isempty(yr)
               [a,j] = min([min(xr) min(yr)]);
               % [s_i(t) = Xraw, f = threshold, deltat = deltat]
               % -deltat < t' <= deltat
               if ~isempty(xr) && j == 1
                  % CASE WHERE s_i(t) < -f
                  [a,j] = min(Xraw(xr,i));
                  a = max([xr(j)-deltat 1]):min([xr(j)+deltat 48]);
                  % [f_val = Xraw(xr(j),i), i.e., the valley of the waveform]
                  % Check if s_i(t + t') > f_val + 2*f
                  iXraw(i) = ...
                     any(Xraw(a,i) > Xraw(xr(j),i)+2*abs(threshold));
               elseif j == 2 || (~isempty(yr) && j == 1)
                  % CASE WHERE s_i(t) >  f
                  [a,j] = max(Xraw(yr,i));
                  a = max([yr(j)-deltat 1]):min([yr(j)+deltat 48]);
                  % [f_pk  = Xraw(yr(j),i), i.e., the peak   of the waveform]
                  % Check if s_i(t + t') < f_pk  - 2*f
                  iXraw(i) = ...
                     any(Xraw(a,i) < Xraw(yr(j),i)-2*abs(threshold));
               end
%             else % if xr and yr are both empty
%                iXraw(i) = 0;
            end
         end
         % Make sure iXraw is proper length
         if length(iXraw) < i
            iXraw(i) = false;
         end
            
         % End DYNAMIC MULTIPHASIC FILTER Method

         Xrealigned = Xraw(:,iXraw);
      else
         Xrealigned = Xraw;
      end
      
      %% REALIGN SUPRA-THRESHOLD TRACES
      if size(Xrealigned,2) > 1
         if ~isempty(find(iXraw, 1))
            % %% AZ20090310 DEBUG %%
            % xr = []; % xr = Xrealigned, horizontally concatenated
            % yr = []; % yr = threshold for i'th window of Xrealigned, concatenated
            % for i = 1:size(Xrealigned,2)
            %    xr = [xr Xrealigned(:,i)'];
            %    yr = [yr threshold(i)'*ones(1,size(Xrealigned,1))];
            % end
            % 
            % figure(6); hold on
            % % fill([1:1:size(xr,2),size(xr,2):-1:1],[yr,fliplr(-yr)],[1 0.8 0.8],'EdgeColor','none')
            % plot( xr,'-r')
            % plot( yr,'-g')
            % plot(-yr,'-g')
            % 
            % figure(2)
            % plot(Xrealigned)
            % figure(3)
            % plot(Xraw(:,find(~iXraw)))
            % %% END AZ20090310 DEBUG %%

            % Re-align offset traces above threshold. All global minima occur on
            % sample 12 (t = 0.3667 ms).
            % [xr,yr] = temporary variables (values & indices of minima)
            [xr,yr] = min(Xrealigned,[],1);
   % % AZ 20090311 Debug
   % figure(4)
   % hist(yr,48)
   % % End AZ 20090311 Debug

            % Reincorporate traces that didn't pass threshold
            % Re-expand Xrealigned to size of Xraw
            Xrealigned(:,iXraw) = Xrealigned;
            % Add missing rows to end
            Xrealigned = [Xrealigned zeros(48,size(Xraw,2) - size(Xrealigned,2))];
            % Add in traces that didn't pass threshold
            Xrealigned(:,~iXraw) = Xraw(:,~iXraw);

            % Re-expand yr to size of Xraw
            yr(:,iXraw) = yr;
            % REALIGN traces that did pass threshold
            for toRealign = intersect(find(yr'-alignToSampleNum ~= 0),find(iXraw))'
               % Must accomodate left- and right-ward shifts
               Xrealigned(:,toRealign) = ...
                  [ zeros(1,(max(      -(yr(toRealign)-alignToSampleNum ),0) ))      ... % Pad left , if needed
                  Xrealigned(max(        yr(toRealign)-alignToSampleNum+1,1        ):...
                             min(end,end+yr(toRealign)-alignToSampleNum ),toRealign)'... % Snippet
                    zeros(1, max(        yr(toRealign)-alignToSampleNum  ,0  ))        ];% Pad right, if needed
            end
         end
      else
         Xrealigned = [];
         if ~isempty(f1h)
            set(findobj(gcf,'Tag','txt.status'),'String','Threshold too high!');
            set(gcf,'Pointer','arrow');
         end
      end
   else
      disp('No .nev file loaded.')
   end % end if .nev file loaded
end