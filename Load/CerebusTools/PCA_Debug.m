function f1h = PCA_Debug(f1h,xAxisInMsec,U,S)
   if ~isfield(f1h.fig,'PCA') || isempty(f1h.fig.PCA) || ~ishandle(f1h.fig.PCA)
      f1h.fig.PCA = figure('Name','sortnev2: PCAs',...
    'Toolbar','none','MenuBar','none','NumberTitle','off',...
    'Color',get(0,'defaultUicontrolBackgroundColor'));
   else
      figure(f1h.fig.PCA);
   end
   set(gcf,'DefaultAxesColorOrder',[[0.7 0   0]; [0 0.7 0  ]; [0 0   0.7]]);
   subplot(1,4,1:3)
   plot(xAxisInMsec,U);
   axis tight;
   title('First 3 PCAs')
   subplot(1,4,4)
   xr = bar(S);
   set(xr(1),'facecolor',[0.7 0   0  ])
   set(xr(2),'facecolor',[0   0.7 0  ])
   set(xr(3),'facecolor',[0   0   0.7])
   title('First 3 PCAs'' weights')
end