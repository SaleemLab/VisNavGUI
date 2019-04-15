% RETHRESH.M

% Get Threshold
if ~isempty(get(f1h.txtbox.thresh,'String'))
   threshold = str2double(get(f1h.txtbox.thresh,'String'));
else
   threshold = [];
end
[Xrealigned,iXraw] = changeThreshold(threshold,Xraw,f1h,nevopen_outcome);

% RESIZE ( take only samples (7:34) (t = -0.1:1/30:0.8) )
if size(Xrealigned ,1) == nSamplesPerPacket
   Xrealigned = Xrealigned(x(7:end-14),:);
end
if size(xAxisInMsec,2) == nSamplesPerPacket
   xAxisInMsec = (x(7:end-14)-10)/SamplingRateInKHZ; % center about thresholding point (10)
end

% RECALCULATE PCA
[U,S,V] = svds(Xrealigned,3);