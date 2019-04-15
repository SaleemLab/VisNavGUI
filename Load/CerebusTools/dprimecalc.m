function d = dprimecalc(V,signal_idx,noise_idx)
% DPRIMECALC.M: Run actual calulations for d' measurement
%               Plot diagnostics are commented out, can re-enable

Vsignal_mean = mean(V(signal_idx,:),1);
 Vnoise_mean = mean(V(noise_idx ,:),1);


% col    = {[1 0 0] [0 0 1] [0 1 0] [1 0 1] [0.5 0.5 0.5]};
% figure(5); hold on
% plot3(V(signal_idx,1),V(signal_idx,2),V(signal_idx,3),'Color',col{1},'LineStyle','.','markersize',4)
% plot3(V(noise_idx ,1),V(noise_idx ,2),V(noise_idx ,3),'Color',col{2},'LineStyle','.','markersize',4)
% plot3([Vsignal_mean(1); Vnoise_mean(1)],[Vsignal_mean(2); Vnoise_mean(2)],[Vsignal_mean(3); Vnoise_mean(3)], 'k-')

a = mean([Vsignal_mean; Vnoise_mean]);

Vn(signal_idx,:) = V(signal_idx,:) - ones(size(signal_idx))*a;
Vn(noise_idx ,:) = V(noise_idx ,:) - ones(size(noise_idx ))*a;


% figure(6); hold on
% plot3(Vn(signal_idx,1),Vn(signal_idx,2),Vn(signal_idx,3),'Color',col{1},'LineStyle','.','markersize',4)
% plot3(Vn(noise_idx ,1),Vn(noise_idx ,2),Vn(noise_idx ,3),'Color',col{2},'LineStyle','.','markersize',4)
% plot3([Vsignal_mean(1); Vnoise_mean(1)]-a(1),[Vsignal_mean(2); Vnoise_mean(2)]-a(2),[Vsignal_mean(3); Vnoise_mean(3)]-a(3), 'k-')

b = Vsignal_mean - a;

thetay = atan(b(1)/b(3));
RY = [ cos(thetay), 0 , sin(thetay);
                 0, 1 , 0;
      -sin(thetay), 0 , cos(thetay)];

Vsig = Vn(signal_idx,:);
Vnoi = Vn(noise_idx ,:);
Vsig = Vsig*RY;
Vnoi = Vnoi*RY;

VsigM = mean(Vsig,1);
VnoiM = mean(Vnoi,1);

thetaz = atan(VsigM(2)/VsigM(1));
RZ = [ cos(thetaz), -sin(thetaz), 0;
       sin(thetaz),  cos(thetaz), 0;
                 0,            0, 1];
				 
Vsig = Vsig*RZ;
Vnoi = Vnoi*RZ;

VsigM = mean(Vsig,1);
VnoiM = mean(Vnoi,1);
VsigD = std(Vsig(:,1));
VnoiD = std(Vnoi(:,1));
d = abs(VsigM(1) - VnoiM(1))/sqrt((VsigD^2 + VnoiD^2)/2);

% figure(7); hold on
% plot3([VsigM(1) VnoiM(1)],[VsigM(2) VnoiM(2)],[VsigM(3) VnoiM(3)], 'k-')
% plot3(Vsig(:,1),Vsig(:,2),Vsig(:,3),'Color',col{1},'LineStyle','.','markersize',4)
% plot3(Vnoi(:,1),Vnoi(:,2),Vnoi(:,3),'Color',col{2},'LineStyle','.','markersize',4)

if length(Vsig) < length(Vnoi)
   Vsig = [Vsig; nan(size(Vnoi,1)-size(Vsig,1),3)];
elseif length(Vsig) > length(Vnoi)
   Vnoi = [Vnoi; nan(size(Vsig,1)-size(Vnoi,1),3)];
end
% x = -0.06:0.001:0.06;
% figure(8); hist([Vsig(:,1) Vnoi(:,1)],x)
% axis([-0.06 0.04 0 160])

end