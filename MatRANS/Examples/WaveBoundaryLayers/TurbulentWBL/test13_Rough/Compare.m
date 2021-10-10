% Compare laminar wave boundary layer result with theory
clear all; close all;

% Load output file
load 'out_MatRANS.mat';

% Input
T = 9.72; % Period (s)
omega_w = 2*pi/T; % Angular frequency (1/s)
a = U_1m/omega_w;
dt = t(2) - t(1);
N = T/dt; % Number of saved times per period


% Compare velocity profiles for acceleration phases (15 degree increments)
inc = N*15/360; % Degree increments
figure(1), clf
set(gcf,'Position',[300 300 600 270]);
subplot(1,2,1);
i0 = length(t)-72;
for i = i0:inc:i0+6*inc
  plot(U(i,:)/U_1m,y./a); % Computed
  xlim([-0.5 1.201]); ylim([0 y(end)/a])
  %title(['t/T = ' num2str(t(i)/T)])
  drawnow
  hold on;
end
set(gca,'ytick',[0:0.01:1])
xlabel('u/U_{1m}'); ylabel('y/a');
title('Acceleration phases, \omegat=0-90 deg.')

% Add measurements
Umeas = load('../Data/P13MEANU.DAT');
ymeas = Umeas(:,1)./1000; % y (m)
Umeas = Umeas(:,2:end)/1000;
hold on;
for i = 1:7;%length(Umeas(1,:))
  plot(Umeas(:,i)/U_1m,ymeas/a,'ro');
  drawnow;
end

% Also for deceleration phases
subplot(1,2,2);
for i = i0+6*inc:inc:i0+12*inc
  plot(U(i,:)/U_1m,y./a); % Computed
  xlim([-0.5 1.201]); ylim([0 y(end)/a])
  %title(['t/T = ' num2str(t(i)/T)])
  drawnow
  hold on;
end
set(gca,'ytick',[0:0.01:1])
xlabel('u/U_{1m}'); ylabel('y/a');
title('Deceleration phases, \omegat = 90-180 deg.')

% Add measurements
hold on;
for i = 7:13;
  plot(Umeas(:,i)/U_1m,ymeas/a,'ro');
  drawnow;
end


% Compare turbulence quantities
uu = load('../Data/P13RMSUU.DAT');
yuu = uu(:,1)./1000; % y (m)
uu = uu(:,2:end)./1000; % RMS uu (m/s)
vv = load('../Data/P13RMSVV.DAT');
yvv = vv(:,1)./1000; % y (m)
vv = vv(:,2:end)./1000; % RMS vv (m/s)
km = 0.65.*(uu.^2 + vv.^2); % Approximate of k from Justesen (1991), (m^2/s^2)

figure(2), clf
set(gcf,'Position',[300 100 800 500]);
count = 0;
for i = i0:inc:i0+12*inc
  count = count + 1;
  subplot(2,7,count);
  plot(K(i,:)/U_1m^2,y./a); % Computed
  xlim([0 0.01]); ylim([0 y(end)/a])
  hold on; plot(km(:,count)./U_1m^2,yuu/a,'ro'); hold off;
  set(gca,'xtick',[0 0.01],'xticklabel',[0 0.01])
  
  % Add titles to subplots
  phase = (count-1).*15; 
  if i == i0;
    title(['\omegat=' num2str(phase) ' deg.'])
  else
    title([num2str(phase) ' deg.'])
  end
  if count > 7
    xlabel('k/U_{1m}^2')
  end
  drawnow
end
subplot(2,7,1); ylabel('y/a');
subplot(2,7,8); ylabel('y/a');

% Reynolds stress comparison
uv = load('../Data/P13REYUV.DAT');
yuv = uv(:,1)./1000; % y (m)
uv = uv(:,2:end)./1000^2; % uu (m^2/s^2)
figure(3), clf
set(gcf,'Position',[300 100 800 500]);
count = 0;
for i = i0:inc:i0+12*inc
  count = count + 1;
  subplot(2,7,count);
  uv_model = -nu_T(i,:).*gradient(U(i,:),y); % -nu_T*du/dy
  plot(uv_model./U_1m^2,y./a); % Computed Reynolds stresses
  xlim(0.005.*[-1 1]); ylim([0 y(end)/a])
  hold on; plot(uv(:,count)./U_1m^2,yuv/a,'ro'); hold off;
  set(gca,'xtick',[-1 0 1].*0.005,'xticklabel',[-1 0 1].*0.005)
  
  % Add titles to subplots
  phase = (count-1).*15; 
  if i == i0;
    title(['\omegat=' num2str(phase) ' deg.'])
  else
    title([num2str(phase) ' deg.'])
  end
  if count > 7
    xlabel('$\overline{u^\prime v^\prime}/U_{1m}^2$','Interpreter','Latex')
  end    
  drawnow
end
subplot(2,7,1); ylabel('y/a');
subplot(2,7,8); ylabel('y/a');


% Bed shear stress comparison
figure(4), clf
plot(t/T*360-4*360,sign(tau_b).*sqrt(abs(tau_b)/rho)/U_1m);
xlabel('\omegat (deg)'); ylabel('U_f/U_{1m}');
xlim([0 360])
set(gca,'xtick',[0:45:360])
pos = get(gcf,'position');
set(gcf,'position',[1 1 1 0.5].*pos);

% Add measurements
load '../Data/Uf_test13.mat';
hold on; plot(phasem,Ufm/U_1m/100,'ro'); hold off;
