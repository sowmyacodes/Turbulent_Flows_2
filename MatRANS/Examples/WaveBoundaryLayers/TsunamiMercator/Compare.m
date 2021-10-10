% Compare computed bed shear stresses for solitary waves
clear all; close all;

% Load output file, tests from Sumer et al. (2010)
load 'out_MatRANS.mat';

% Additional parameters
rho = 1000; % Density (kg/m^3);

% Animate velocity profile and turbulent kinetic energy
figure(1), clf
for n = 1:length(t)

  subplot(1,2,1); plot(U(n,:),y); 
  xlim([min(min(U)) max(max(U))]); ylim([0 max(y)]);
  xlabel('u (m/s)'); ylabel('y (m)');
  title(['t = ' num2str(t(n),4) ' s']);
  subplot(1,2,2); plot(K(n,:),y); 
  xlim([0 max(max(K))]); ylim([0 max(y)]);
  xlabel('k (m^2/s^2)');
  %pause(0.1);
  drawnow;

end

% Plot the computed free stream velocity time series
figure(2), clf
subplot(3,1,1)
plot(t,U(:,end));
ylabel('u_0 (m/s)');
xlim([0 t(end)]); 

% Add the digitized velocity
load Data/MercatorVelocity.mat;
hold on; plot(t_ud,U_ud,'r--'); hold off;

% Add legend
hleg = legend('Computed','Target','Location','NorthWest');
set(hleg,'FontSize',8);

% Plot the computed bed shear stress
subplot(3,1,2)
plot(t,tau_b./rho);
ylabel('\tau_0/\rho (m^2/s^2)');
xlim([0 t(end)]); 
%xlabel('t (s)'); 

% Plot the computed Shields parameter
subplot(3,1,3)
plot(t,Shields);
ylabel('\theta');
xlabel('t (s)'); 
xlim([0 t(end)]); 
drawnow;







