% Compare laminar wave boundary layer result with theory
clear all; close all;

% Input
h_m = 1; % Model height (m)
U_1m = 1; % Free stream velocity magnitude (m/s)
T = 2*pi; % Wave period

% Load output file
load 'out_MatRANS.mat';

% Calculate theoretical friction velocity
omega_w = 2*pi/T; % Angular frequency
delta1 = sqrt(2.*nu/omega_w); % Stokes length

% Compare animated velocity profiles
figure(1), clf
for i = 1:length(t)
  plot(U(i,:)/U_1m,y./delta1); % Computed
  Ulam = sin(omega_w.*t(i)) - exp(-y./delta1).*sin(omega_w.*t(i) - y./delta1); % Laminar theory
  hold on; plot(Ulam,y./delta1,'r--'); hold off;
  xlabel('u/U_{1m}'); ylabel('y/\delta_1');
  xlim(1.2.*[-1 1]); ylim([0 h_m./delta1])
  title(['t/T = ' num2str(t(i)/T)])
  drawnow
  hold on;
end
legend('Computed','Theory');

% Compare bed shear stress time series
figure(2), clf
plot(omega_w.*t,tau_b/rho*delta1/U_1m/nu);
tau_lam = (sin(omega_w.*t)+cos(omega_w.*t)); % tau_b*delta1/(rho*nu*U_1m) from laminar theory
hold on; plot(omega_w.*t,tau_lam,'ro'); hold off;
xlabel('\omega t'); ylabel('\tau_b\delta_1/(\rho\nu U_{1m})')
legend('Computed','Theory');
