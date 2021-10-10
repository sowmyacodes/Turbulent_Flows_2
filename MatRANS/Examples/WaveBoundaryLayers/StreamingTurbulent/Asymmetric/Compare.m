% Compare laminar streaming with theory
clear all; close all;

% Load output file
load 'out_MatRANS.mat';

% Input
h_m = y(end); % Model height (m)
T = 6; % Wave period
nT = T/t(2); % Number of time steps per period
c = 7.5; % Wave celerity (m/s)

% Calculate theoretical friction velocity
omega_w = 2*pi/T; % Angular frequency (1/s)
a = U_1m/omega_w; % Orbital amplitude (m)
k = omega_w/c; % Wavenumber (1/m)

% Compare mean velocity profiles
figure(1), clf
for i = nT+1:10*nT:length(t)

  % Plot computed streaming profile
  ivec = i-nT+1:i; % Incices spanning one period
  Umean = mean(U(ivec,:)); % Average velocity over current period
  plot(Umean.*c./U_1m^2,y/a); % Computed
  xlabel('<u>c/U_{1m}^2'); ylabel('y/a');
  title(['t/T = ' num2str(t(i)/T)])
  ylim([0 h_m/a])
  drawnow
end

% Compute dimensionless streaming velocity
Us = Umean(end).*c./U_1m^2 % Computed

% Add line from Holmedal & Myrhaug
UsHM = 0.26.*3/4 % Holmedal & Myrhaug (2009)
hold on; plot(UsHM,y(end-10:end)/a,'ro'); hold off;

% Add legend
hleg = legend('Computed','Holmedal & Myrhaug (2009)');
set(hleg,'FontSize',8,'location','NorthWest');

% Compute normalized mean bed shear stress 
tau_b_mean = mean(tau_b(ivec))/rho/a/k/U_1m^2 
