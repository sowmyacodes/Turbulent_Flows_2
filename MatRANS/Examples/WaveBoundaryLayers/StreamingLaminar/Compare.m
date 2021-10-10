% Compare laminar streaming with theory
clear all; close all;

% Load output file
load 'out_MatRANS.mat';

% Input
h_m = 1; % Model height (m)
U_1m = 1; % Free stream velocity magnitude (m/s)
T = 2*pi; % Wave period
nT = T/t(2); % Number of time steps per period
c = 20; % Wave celerity (m/s)

% Calculate theoretical friction velocity
omega_w = 2*pi/T; % Angular frequency
delta1 = sqrt(2.*nu/omega_w); % Stokes length
k = omega_w/c; % Wavenumber
xi = y./delta1; % Normalized vertical coordinates

% Compare mean velocity profiles
figure(1), clf
Usmax = 0;
for i = nT+1:10*nT:length(t) % Does so every 10th period

  % Plot computed streaming profile
  ivec = i-nT+1:i; % Indices spanning one period
  Umean = mean(U(ivec,:)); % Average velocity over current period
  if Umean(end) > Usmax % Check for max value
    Usmax = Umean(end);
    nt = i;
  end
  plot(Umean.*c./U_1m^2,xi); % Computed
  xlabel('<u>c/U_{1m}^2'); ylabel('y/\delta_1');
  title(['t/T = ' num2str(t(i)/T)])
  ylim([0 xi(end)])
  
  % Compare with laminar theory
  Ulam = 0.25.*(3 + exp(-xi).*(-4.*cos(xi) + 2.*sin(xi) + ...
			       exp(-xi) - 2.*xi.*sin(xi) - 2.*xi.*cos(xi)));
                               % See Nielsen, P. (1992, p. 56, with corrected typo)
  hold on; plot(Ulam,xi,'r--'); hold off;
  drawnow

end

% Add legend
hleg = legend('Computed','Theory');
set(hleg,'FontSize',8,'location','NorthWest');

% Compute dimensionless streaming velocity
Us = Umean(end).*c./U_1m^2 % Should be 3/4

% Compute normalized mean bed shear stress 
tau_b_mean = mean(tau_b(ivec))/rho/delta1/k/U_1m^2 % should be 1/4
