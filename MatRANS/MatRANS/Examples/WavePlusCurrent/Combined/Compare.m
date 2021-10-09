% Compare laminar streaming with theory
clear all; close all;

% Load output file
load 'out_MatRANS.mat';

% Input
h_m = y(end); % Model height (m)
%U_1m = 1; % Free stream velocity magnitude (m/s)
T = 6; % Wave period
nT = T/t(2); % Number of time steps per period
c = 7.5; % Wave celerity (m/s)

% Calculate theoretical friction velocity
omega_w = 2*pi/T; % Angular frequency (1/s)
a = U_1m/omega_w; % Orbital amplitude (m)
%delta1 = sqrt(2.*nu/omega_w); % Stokes length
k = omega_w/c; % Wavenumber (1/m)
%xi = y./a; % Normalized vertical coordinates

% Compare mean velocity profiles
figure(1), clf
for i = nT+1:10*nT:length(t)

  % Plot computed streaming profile
  ivec = i-nT+1:i; % Incices spanning one period
  Umean = mean(U(ivec,:)); % Average velocity over current period
  plot(Umean.*c./U_1m^2,y/a); % Computed
  %plot(Umean,y);
  xlabel('U_sc/U_{1m}^2'); ylabel('y/a');
  title(['t/T = ' num2str(t(i)/T)])
  ylim([0 h_m/a])
  
  % Compare with laminar theory
%  Ulam = 0.25.*(3 + exp(-xi).*(-4.*cos(xi) + 2.*sin(xi) + ...
%			       exp(-xi) - 2.*xi.*sin(xi) - 2.*xi.*cos(xi)));
                               % See Nielsen, P. (1992, p. 56, with corrected typo)
  %hold on; plot(Ulam,xi,'r--'); hold off;
  drawnow

end

% Compute dimensionless streaming velocity
Us = Umean(end).*c./U_1m^2 % 

% Add line from Holmedal & Myrhaug
UsHM = 0.46.*3/4 % Holmedal & Myrhaug (20??)
hold on; plot(UsHM,y/a,'r--'); hold off;

% Compute normalized mean bed shear stress 
tau_b_mean = mean(tau_b(ivec))/rho/a/k/U_1m^2 
