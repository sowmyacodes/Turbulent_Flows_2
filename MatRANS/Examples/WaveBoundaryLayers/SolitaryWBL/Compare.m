% Compare computed bed shear stresses for solitary waves
% With theoretical (laminar, Liu et al. 2007) and measured 
% values (Sumer et al. 2010). 
%
%Liu, P.L.F. and Orfila, A. 2007 Boundary layer flow and bed shear stress under a solitary
%   wave. J. Fluid Mech. 574, 449-463.

%Sumer, B.B., Jensen, P.M., Sorensen, L.B., Fredsoe, J., Liu, P.L.F. and Carstensen, S. 
%   2010 Coherent structures in wave boundary layers. Part 2. Solitary motion. J. Fluid
%   Mech. 646, 207-231.
%
clear all; close all;

% Load output file, tests from Sumer et al. (2010)
load 'out_MatRANS.mat'; T = 7.9; % re=1.8e6, Sumer et al. (2010), test 15

% Calculate additional parameters
omega_w = 2*pi/T; % Angular frequency (1/s)
a = U_1m/omega_w; % Amplitude of orbital motion (m)
t0 = T; % Offset time (s)
rho = 1000; % Density (kg/m^3);

% Animate velocity profile and turbulent kinetic energy
figure(1), clf
for n = 1:length(t)
  subplot(1,2,1); plot(U(n,:),y); 
  xlim([min(min(U)) max(max(U))]); ylim([0 max(y)]);
  xlabel('u (m/s)'); ylabel('y (m)');
  title(['t = ' num2str(t(n)) ' s']);
  subplot(1,2,2); plot(K(n,:),y); xlim([0 max(max(K))]); ylim([0 max(y)]);
  xlabel('k (m^2/s^2)');
  pause(0.1);
  drawnow;
end

% Plot the free stream velocity time series
figure(2), clf
subplot(2,1,1)
omegat = 180/pi.*omega_w.*(t-t0);
plot(omegat,U(:,end));
ylabel('u_0 (m/s)');
xlim([-180 180]); set(gca,'xtick',[-360:30:360]);

% Plot the computed bed shear stress
subplot(2,1,2)
plot(omegat,tau_b./rho);
ylabel('\tau_0/\rho (m^2/s^2)');
xlabel('\omegat (^o)'); 
xlim([-180 180]); set(gca,'xtick',[-360:30:360]);
drawnow;

% Reynolds number, friction factor
Re_w = a*U_1m/nu
f_w = 2*max(tau_b)/rho/U_1m^2

% Add the laminar theoretical bed shear stress result (Liu et al. 2007)
disp('Computing laminar theory bed shear stress ...')
dt=1/10000; t = 0:dt:2; c=1;
u0 = sech(2*pi*(-c.*(t-1))).^2; % Free stream velocity
u0t = 4*pi.*sech(2*pi.*(1-t)).^2.*tanh(2*pi.*(1-t)); % Free stream acceleration% Perform integral
i = 0; int(1) = 0;
for tt = t % Peform numerical integral of their eq. (2.21)
  i = i + 1; % Update counter
  tau = t(1:i); 
  tmtau = tt - tau(1:end); % (t-tau)
  inttemp=cumtrapz(u0t(1:i)./sqrt(tmtau)).*dt; % Perform cumulative integral
  if i>1, int(i) = inttemp(end-1); end; % Store last value
end
tau0 = int./sqrt(pi); % Dimensionless bed shear stress
h = (3*U_1m*g*T^2/(sqrt(g)*16*pi^2))^(2/3); % Corresponding physical depth (m)
H = 16*pi^2*h^2/(3*g*T^2); % Corresponding physical solitary wave height (m)
tau0 = nu*(3*g^3*H^5/(16*h^4*pi^2*nu^2)).^(1/4).*tau0; % Dimensional tau_b/rho (m^2/s^2)
t = (t-1)*360; % Shift time to phase (in deg.)
inc = 200; hold on; plot(t(1:inc:end),tau0(1:inc:end),'b.','MarkerSize',8); hold off; % Add theoretical bed shear stress to plot

% Add measurements
load Data/Fig12f;
hold on; plot(x,y.*max(tau0).*1.3,'ro','MarkerSize',4); hold off;
leg=legend('Computed','Laminar theory','Sumer et al. (2010)','Location','NorthWest');
set(leg,'FontSize',7);