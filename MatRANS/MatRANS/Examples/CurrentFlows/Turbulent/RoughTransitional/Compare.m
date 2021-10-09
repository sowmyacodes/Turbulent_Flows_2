% Compare rough-turbulent current results
clear all; close all;

% Load output file
load 'out_MatRANS.mat';

% Plot final computed velocity profile
Uf = 0.021; % Friction velocity (m/s)
y_kN = y./k_N(end); % Dimensionless vertical coordinate
subplot(1,2,1);
plot(U(length(U(:,1)),:)/Uf,y_kN); 
xlabel('u/U_f'); ylabel('y/k_N');
xlim([0 15]);
ylim([0 y_kN(end)]);
%set(gca,'yscale','log'); % Makes y-axis logarithmic

% Plot final computed k profile
subplot(1,2,2)
plot(K(length(K(:,1)),:)./Uf^2,y_kN);
xlabel('k/U_f^2'); ylabel('y/k_N');
xlim([0 3.5]);
ylim([0 y_kN(end)]);

% Compare with measurements
load 'Data.mat';
subplot(1,2,1);
hold on; plot(u_Uf,y_kN_u,'ro'); hold off; % Add velocity measurements
subplot(1,2,2);
hold on; plot(k_Uf2,y_kN_k,'ro'); hold off;

% Add logarithmic velocity profile
kappa = 0.4; % von Karman constant
u_log = 1/kappa*log(30*y_kN); % u/U_f
subplot(1,2,1)
hold on; plot(u_log,y_kN,'k--'); hold off;

% Legend
hl = legend('Computed','Fuhrman et al. (2010)','Logarithmic');
set(hl,'FontSize',6,'location','NorthWest');

% Reference
%
% Fuhrman, D.R., Dixen, M. & Jacobsen, N.G. (2010) Physically-consistent
%   wall boundary conditions for the k-omega turbulence model. 
%   J. Hydraul. Res. 48, 793-800.
