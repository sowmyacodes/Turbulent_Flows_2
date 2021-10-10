% Compare smooth-turbulent current results
clear all; close all;

% Load output file
load 'out_MatRANS.mat';

% Plot final computed velocity profile
Uf = 0.016; U_f(end); % Final friction velocity
yp = y.*Uf/nu; % y^+ = y*U_f/nu
subplot(1,2,1);
plot(U(length(U(:,1)),:)/Uf,yp); 
xlabel('u/U_f'); ylabel('y^+');
ylim([0 yp(end)]);
%set(gca,'yscale','log'); % Makes y-axis logarithmic

% Plot final computed k profile
subplot(1,2,2)
plot(K(length(K(:,1)),:)./Uf^2,yp);
xlabel('k/U_f^2'); ylabel('y^+');
ylim([0 yp(end)]);

% Compare with measurements
load 'DataSmooth.mat';
subplot(1,2,1);
hold on; plot(u_Uf,yp_u,'ro'); hold off; % Add velocity measurements
subplot(1,2,2);
hold on; plot(k_Uf2,yp_k,'ro'); hold off;

% Add logarithmic velocity profile
kappa = 0.4; % von Karman constant
u_log = 1/kappa*log(yp) + 5.1; % u/U_f
subplot(1,2,1)
hold on; plot(u_log,yp,'k--'); hold off;

% Legend
hl = legend('Computed','Fuhrman et al. (2010)','Logarithmic');
set(hl,'FontSize',6,'location','NorthWest');

% Reference
%
% Fuhrman, D.R., Dixen, M. & Jacobsen, N.G. (2010) Physically-consistent
%   wall boundary conditions for the k-omega turbulence model. 
%   J. Hydraul. Res. 48, 793-800.


