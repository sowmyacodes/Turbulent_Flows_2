% Compare with the Vanoni distribution for suspended sediment 
clear all; close all;

% Load output file
load 'out_MatRANS.mat';
k_N = k_N(end);

% Plot final computed velocity profile
Uf = 0.1; % Friction velocity (m/s)
%Uf = U_f(end); % Friction velocity (m/s)
h_m = y(end);
subplot(1,3,1);
plot(U(length(U(:,1)),:)/Uf,y./h_m); 
xlabel('U/U_f'); ylabel('y/h');
xlim([0 30]);
ylim([0 y(end)/h_m]);
%set(gca,'yscale','log'); % Makes y-axis logarithmic

% Add logarithmic velocity profile
kappa = 0.4; % von Karman constant
u_log = 1/kappa*log(30*y/k_N); % u/U_f
subplot(1,3,1)
hold on; plot(u_log,y/h_m,'k--'); hold off;


% Plot eddy viscosity
subplot(1,3,2);
plot(nu_T(end,:)./Uf/h_m,y./h_m);
xlabel('\nu_T/(U_f h)'); ylabel('y/h');
xlim([0 0.15])
%set(gca,'yscale','log'); % Makes y-axis logarithmic

% Add analytical result
kappa = 0.4;
nuT = kappa.*y./h_m.*(1 - y./h_m);
hold on; plot(nuT,y./h_m,'k--'); hold off;


% Plot suspended sediment distribution
subplot(1,3,3);
c_b = C(end,1);
plot(C(end,:)./c_b,yc./h_m);
set(gca,'yscale','log')
xlabel('c/c_b')
ylabel('y/h');

% Add Vanoni distribution
Rouse = w_s./kappa./Uf; % Rouse parameter
b = yc(1); % Reference level
cV = ( (h_m-yc)./yc.*b./(h_m-b) ).^Rouse;
hold on; plot(cV,yc./h_m,'k--'); hold off;


% Legend
hl = legend('Computed','Vanoni');
set(hl,'FontSize',6);
