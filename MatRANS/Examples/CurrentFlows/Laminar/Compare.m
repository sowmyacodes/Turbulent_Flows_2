% Compare laminar current result with theory
clear all; close all;

% Input
h_m = 1;
Px = -0.001;

% Load output file
load 'out_MatRANS.mat';

% Calculate theoretical friction velocity
Uf = sqrt(abs(Px)*h_m);
max(U_f)./Uf % Ratio should be close to 1

% Plot final computed velocity profile
plot(U(length(U(:,1)),:)/Uf,y./h_m); 
xlabel('u/U_f'); ylabel('y/h');

% Compare with laminar solution
Ulam = 1./nu.*(0.5.*Px*y.^2 - Px.*h_m.*y); % Analytical solution
hold on; plot(Ulam./Uf,y./h_m,'ro'); hold off;
legend('Computed','Theory','location','NorthWest')
