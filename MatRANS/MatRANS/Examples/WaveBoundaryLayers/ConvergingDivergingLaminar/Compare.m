clear all; close all;

% Load model output
load out_MatRANS;

% Input
nT = 36; % No. of saved time steps per period
T = 2*pi; omega = 2*pi/T; % Period (s); Angular frequency (1/s)
a = 1; % Amplitude of orbital bottom motion (m)
dadx = S*a/y(end); % da/dx = S*a/h

% Normalize vertical coordinates
delta1 = (nu*T/pi)^(0.5); % Stoke's length (m)
mu=y./delta1;

% Animate period-averaged velocity profile
np = 0; % Period counter
for i = 1:nT:length(t)-nT
   np = np + 1; % Update period counter   
   Um = mean(U(i:i+nT-1,:)); % Period-averaged horizontal velocity (m/s)
   plot(Um.*omega./(a*dadx),mu);
   xlim([-1 1].*1); ylim([0 mu(end)]);
   xlabel('<u>\omega/(a da/dx)'); 
   ylabel('y/\delta_1')
   title(['Period no. ' int2str(np)])
   drawnow;
end

% Add theoretical solution from Bijker et al. (1974)
uB = ( 0.5.*mu.*exp(-mu).*(sin(mu)-cos(mu)) + 2.*exp(-mu).*sin(mu) ...
+ 0.5.*exp(-mu).*cos(mu) + 0.25.*exp(-2.*mu) - 0.75 );
hold on; plot(uB,mu,'r--'); hold off;
title('');
hleg = legend('Computed','Bijker et al. (1974)'); set(hleg,'FontSize',8);
set(gcf,'position',[620 300 275 240])
xlim([-1 0.5]); ylim([0 10])

% Reference
%
% Bijker, E.W., Kalkwijk, J.P.T. & Pieters, T. (1974) Mass transport in
% gravity waves on a sloping bottom.  Proc. 14th International
% Conference on Coastal Engineering, p. 447-465.
