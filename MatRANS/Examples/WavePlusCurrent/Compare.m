% Compares wave plus current case
clear all; close all;

% Input
nT = 36; % Saved time steps per period
Ufc0 = 0.0157; % Current-only friction velocity (m/s)
Ufc = 0.0267; % Current friction velocity (m/s)
Ufw = 0.0662; % Wave friction velocity (m/s)
kw = 0.114; % Apparent wave roughness (m)
kappa = 0.4; % von Karman constant

% Animate pure current, used to generate initial conditions
load('CurrentOnly/out_MatRANS');
for n = 1:10:length(t)
  plot(U(n,:),y);
  xlabel('u (m/s)'); ylabel('y (m)');
  title(['Pure current: t = ' num2str(t(n)) ' s'])
  xlim([0 0.6]);
  drawnow;
end
uend = U(end,:); % Save final value

% Animate combined wave plus current
load('Combined/out_MatRANS');
for n = 1:nT:length(t)-36
  plot(mean(U(n+[0:35],:)),y);
  %hold on; plot(uend,y,'k--'); hold off;
  xlabel('<u> (m/s)'); ylabel('y (m)');
  title(['Wave plus current: t = ' num2str(t(n)) ' s'])
  xlim([0 0.6]);
  drawnow;  
end

% Add Fredsoe & Deigaard (1992), p. 61, analytical curves
k_N = k_N(end);
U320 = Ufc/kappa*Ufc/(Ufw+Ufc).*log(30.*y/k_N); % Inner solution
U313 = Ufc/kappa*log(30.*y/kw); % Outer solution
Ufd = max(U320,U313); % Take maximum of the two
hold on; plot(Ufd,y,'r--'); hold off;
set(gca,'yscale','log');

% Add logarithmic profile for pure current (assuming rough bed)
ulog = Ufc0./kappa.*log(30.*y./k_N);
hold on; plot(ulog,y,'k--'); hold off

% Add legend
legend('Wave plus current','Fredsoe & Deigaard (1992)','Current (log profile)', ...
    'location','SouthEast');
