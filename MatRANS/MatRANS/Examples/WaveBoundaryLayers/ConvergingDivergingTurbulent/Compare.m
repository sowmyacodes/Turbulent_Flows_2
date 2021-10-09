% Compare against Test 16 of Sumer et al. (1993)
clear all; close all;

% Input
nT = 36; % Saved time levels per period

% Load data
load out_MatRANS;

% Animate period averaged velocity
ip = 0;
for i = 1:nT:length(t)-nT+1
  ip = ip + 1; % Count no. of periods
  ivec = i + [0:35];
  ubar = mean(U(ivec,:)); % Period-averaged velocity profile
  plot(ubar/U_1m,y/y(end),'b-');
  xlabel('<{\itu}>/{\itU}_{1m}'); ylabel('{\ity}/{\ith}');
  uavg = trapz(y,ubar)./y(end); % Depth-averaged velocity
  title(['t/T = ' int2str(ip) ', q/h/U_{1m} = ' num2str(uavg/U_1m)]);
  %xlim([-1 1].*0.06); ylim([0 2])  
  xlim([-1 1].*0.3); ylim([0 2])  
  pause(0.1);
  drawnow;
end

% Add data
load Data/Data16.mat; % Digitized from Sumer et al. (1993), Fig. 7b
hold on; plot(xdata/100/U_1m,ydata/100/y(end),'ro','MarkerSize',4); hold off; % Experiments

% Adjust figure
set(gcf,'position',[200 200 300 300]) % Change size
title(''); % Remove the title
xlim([-1 1].*0.05);
set(gca,'xtick',[-0.05:0.025:0.05],'ytick',[0:0.5:2]);
grid on;

% Reflect computed data across vertical centerline
hold on; plot(ubar(end:-1:1)/U_1m,2-y(end:-1:1)./y(end),'b-'); hold off;

% Add legend
hl = legend('Computed','Sumer et al. (1993)');
set(hl,'FontSize',6,'location','SouthEast');