clear all; close all;

% Input
T = 5; % Period

% Load data
load out_MatRANS.mat;

% Mean velocity profile
nT = 36; 
n=12; ivec=n*nT:(n+1)*nT-1; % Use n-th period

% Mean suspended sediment flux profile
figure(1); clf;
ms=4; ytmp = y; load ../Data/Fig12a.mat; 
hold on; plot(x,y,'ro','MarkerSize',ms); hold off; y=ytmp;% Data
UC = C(ivec,:).*U(ivec,ib:end);
hold on; plot(mean(UC),yc); hold off;
xlim([-1 1].*0.03); ylim([-0.005 0.03]); 
xlabel('<{\ituc}> (m/s)');
ylabel('{\ity} (m)');
set(gca,'ytick',[0:0.01:y(end)])
box on;
set(gcf,'position',[879   292   267   225]);
UCint = trapz(yc,mean(UC))*1e6 % Integrate period-averaged suspended sediment flux

% Instantanous velocity profiles
omega = 2*pi/T; 
omegat = omega.*t(ivec).*180/pi; omegat = omegat - omegat(1); % + omega.*t0*180/pi;

% Plot instantaneous sediment flux profiles
% Add ensemble averaged measurements data
figure(2); ms = 4;
set(gcf,'position',[520   160   650   400])
subplot(2,4,1); load ../Data/Fig9a.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,2); load ../Data/Fig9b.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,3); load ../Data/Fig9c.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,4); load ../Data/Fig9d.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,5); load ../Data/Fig9e.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,6); load ../Data/Fig9f.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,7); load ../Data/Fig9g.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;
subplot(2,4,8); load ../Data/Fig9h.mat; hold on; plot(x,y,'ro','MarkerSize',ms); hold off;

% Computed results
phaseOW = [0 0.08 0.21 0.33 0.42 0.56 0.72 0.89].*360;
   % Phases from O'Donoghue & Wright (2004b), Fig. 10
figure(2);
for j = 1:length(phaseOW)
  uc = interp1(omegat,UC,phaseOW(j));
  subplot(2,4,j);
  hold on; plot(uc,yc); hold off;
  xlim([-1 1].*0.2); 
  ylim([-0.005 0.03]);
  set(gca,'ytick',[0:0.01:0.04]);
  if j > 4
    xlabel('{\ituc} (m/s)');
  end
  if j == 1 | j == 5
    ylabel('{\ity} (m)');
  end
  title(['{\it\omegat} = ' num2str(phaseOW(j)) '\circ']);
  box on;  
end

% Add legend
figure(2), legend('O''Donoghue & Wright (2004)','Model','Location','North','FontSize',6)
figure(1), legend('O''Donoghue & Wright (2004)','Model','Location','North','FontSize',6)
title('Period-averaged suspended sediment flux','FontSize',6)

% Calculate average sediment fluxes
qBm = mean(qB(ivec))*1e6 % Bed load
qSm = mean(qS(ivec))*1e6 % Suspended load
qTm = qBm + qSm % Total