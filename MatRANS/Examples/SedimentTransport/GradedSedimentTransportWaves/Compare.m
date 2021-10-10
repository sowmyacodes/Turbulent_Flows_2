% Compute sediment fluxes, compare transport rates
clear all; close all;

% Input
nT = 12; % No. of period
T = 6.5; % Period (s)
marker = 'k-'; % Used in plots

% Load the data
load out_MatRANS;

% Plot Shields parameters
figure(100); set(gcf,'units','centimeter','position',[3 15 9 4]);
ivec = nT*36 + [1:36]; ivecp = [ivec ivec(end)+1];
omegat = t(ivecp).*2*pi/T; omegat = omegat - nT*2*pi; omegat = omegat/pi*180;
count = 0;
for nn = find(w_f)
    count = count + 1; hold on;
    if count == 1
        plot(omegat,Shields(ivecp,nn),'k-'); 
    else
        plot(omegat,Shields(ivecp,nn),'r--'); 
    end; hold off; 
end
xlim([0 360]); set(gca,'xtick',[0:90:360])
xlabel('\omegat (^o)','FontSize',8); ylabel('\theta_i','FontSize',8)
hleg1=legend('d_i=0.13 mm','d_i=0.34 mm'); set(hleg1,'FontSize',8)
box on;

% Plot suspended load
figure(101); set(gcf,'units','centimeter','position',[3 3 9 9]);
nplots = sum(w_f>0); count = 0;
for nn = find(w_f)
    count = count + 1;
   subplot(nplots,1,count);
   hold on; plot(omegat,qS(ivecp,nn)*1e6,marker); hold off;
   xlim([0 360]); set(gca,'xtick',[0:90:360])
   ylabel('q_{S,i}\times 10^6 (m^2/s)','FontSize',8); box on;
   yl = get(gca,'ylim'); 
   if count == 1, let='a'; else, let='b'; end
   title(['(' let ') d_i=' num2str(d(nn)*1e3) ' mm'],'FontSize',8)
   box on;
end
xlabel('\omegat (^o)','FontSize',8);  

% Plot bedload
figure(102); set(gcf,'units','centimeter','position',[14 3 9 9]);
nplots = sum(w_f>0); count = 0;
for nn = find(w_f)
   count = count + 1;
   subplot(nplots,1,count);
   hold on; plot(omegat,qB(ivecp,nn)*1e6,marker); hold off;
   xlim([0 360]); set(gca,'xtick',[0:90:360])
   ylabel('q_{B,i}\times 10^6 (m^2/s)','FontSize',8); box on;
   yl = get(gca,'ylim'); 
   if count == 1, let='a'; else, let='b'; end
   title(['(' let ') d_i=' num2str(d(nn)*1e3) ' mm'],'FontSize',8)
   box on;
end
xlabel('\omegat (^o)','FontSize',8); 

% Plot concentration profiles
imax = find(U(:,end)==max(U(ivec,end))); % Index where flow is max
figure(104); set(gcf,'units','centimeter','position',[26 3 9 7]);
count = 0;
for nn = find(w_f)
    count = count + 1;
    subplot(1,nplots,count);
    jj = (nn-1).*length(yc) + [1:length(yc)];
    hold on; plot(C(imax,jj),yc,marker); set(gca,'yscale','log'); hold off;
    xlabel('c_i','FontSize',8); 
    xlim([0 0.12]); ylim([3e-4 2e-2])
    if count == 1, let='a'; else, let='b'; end
   title(['(' let ') d_i=' num2str(d(nn)*1e3) ' mm'],'FontSize',8)
end;
subplot(1,nplots,1); ylabel('y (m)','FontSize',8);

% Mean transport rates 1e6*m^2/s
qSm = mean(qS(ivec,:))*1e6;
qBm = mean(qB(ivec,:))*1e6;
qTotm = qBm + qSm 
sum(qTotm); % Total transport

% Compare with those measured by Hassan & Ribberink (2005)
qTotHR = [4.0 0 13.4 0] % Case K2
