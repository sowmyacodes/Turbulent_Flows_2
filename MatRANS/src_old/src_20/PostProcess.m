% Post-processor for MaRANS Matlab k-omega model
% Programmed by David R. Fuhrman
clear all;

% Load data from simulation
OutFileName = 'Case7.mat'; 
load(OutFileName);

% Choose which variables to plot/animate (0: Off, 1: On)
pvec_ts = [1   0         1     1  1  0  0  0   ]; % Plot control for time series 
%         [U_0 tau_b/rho theta qB qS qT cb cavg] 
ihold = 1; % Time series: 0: hold off, 1: hold on
pvec_an = [1 0 1 0 0    0       1 1  ]; % Plot control for animation
%         [U V K W nu_T tau/rho C U*C]
yscale = 'linear'; % y-scale for animation (either 'linear' or 'log');
ymax = max(y)/2; % Maximum vertical elevation to plot
marker = 'b-'; % Marker for time series plots
nrepeat = 3; % No. of times to repeat animation


% Choose time indices to plot / animate
%ivec = 1:length(t); % Plot all saved time indices
nT = 36; ivec = length(t)-nT:length(t)-1; % Plot last period
%ivec = ivec(1:10); % Uncomment to stop at omega*t=90 deg.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot selected time series
np_ts = sum(pvec_ts); % Number of sub-plots
figure(1)
sp = 0;
if pvec_ts(1) % Free stream velocity
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),U(ivec,length(U(1,:))),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('U_0 (m/s)')
end
if pvec_ts(2) % Bed shear stress
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),tau_b(ivec)/rho,marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('\tau_b/\rho (m^2/s^2)')
end
if pvec_ts(3) % Shields parameter 
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),Shields(ivec),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('\theta')
end
if pvec_ts(4) % Bed load transport 
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),qB(ivec),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('q_B (m^2/s)')
end
if pvec_ts(5) % Suspended load transport 
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),qS(ivec),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('q_S (m^2/s)')
end
if pvec_ts(6) % Total transport qT = qB + qS
  qT = qB + qS;
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),qT(ivec),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('q_T (m^2/s)')
end
if pvec_ts(7) % Reference concentration cb
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),C(ivec,1),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('c_b')
end
if pvec_ts(8) % Average of c over column
  for i = 1:length(t)
    Cavg(i) = trapz(yc,C(i,:))./yc(end);
  end
  if ihold == 1; hold on; else; hold off; end
  sp = sp+1; subplot(np_ts,1,sp); plot(t(ivec),Cavg(ivec),marker); 
  xlim([t(ivec(1)) t(ivec(end))]); ylabel('c_{avg}')
end
xlabel('t (s)');
hold off;

% Determine bounds
Umax = max(max(abs(U(ivec,:))));
Vmax = max(max(abs(V(ivec,:))));
Kmax = max(max(abs(K(ivec,:))));
Wmax = max(max(abs(W(ivec,:))));
nu_Tmax = max(max(abs(nu_T(ivec,:))));
taumax = max(max(abs(tau(ivec,:)./rho)));
Cmax = max(max(abs(C(ivec,:))));
for i = 1:length(t)
  UC(i,:) = U(i,ib:end).*C(i,:);
end
UCmax = max(max(abs(UC(ivec,:))));


% Animate
np_an = sum(pvec_an); % Number of subplots
figure(2), clf
for r = 1:nrepeat
  disp(['Repeat of animation no. ' int2str(r)]);
  for i = ivec
  
    % Current plot
    sp = 0; % Initialize sub-plot counter
    if pvec_an(1) % U
      sp=sp+1; subplot(1,np_an,sp); plot(U(i,:),y); 
      xlim(Umax.*[-1 1]); ylim([0 ymax]); xlabel('U (m/s)');
      set(gca,'yscale',yscale);
    end
    if pvec_an(2) % V 
      sp=sp+1; subplot(1,np_an,sp); plot(V(i,:),y); 
      xlim((Vmax+1e-16).*[-1 1]); ylim([0 ymax]); xlabel('V (m/s)');
      set(gca,'yscale',yscale);
    end
    if pvec_an(3) % K
      sp=sp+1; subplot(1,np_an,sp); plot(K(i,:),y); 
      xlim(Kmax.*[0 1]); ylim([0 ymax]); xlabel('K (m^2/s^2)');
      set(gca,'yscale',yscale);
    end  
    if pvec_an(4) % omega
      sp=sp+1; subplot(1,np_an,sp); plot(W(i,:),y); 
      xlim(Wmax.*[0 1]); ylim([0 ymax]); xlabel('\omega (1/s)');
      set(gca,'yscale',yscale);
    end 
    if pvec_an(5) % nu_T
      sp=sp+1; subplot(1,np_an,sp); plot(nu_T(i,:),y); % nu_T 
      hold on; plot(K(i,:)./W(i,:),y,'r--'); hold off; legend('\nu_T','K/\omega',2) % Also add K/omega
      xlim(nu_Tmax.*[0 1]); ylim([0 ymax]); xlabel('\nu_T (m^2/s)');
      set(gca,'yscale',yscale);
    end  
    if pvec_an(6) % tau/rho
      sp=sp+1; subplot(1,np_an,sp); plot(tau(i,:)./rho,y); 
      xlim(taumax.*[-1 1]); ylim([0 ymax]); xlabel('\tau/\rho (m^2/s^2)');
      set(gca,'yscale',yscale);
    end
    if pvec_an(7) % C
      sp=sp+1; subplot(1,np_an,sp); plot(C(i,:),yc); 
      xlim(Cmax.*[0 1]); ylim([0 ymax]); xlabel('C');
      set(gca,'yscale',yscale);
    end   
    if pvec_an(8) % C*U
      sp=sp+1; subplot(1,np_an,sp); plot(UC(i,:),yc); 
      xlim(UCmax.*[-1 1]); ylim([0 ymax]); xlabel('C*U');
      set(gca,'yscale',yscale);
    end
    
    % Display y-axis label, and time
    subplot(1,np_an,1); ylabel('y (m)'); title(['t = ' num2str(t(i)) ' s']);  
    pause(0.15);
    drawnow % Update plots
  end
end


% Plot period-averaged profiles
np_an = sum(pvec_an); % Number of subplots
figure(3), clf
   
% Current plot
sp = 0; % Initialize sub-plot counter
if pvec_an(1) % U
  sp=sp+1; subplot(1,np_an,sp); plot(mean(U(ivec,:)),y); 
  ylim([0 ymax]); xlabel('U (m/s)');
  set(gca,'yscale',yscale);
end
if pvec_an(2) % V 
  sp=sp+1; subplot(1,np_an,sp); plot(mean(V(ivec,:)),y); 
  ylim([0 ymax]); xlabel('V (m/s)');
  set(gca,'yscale',yscale);
end
if pvec_an(3) % K
  sp=sp+1; subplot(1,np_an,sp); plot(mean(K(ivec,:)),y); 
  ylim([0 ymax]); xlabel('K (m^2/s^2)');
  set(gca,'yscale',yscale);
end  
if pvec_an(4) % omega
  sp=sp+1; subplot(1,np_an,sp); plot(mean(W(ivec,:)),y); 
  ylim([0 ymax]); xlabel('\omega (1/s)');
  set(gca,'yscale',yscale);
end 
if pvec_an(5) % nu_T
  sp=sp+1; subplot(1,np_an,sp); plot(mean(nu_T(ivec,:)),y); % nu_T 
  hold on; plot(mean(K(ivec,:))./mean(W(ivec,:)),y,'r--'); hold off; legend('\nu_T','K/\omega',2) % Also add K/omega
  ylim([0 ymax]); xlabel('\nu_T (m^2/s)');
  set(gca,'yscale',yscale);
end  
if pvec_an(6) % tau/rho
  sp=sp+1; subplot(1,np_an,sp); plot(mean(tau(ivec,:))./rho,y); 
  ylim([0 ymax]); xlabel('\tau/\rho (m^2/s^2)');
  set(gca,'yscale',yscale);
end
if pvec_an(7) % C
  sp=sp+1; subplot(1,np_an,sp); plot(mean(C(ivec,:)),yc); 
  ylim([0 ymax]); xlabel('C');
  set(gca,'yscale',yscale);
end   
if pvec_an(8) % C*U
  sp=sp+1; subplot(1,np_an,sp); plot(mean(UC(ivec,:)),yc); 
  ylim([0 ymax]); xlabel('C*U');
  set(gca,'yscale',yscale);
end
subplot(1,np_an,1); ylabel('y (m)'); title(['Averaged']);  % Display y-axis label


% Contour plot of suspended sediment concentration
figure(4), clf
[cs,hc] = contour( t(ivec), yc, C(ivec,:)',[0:0.025:0.6]);
set(gca,'yscale','log');
clabel(cs,hc); % Label contours
xlabel('t (s)'); ylabel('y (m)'); title('C(y,t)');
ylim([0 ymax]);
drawnow % Update plots

% Calculate average sediment fluxes
qBm = mean(qB(ivec))
qSm = mean(qS(ivec))
qTm = qBm + qSm;


