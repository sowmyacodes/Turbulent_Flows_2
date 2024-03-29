clear all; close all; clc

load('out_MatRANS.mat');
%% Initialising

if ~isfield(MatRANS, 'Uf')
    MatRANS.Uf = sqrt(MatRANS.tau0/MatRANS.rho);
end
%% 

% Constants
g = 9.81; % Gravitational acceleration (m/s^2)
rho = 1000; % Fluid density (kg/m^3);
nu = 1e-6; % Fluid kinematic viscosity (m^2/s)
beta = 9/100;

%this is a test change for github

%import values from MatRANS
y = MatRANS.y;
h = y(end);
u_mean = MatRANS.u(end,:); %mean velocity %last time step only
k = MatRANS.k(end,:); %turbulent kinetic energy
k(1) = 10e-4;
w = MatRANS.omega(end,:);
Uf = MatRANS.Uf(end);
L = beta^(-1/4).*(sqrt(k)./w);
%vprime = 1./(k.*3);
vprime = k./3;

RND_mean = 0;
RND_std_dev = 1;


%% Non-dimensional values
%w_nd = w*h/Uf;
Y = y./h;
U_nd = u_mean./Uf;
k_nd = k./(Uf^2);

%vprime_nd = vprime./Uf; %also wrong
%vprime_nd = k_nd./3; %WRONG

vrms_nd = sqrt(vprime)./Uf;

L_nd = L./h;
dT = L_nd./vrms_nd;

%% Parameters and Arrays for storing data

Np = 5; % number of particles
N_it = 5000; %number of iterations for break

NT = 25; %non-dimensional time limit

P = struct('Xp',[],'Yp',[],'Tp',[]); % structure to save particle tracks


%do not prealloacte because then you cant see the final position

% % % % % % Vertical and horizontal positions
% % % % % Yp  = zeros(N_it,1);
% % % % % Xp  = zeros(N_it,1);
% % % % % % Time step size and time values
% % % % % del_T  = zeros(N_it,1);
% % % % % Tp  = zeros(N_it,1);
% % % % % % Streamwise velocity
% % % % % Up_y = zeros(N_it,1);
% % % % % % Vertical and horizontal displacements
% % % % % dY  = zeros(N_it-1,1);
% % % % % dX  = zeros(N_it-1,1);

%boundary conditions
y_bottom = 70*nu/Uf/h;
y_top = 1;

% fig1 = figure();
% xlabel('X position');
% ylabel('Y position');
% hold on
% yline(y_bottom);
% yline(y_top);

fig2 = figure();
xlabel('X position');
ylabel('Y position');
hold on
yline(y_bottom);
yline(y_top);

% fig3 = figure();
% hold on
% xlabel('Mean X position')
% ylabel('Time')

rng(1);


% Yp = zeros();
Y0 = linspace(y_bottom,y_top,Np); 
%% loop over the number of particles
for jj = 1 : Np
    
    ii = 1; %time steps
    
    Xp = zeros();
    Yp = zeros();
    Tp = zeros();
    
    
    %initial conditions
    Xp(1) = 0;
    Tp(1) = 0;
    Yp(1) = Y0(jj);
%     Yp(1) = y_bottom + (y_top-y_bottom)*rand(1);
    
    Up_y(1) = interp1(Y, U_nd, Yp(1)); %interpolate the non-dimensional u against y to find it at the step ii

    
    % loop to calculate one particle track based on the initial position
    %while Tp(ii)<=N_it && ii<=NT %why is the limit like this i dont understand
    while ii<=N_it && Tp(ii)<=NT %this is wrong but why? i do not
    %understand what these limits mean

        ar = randn;
        %ar = normrnd(RND_mean, RND_std_dev); %normal distribution
        del_T = interp1(Y, dT, Yp(ii)); %interpolate the time step across Y
        Tp(ii+1) = Tp(ii) + del_T; %add the time step to the initial time
        
        %initialise the velocity of the particle for the first time step
        if ii > 1
        %if it is not the initial value, it just stays the same --> run later
            Up_y(ii) = up_dy;
        end
        
        %find rms v prime
        vrms_p = interp1(Y, vrms_nd, Yp(ii)); %v rms as a vector in time
        dY = ar*vrms_p*del_T; %change in y velocity

        
        %boundary conditions for rms v prime (reflect the y vector with the
        %same angle as before)
        if (Yp(ii)+ dY) > y_top
            Yp(ii+1) = y_top - (dY -(y_top - Yp(ii)));
        elseif (Yp(ii) + dY) < y_bottom
            Yp(ii+1) = y_bottom + (-dY - (Yp(ii)-y_bottom));
        else
            Yp(ii+1) = Yp(ii) + dY;
        end
        
        % Obtain the streamwise velocity at the next vertical position
        up_dy    = interp1(Y, U_nd, Yp(ii+1));
        % Calculate the streamwise displacement and position
        dX(ii)   = 0.5 * del_T *(up_dy + Up_y(ii));
        %dx = (U_nd.*Yp(ii) + U_nd.*(Yp(ii)+up_dy(ii))).*del_T(ii)./2;
        Xp(ii+1) = Xp(ii) + dX(ii);
        
        
% % %         if Time<=N_it %number of iterations as a break
% % %             break
% % %         end
        
        ii = ii + 1; % update time step at the end of the while loop
        
    end
    
    
    % Store particle track in variable P using e.g. P(jj).Xp = ...;
    P(jj).Yp = Yp(1:(ii));
    P(jj).Xp = Xp(1:(ii));
    P(jj).Tp = Tp(1:(ii));
% % %     X_mean(jj) = sum(Xp)/length(Xp);
% % %     X_mean(ii) = mean(Xp(ii));
   

%     figure(fig1);
%     plot(Xp(1:(ii)),Yp(1:(ii)),'-');

%     figure(fig2);
%     plot(Xp(end), Yp(end),'o');

%     figure(fig3);
%     plot(Xp,Tp);
end


%% Analyse results
times = 5:5:25;
Ntimes = size(times,2);

% Preallocate for auxiliary calculations
XP_store = zeros(Np,Ntimes);
YP_store = zeros(Np,Ntimes);

% Interpolate for each particle at the target time
for jj = 1:Ntimes
    for ii = 1:Np

        XP_store(ii,jj) = interp1(P(ii).Tp, P(ii).Xp, times(jj));
        YP_store(ii,jj) = interp1(P(ii).Tp, P(ii).Yp, times(jj));

    end

end

%%
% Calculate the mean
Xmean = mean(XP_store);
Ymean = mean(YP_store);

% Calculate the variance
Xvar  = var(XP_store);
Yvar  = var(YP_store);
%Xvar = (XP_store-Xmean).^2;
%%

figure(4)
plot(times,Xmean,'k--');
hold on
plot(times,Xvar,'k-o');
% ylabel('$\overline{X}$ [-]','Interpreter','latex');
xlabel('Time [-]','Interpreter','latex');

p1 = polyfit(times,Xmean,1);
p2 = polyfit(times,Xvar,1);
xmean_fit = polyval(p1,times);
xvar_fit = polyval(p2,times);
plot(times,xmean_fit);
plot(times,xvar_fit);

legend ('$\overline{X}$ [-]','$(X - \overline{X})^{2}$ [-]','$\overline{X}$ fit','$(X - \overline{X})^{2}$ fit','Interpreter','latex');


% figure(5)
% plot(times,Xvar);
% ylabel('$(X - \overline{X})^{2}$ [-]','Interpreter','latex');
% xlabel('Time [-]','Interpreter','latex');



