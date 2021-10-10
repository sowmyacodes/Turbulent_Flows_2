clear all; close all; clc
load('out_MatRANS.mat');
%% Initialising

if ~isfield(MatRANS, 'Uf')
    MatRANS.Uf = sqrt(MatRANS.tau0/MatRANS.rho);
end

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

% vprime = 1./(k.*3);
vprime = k./3;
RND_mean = 0;
RND_std_dev = 1;


%% Non-dimensional values
%w_nd = w*h/Uf;
Y = y./h;
U_nd = u_mean./Uf;
k_nd = k./(Uf^2);

vprime_nd = vprime./Uf;
%vprime_nd = k_nd./3; %WRONG

vrms = sqrt(vprime_nd);

L_nd = L./h;
dT = L_nd./vrms;

%% Parameters and Arrays for storing data

Np = 100; % number of particles
N_it = 50000; %number of iterations for break
NT = 25; %time

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


%% loop over the number of particles
for jj = 1 : Np
    
    ii = 1; %time steps
    
    %initial conditions
    Xp(1) = 0;
    Tp(1) = 0;
    Yp(1) = y_bottom + (y_top-y_bottom)*rand(1);
    
    Up_y(1) = interp1(Y, U_nd, Yp(1)); %interpolate the non-dimensional u against y to find it at the step ii

    
    % loop to calculate one particle track based on the initial position
    while Tp(ii)<=NT && ii<=N_it 
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
        vrms_p = interp1(Y, vrms, Yp(ii)); %v rms as a vector in time
        
        dY = ar.*vrms_p.*del_T; %change in y velocity
        
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
    X_mean(jj) = sum(Xp)/length(Xp);
    
    figure(1)
    plot(Xp(1:(ii)),Yp(1:(ii)),'-');
    xlabel('X position');
    ylabel('Y position');
    hold on
    yline(y_bottom);
    yline(y_top);

    figure(2)
    plot(Xp(end-1), Yp(end-1),'o');
    xlabel('X position');
    ylabel('Y position');
    hold on
    yline(y_bottom);
    yline(y_top);


%    figure(3)
%    plot(Xp,Tp);
%    hold on
%    xlabel('Mean X position')
%    ylabel('Time')
    
end


%% Analyse results




