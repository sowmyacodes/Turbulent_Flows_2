clear all; close all; clc
%% MatRANS Simulations for the turbulent flow in a flume
% Version 1
% Author: Carlos Perez Moreno
% Date: 08/10/2021
% Denmark Technical University
% 41129 Turbulent Flows · Assignment 2
% -------------------------------------------------------------------------
% MODIFICATIONS
% Version 2. 09/10/2021
% · Correct the mirroring of particles 
% · Store the data in the structure and plot it in the loop
% _________________________________________________________________________
% Version 3. dd/mm/yyyy
% · 
% -------------------------------------------------------------------------
%% FLOW REPRESENTATION
%
%  y/\ 
%   |
%   |_____________________________________________________________
%   |      |------|                            |
%   |      |      /                            |
%   |      |-----/                             |
%   |      |    /                              | h
%   |      |---/                               |
%   |      |  /                                |
%   |      |_/                                 |
%   |------------------------------------------|------------------ y+
%   |                                          |
%   |_______________________________________________________________
%   --------------------------------------------------------------> x 
%
%% PARAMETERS DICTIONARY
%
% beta_star0      ----> Length scale parameter
% file_outmatrans ----> Path to file containing output of RANS model
% max_Tit         ----> Maximum number of iterations in the Time loop
% max_T           ----> Non-dimensional Time limit for the Time loop
%
%
%% VARIABLES DICTIONARY
%
% n_t             ----> Number of time steps until RANS convergence [-]
% ny              ----> Number of vertical grid points [-]
% y               ----> Dimensional vertical position [m]
% h               ----> Height of the flume [m]
% Y               ----> Non-dimensionalise vertical position [-]
% Uf_vec          ----> Frictional speed at each RANS time iteration [m/s]
% Uf              ----> Frictional speed at last RANS time iteration [m/s]
%
%

%% Parameters
% Preprocess of data
beta_star0      = 0.09;
file_outmatrans = 'out_MatRANS.mat';
% Track particles
max_Tit = 50;
max_T   = 25;
Np      = 1000; % number of particles

%% Load the results file
load(file_outmatrans);

% If the frictional velocity Uf is not included in the output calculate it
% Probably unnecesary for our case but nevermind
if ~isfield(MatRANS, 'Uf')
    MatRANS.Uf = sqrt(MatRANS.tau0/MatRANS.rho);
end

%% Load the MatRANS Output Data

% Number of time simulations and number of grid points
n_t = MatRANS.n_t; 
ny  = MatRANS.n_y;

% Vertical position
y = MatRANS.y;      % Dimensional vertical position
h = y(end);         % Height of the flume [m]
Y = y/h;            % Non-dimensionalise vertical direction

% Frictional speed
Uf_vec = MatRANS.Uf;
Uf     = Uf_vec(end);

% Streamwise speed
u_total = MatRANS.u; 
u       = u_total(end,:);
U       = u ./ Uf;

% Turbulent kinetic energy
k           = MatRANS.k; 
k_last      = k(end, :);
k_last(1,1) = 10e-4;                % Change the 0 to a small value
kp_last     = k_last ./ (Uf^2);     % Non-dimensionalise with Uf

% Turbulent dissipation rate
omega            = MatRANS.omega; 
omega_last       = omega(end,:);
% omega_last(1,1)  = 10^-8; % I dont know if this is needed

nu_t = MatRANS.nu_t;
tau0 = MatRANS.tau0;

nu   = MatRANS.nu; 
rho  = MatRANS.rho; 
k_s  = MatRANS.k_s;
h_m  = MatRANS.h_m; 
T    = MatRANS.T; 
U0m  = MatRANS.U0m;

% ... setup parameters other relevant parameters



P = struct('Xp',[],'Yp',[],'Tp',[]); % structure to save particle tracks

yplus_bottom = 70;  % Non-dimensional height of the bottom limit

y_bottom = yplus_bottom *  nu/ Uf(end); % Dimensional bottom height
y_top = h;                              % Dimensional top height

yp_bottom = y_bottom/h;                 % Non-dimensional bottom height
yp_top    = y_top/h;                    % Non-dimensional top height



% Calculate the length scale 
l = beta_star0^(-0.25) * sqrt(k_last) ./ omega_last;

Lp = l ./ h; 

rms_vprime = sqrt(1 / 3 .* kp_last);  % Calculate the v' 

% Calculate the time increment vector (the time step size depends on the
% vertical location)
Delta_t = Lp ./ rms_vprime;

% Plot the time step size against the vertical dimension
plot_me_(Delta_t,y,1, 'Step Size for each vertical position', ...
         '$\Delta t$', '$y/h$', 'k-o');

%% Preallocate arrays to store the data in each loop iteration
% Vertical and horizontal positions
Yp  = zeros(max_Tit,1);
Xp  = zeros(max_Tit,1);
% Time step size and time values
Dt  = zeros(max_Tit,1);
Tp  = zeros(max_Tit,1);
% Streamwise velocity
upy = zeros(max_Tit,1);
% Vertical and horizontal displacements
Dy  = zeros(max_Tit-1,1);
Dx  = zeros(max_Tit-1,1);

%% Open figure to plot particles trajectories
figure(50)
hold on

%% Calculate the particles' trajectories

% loop over the number of particles
for jj = 1 : Np

    % Time counter
    ii = 1; 
    % initialize particle start position
    ran_num   = rand(1);    % Generate a random number to get the starting 
                            % point of the particle
    yp0 = yp_bottom + (yp_top-yp_bottom)*ran_num; % Scale the random number 
                                                  % to the flume 
                                                  
    % Store the first value
    Yp(1) = yp0;
    
    % Initialize time variable 
    time = 0;

    % loop to calculate one particle track based on the initial position
    while ii <= max_Tit && time <= max_T % Limit the number of iterations and 
                                       % the maximum time
        
        % Interpolate to get the time step
        Dt(ii) = interp1(Y, Delta_t, Yp(ii));
        Tp(ii) = time + Dt(ii);
        time   = Tp(ii);

        % Interpolate to get the non-dimensional streamwise velocity at y
        if ii == 1      % Only for the first time step
            upy(ii) = interp1(Y, U, Yp(ii));
        else            % Otherwise just assign the previous value
            upy(ii) = up_dy;
        end

        % Generate random number between -1 and 1 using a normal
        % distribution with mean 0 and stddev 1
        %a_r = normrnd(0,0.1);
        a_r = randn;

        % Do calculations required in each time step 

        % Interpolate to get the V' rms
        vpri_rms = interp1(Y,rms_vprime, Yp(ii));
        % Calculate the vertical displacement and position
        Dy(ii)   = a_r * sqrt(vpri_rms) * Dt(ii);
        Yp(ii+1) = Yp(ii) + Dy(ii); 
        % Check that the vertical position is not out of bounds
        if Yp(ii+1) > y_top/h
            auxY = Yp(ii+1) - y_top/h;
            Yp(ii+1) = Yp(ii+1) - auxY;
            Dy(ii) = Yp(ii+1) - Yp(ii);
        elseif Yp(ii+1) < y_bottom/h
            auxY = y_bottom/h - Yp(ii+1);
            Yp(ii+1) = Yp(ii+1) + auxY;
            Dy(ii) = Yp(ii+1) - Yp(ii);
        end

        % Obtain the streamwise velocity at the next vertical position
        up_dy    = interp1(Y, U, Yp(ii+1));
        % Calculate the streamwise displacement and position
        Dx(ii)   = 0.5 * Dt(ii)*(up_dy + upy(ii));
        Xp(ii+1) = Xp(ii) + Dx(ii);
        
        ii = ii + 1; % update time step at the end of the while loop
        
    end

    % Plot particle track
    figure(50)
    plot(Xp(1:ii-1),Yp(1:ii-1),'--')

    % Store particle track in variable P using e.g. P(jj).Xp = ...; 
    P(jj).Yp = Yp(1:ii-1); 
    P(jj).Xp = Xp(1:ii-1);
    P(jj).Tp = Tp(1:ii-1);

    
end

% Analyse results
