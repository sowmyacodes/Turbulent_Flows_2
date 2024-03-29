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
random_selection = 0;
max_Tit = 1000;
max_T   = 25;
Np      = 1000; % number of particles
yplus_bottom = 70;  % Non-dimensional height of the bottom limit
% Plotting
plot_trajectories = 0;
plot_dispersant_cloud = 1;

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
k_last(1,1) = 10e-10;                % Change the 0 to a small value
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


P = struct('Xp',[],'Yp',[],'Tp',[]); % structure to save particle tracks

y_bottom = yplus_bottom *  nu/ Uf(end); % Dimensional bottom height
y_top = h;                              % Dimensional top height

yp_bottom = y_bottom/h;                 % Non-dimensional bottom height
yp_top    = y_top/h;                    % Non-dimensional top height

% Calculate the length scale 
l = beta_star0^(-0.25) * sqrt(k_last) ./ omega_last;

Lp = l ./ h; % Non-dimensionalise the length scale

rms_Vprime = sqrt(1 / 3 .* kp_last);  % Calculate the v' for each height
rms_vprime = sqrt(1 / 3 .* k_last);   % Same but dimensional

% Calculate the time increment vector (the time step size depends on the
% vertical location)
Delta_t = Lp ./ rms_Vprime; % Non-dimensional
delta_t = l  ./ rms_vprime; % Dimensional

% Plot the time step size against the vertical dimension
plot_me_(Delta_t,y,1, 'Step Size for each vertical position', ...
         '$\Delta t$', '$y/h$', 'k-o');

%% Preallocate arrays to store the data in each loop iteration
% Vertical and horizontal positions
Yp  = zeros(max_Tit,1); % Non-dimensional
Xp  = zeros(max_Tit,1);
yp  = zeros(max_Tit,1); % Dimensional
xp  = zeros(max_Tit,1);
% Time step size and time values
Dt  = zeros(max_Tit,1); % Non-dimensional
Tp  = zeros(max_Tit,1);
dt  = zeros(max_Tit,1); % Dimensional
tp  = zeros(max_Tit,1);
% Streamwise velocity
Upy = zeros(max_Tit,1); % Non-dimensional
upy = zeros(max_Tit,1); % Dimensional
% Vertical and horizontal displacements
Dy  = zeros(max_Tit-1,1); % Non-dimensional
Dx  = zeros(max_Tit-1,1);
dy  = zeros(max_Tit-1,1); % Dimensional
dx  = zeros(max_Tit-1,1);

%% Open figure to plot particles trajectories
if plot_trajectories
    figure(50)
    hold on
end

if plot_dispersant_cloud
    figure(51)
    hold on
end


%% Calculate the particles' trajectories
% Set a seed for the random number generator in order to compare results
rng(1)

if ~random_selection
    y_initial = linspace(y_bottom, y_top, Np);
end

% loop over the number of particles
for jj = 1 : Np

    % Time counter
    ii = 1; 

    if random_selection
        % initialize particle start position with a random distribution
        ran_num   = rand(1);    % Generate a random number to get the starting 
                                % point of the particle
        % Scale the random number to the flume 
        y0  = y_bottom + (y_top - y_bottom)*ran_num; % Dimensional
    else
        y0  = y_initial(jj);
    end
                                                  
    % Store the first value
    yp(1) = y0;
    % Initialize time variable 
    time = 0;

    % loop to calculate one particle track based on the initial position
    while ii <= max_Tit && time <= max_T * h / Uf 
                                % Limit the number of iterations and 
                                         % the maximum time
        
        % Interpolate to get the time step
        dt(ii)   = interp1(y, delta_t, yp(ii));
        tp(ii+1) = time + dt(ii);
        time     = tp(ii+1);

        % Interpolate to get the non-dimensional streamwise velocity at y
        if ii == 1      % Only for the first time step
            upy(ii) = interp1(y, u, yp(ii));
        else            % Otherwise just assign the previous value
            upy(ii) = up_dy;
        end

        % Generate random number between -1 and 1 using a normal
        % distribution with mean 0 and stddev 1
        
        a_r = normrnd(0,1);
%         a_r = randn;
        % Do calculations required in each time step 

        % Interpolate to get the V' rms
        vpri_rms = interp1(y,rms_vprime, yp(ii));
        % Calculate the vertical displacement and position
        dy(ii)   = a_r * vpri_rms * dt(ii);
        yp(ii+1) = yp(ii) + dy(ii); 
        % Check that the vertical position is not out of bounds
%         if jj == 991
%             disp('something')
%         end
        
        while yp(ii+1) > y_top || yp(ii+1) < y_bottom
            if yp(ii+1) > y_top
                auxY = yp(ii+1) - y_top;
                yp(ii+1) = y_top - auxY;
                dy(ii) = yp(ii+1) - yp(ii);
            elseif yp(ii+1) < y_bottom
%                 if jj == 991
%                     disp('im in the if branch')
%                 end
                auxY = y_bottom - yp(ii+1);
                yp(ii+1) = y_bottom + auxY;
                dy(ii) = yp(ii+1) - yp(ii);
            end
        end

        % Obtain the streamwise velocity at the next vertical position
        up_dy    = interp1(y, u, yp(ii+1));
        % Calculate the streamwise displacement and position
        dx(ii)   = 0.5 * dt(ii)*(up_dy + upy(ii));
        xp(ii+1) = xp(ii) + dx(ii);
        
        ii = ii + 1; % update time step at the end of the while loop
        
    end

    % Plot particle track
    if plot_trajectories
        figure(50)
        plot(xp(1:ii),yp(1:ii),'--')
    end

    if plot_dispersant_cloud
        figure(51)
        plot(xp(ii),yp(ii), 'o')
    end

    % Store particle track in variable P using e.g. P(jj).Xp = ...; 
    % Using non-dimensional variables
    P(jj).Yp  = yp(1:ii)/h; 
    P(jj).Xp  = xp(1:ii)/h;
    P(jj).Tp  = tp(1:ii)*Uf/h;
    P(jj).Dx  = dx(1:ii)/h;
    P(jj).Dy  = dy(1:ii)/h;
    P(jj).Dt  = dt(1:ii)*Uf/h;
    P(jj).upy = upy(1:ii)/Uf;

    
end

%% Analyse results for T = 10s

Tobj = 10;
[XP_mean, YP_mean, XP_var, YP_var] = get_particle_stats(P,Tobj);


%% Do the same for a series of time targets

Tobj   = [5 10 15 20 25];
NTimes = size(Tobj,2);

% Preallocate arrays to store the mean values
XP_means = zeros(NTimes,1);
YP_means = zeros(NTimes,1);
XP_vars  = zeros(NTimes,1);
YP_vars  = zeros(NTimes,1);

for ii = 1:NTimes

    [XP_means(ii), YP_means(ii), XP_vars(ii), YP_vars(ii)] = ...
        get_particle_stats(P,Tobj(ii));

end

plot_me_(Tobj, XP_means, 1,'Streamwise mean position of particles', ...
    'Time [-]', '$\overline{X}$', 'k-')
plot_me_(Tobj, YP_means, 1,'Verticle mean position of particles', ...
    'Time [-]', '$\overline{Y}$', 'k-')
plot_me_(Tobj, XP_vars, 1,'Streamwise position variance of particles', ...
    'Time [-]', 'Variance $\overline{(X - \overline{X})^{2}}$', 'k-')
plot_me_(Tobj, YP_vars, 1,'Verticle position variance of particles', ...
    'Time [-]', 'Variance $\overline{(Y - \overline{Y})^{2}}$', 'k-')
