%% MatRANS Simulations for the turbulent flow in a flume
% Denmark Technical University
% 41129 Turbulent Flows · Assignment 2
% Authors: 
%  · Sowmya Srinivasan Iyer
%  · Videep Goverdhan Kamath
%  · Yann Birnie Scott 
%  · Carlos Perez Moreno
% 
% Version 1
% Date: 08/10/2021
% -------------------------------------------------------------------------
% MODIFICATIONS
% Version 2. 09/10/2021
% · Correct the mirroring of particles 
% · Store the data in the structure and plot it in the loop
% _________________________________________________________________________
% Version 3. 13/10/2021
% · Include comparison with number of particles
% · Correct error in DY formula
% · Gather one-particle-analysis in an external function:
%                                              'multi_particle_walk.m'
% _________________________________________________________________________
% Version 4. 14/10/2021
% · Auto-print figures in pdf
% -------------------------------------------------------------------------
%% FLOW PARAMETERS AND COORDINATES REPRESENTATION
%
%  y/\ 
%   |
%   |_____________________________________________________________ y_top
%   |      |----->|                            |
%   |      |      /                            |
%   |      |---->/                             |
%   |      |    /                              | h
%   |      |-->/                               |
%   |      |  /                                |
%   |      |_/                                 |
%   |------------------------------------------|------------------ y+
%   |      ---> Uf                             |
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
%% VARIABLES DICTIONARY
%
% Np              ----> Number of particles for the analysis
% Np_min          ----> Exponent of the min 10th power for comparison
% Np_max          ----> Exponent of the max 10th power for comparison
% NumNp           ----> Number of different Np for no of particle analysis
% Np_vec          ----> Array with the numbers of particles per case
% yplus_bottom    ----> Non-dimensional bottom limit for the calculations
% n_t             ----> Number of time steps until RANS convergence [-]
% ny              ----> Number of vertical grid points [-]
% y               ----> Dimensional vertical position [m]
% h               ----> Height of the flume [m]
% Y               ----> Non-dimensionalise vertical position [-]
% Uf_vec          ----> Frictional speed at each RANS time iteration [m/s]
% Uf              ----> Frictional speed at last RANS time iteration [m/s]
% T_comp          ----> Non-dimensional time for 
%

clear; close all; clc

%% Parameters
% Preprocess of data
beta_star0      = 0.09;
file_outmatrans = 'out_MatRANS.mat';
% Track particles
random_selection = 0; % Flag for using random selection of starting points
max_Tit = 5000; % Limit of the time iterations to reach maximum allowed time
max_T   = 25;   % Maximum allowed time (non-dimensional)
Np_min  = 1;   % Exponent of the min 10th power for no of particle analysis
Np_max  = 5;   % Exponent of the max 10th power for no of particle analysis
NumNp   = 25;  % Number of different Np for no of particle analysis
Np      = 680;     % number of particles for standard analysis
yplus_bottom = 70;  % Non-dimensional height of the bottom limit
Tobj   = [5 10 15 20 25]; % Target non-dimensional times for the statistical analysis
T_comp = 10;   % Non-dimensional time to analyse the influence of Np
% Plotting
plot_trajectories = 0;     % Flag for plotting trajectories
plot_dispersant_cloud = 0; % Flag for plotting the dispersant cloud at maxT
% File comparison Np data and re-run
file_comparison_np = 'comparison_np_data.mat';
rerun = 0;
Np_showpath = 5; % Exponent of the number of particles to show the path

%% Load the results file
load(file_outmatrans);

%% Show the path of some few particles 

[Ppaths,~] = multi_particle_walk(MatRANS, Np_showpath, beta_star0, ...
                    max_Tit, max_T, yplus_bottom, T_comp, ...
                    random_selection, 1, 0);

%% Calculate the results

% Define an array with the target number of particles
Np_vec = logspace(Np_min, Np_max, NumNp);

if rerun == 1
    % Analysis for a series of number of particles
    % Preallocate arrays to store the values to plot
    Xmeans = zeros(NumNp,1);
    Xvars  = zeros(NumNp,1);
    
    
    for ii = 1:NumNp
        % Select the number of particles
        NPart = Np_vec(ii); 
        % Call the particle walker to obtain the statistics
        [~, stats] = multi_particle_walk(MatRANS, NPart, beta_star0, ...
                    max_Tit, max_T, yplus_bottom, T_comp, ...
                    random_selection, plot_trajectories, ...
                    plot_dispersant_cloud);
        % Store the mean and variance
        Xmeans(ii) = stats.Xmeans;
        Xvars(ii)  = stats.Xvars;
    
    end
    % Save the file to avoid rerun unnecessarily
    save('comparison_np_data.mat', 'Xmeans','Xvars');

else % If you don't want to rerun

    % Load the existing data to compare
    load('comparison_np_data.mat');

end

%% Compare the values of the maximum number of points to the different

XMeanLast = Xmeans(end);
XVarLast  = Xvars(end);

RelError_mean = abs(Xmeans - XMeanLast) / XMeanLast;
RelError_var  = abs(Xvars - XVarLast) / XVarLast;
RelError      = table(Np_vec', RelError_mean, RelError_var);

table2latex(RelError,'table_relerror_latex');

%% Open figures for plotting the results of the comparison

figure("Name","Comparison of mean values for different Np");
semilogx(Np_vec, Xmeans,'k-o');
grid on
xlabel('Number of particles, $N_{P}$', 'Interpreter','latex')
ylabel('$\overline{X}$', 'Interpreter','latex')
title('\textbf{Comparison of mean values at $T = 10$}', ...
    'Interpreter','latex')
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('mean_over_np_T10','-dpdf');


figure("Name","Comparison of variance for different Np");
semilogx(Np_vec, Xvars,'k-o')
grid on
xlabel('Number of particles, $N_{P}$', 'Interpreter','latex')
ylabel('$\overline{(X - \overline{X})^{2}}$', 'Interpreter','latex')
title('\textbf{Comparison of variance values at $T = 10$}', ...
    'Interpreter','latex')
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('pseudovariance_over_np_T10','-dpdf');



%% Analysis along time

[P, stats] = multi_particle_walk(MatRANS, Np, beta_star0, max_Tit, max_T...
             , yplus_bottom, Tobj, random_selection, plot_trajectories, ...
             plot_dispersant_cloud);

%% Plot the results

XP_means = stats.Xmeans;
XP_vars  = stats.Xvars;
YP_means = stats.Ymeans;
YP_vars  = stats.Yvars;

fig_xmean = plot_me_(Tobj, XP_means, 1,'Streamwise mean position of particles', ...
                    'Time [-]', '$\overline{X}$', 'k-');

% plot_me_(Tobj, YP_means, 1,'Verticle mean position of particles', ...
%     'Time [-]', '$\overline{Y}$', 'k-'); % Not so interesting to plot 

fig_xvar  = plot_me_(Tobj, XP_vars, 1,'Streamwise position variance of particles', ...
    'Time [-]', '$\overline{(X - \overline{X})^{2}}$', 'k-');

% plot_me_(Tobj, YP_vars, 1,'Verticle position variance of particles', ...
%     'Time [-]', 'Variance $\overline{(Y - \overline{Y})^{2}}$', 'k-');

%% Calculate the dispersion coefficient

% Ultimate mean velocity of the cloud

p1 = polyfit(Tobj,XP_means,1);
p2 = polyfit(Tobj,XP_vars,1);
xmean_fit = polyval(p1,Tobj);
xvar_fit = polyval(p2,Tobj);

% Plot the fit of the mean
figure(fig_xmean);
hold on
plot(Tobj,xmean_fit, 'k--');
legend('$\overline{X}$','$\overline{X}$ fit', 'interpreter', 'latex', ...
    'Location','best')
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('mean_over_time_n_fit','-dpdf');

% Plot the fit of the variance
figure(fig_xvar);
hold on
plot(Tobj,xvar_fit, 'k--');
legend('$(X - \overline{X})^{2}$','$(X - \overline{X})^{2}$ fit', ...
    'interpreter', 'latex', 'Location','best')
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('variance_over_time_n_fit','-dpdf');

% Calculate the longitudinal dispersion coefficient
D1 = 1 / 2 * p2(1);
fprintf('The longitudinal dispersion coefficient is D1 = %4.4f\n', D1);
