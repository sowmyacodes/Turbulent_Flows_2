% Convert the Mercator elevation time series to a velocity time series
clear all; close all;

% Input
%h0 = 14; % Average depth (m), computed below based on first measurement
keel = 2.1; % Keel depth (m)
g = 9.81; % Gravitational acceleration (m/s^2)

% Load the digitized elevation data (obtained using GetData.m)
load Mercator.mat;

% Create time and true depth
t = x; % Time (s)
h = y + keel; % Account for the additioanal keel depth
subplot(3,1,1); plot(t,h); ylabel('depth (m)');
xlim([0 t(end)]); set(gca,'xtick',[0:1800:5000]);
grid on; box on;

% Convert to free surface elevation
h0 = h(1); % Estimated depth (m)
eta = h - h0; % Free surface elevation
subplot(3,1,2); plot(t,eta); ylabel('\eta (m)');
xlim([0 t(end)]); set(gca,'xtick',[0:1800:5000]);
grid on; box on;

% Convert free surface elevation to velocities
u = sqrt(g/h0).*eta; % Velocity (m/s), based on shallow-water scaling
subplot(3,1,3); plot(t,u); ylabel('u (m/s)');
xlim([0 t(end)]); set(gca,'xtick',[0:1800:5000]);
xlabel('t (s)');
grid on; box on;

% Output key variables to screen
h0
umax = max(abs(u)) % Returns maximum velocity

% Save data (also add zeros for initial conditions)
t_ud = [0 t']; U_ud = [0 u'];  % Save model input variables (for iwave=10)
save MercatorVelocity.mat t_ud U_ud
