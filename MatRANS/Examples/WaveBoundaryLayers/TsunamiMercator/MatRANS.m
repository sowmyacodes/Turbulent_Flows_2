% Input file for MatRANS simulation

% Clear workspace, figures, and command line
clear all; close all; clc;
GlobalVars; % Declare global variables

% MODEL INPUT ---------------------------------------------------------------

% Choose which terms to include in the model
turb = 1; % 1: turbulence model 0: Off, 1: On (if 0 then a laminar simulation is made), 2: Transitional
streaming = 0; % 0: Leading-order terms only i.e. d/dx=0, 1: second order (streaming) terms added
susp = 0; % Suspended sediment calculation 0: Off, 1: On

% Constants
g = 9.81; % Gravitational acceleration (m/s^2)
rho = 1000; % Fluid density (kg/m^3)
nu = 1e-6; % Fluid kinematic viscosity (m^2/s)

% Output file name
OutFileName = 'out_MatRANS.mat';

% Grid input
h_m = 14.2377; % Height of model depth (m)
ny = 100; % Numbers of grid points in the vertical direction
%y = linspace(0,h_m,ny); % Define grid, linear spacing
y = h_m.*logspace(-5.4,0,ny-1); y = y - y(1); y = [y h_m]; % Define grid, logarithmic stretching


% Wave input 
T = 14*60; % Wave period (s)
c = 1.0; % Celerity (if streaming=1, used to calculate d/dx = -1/c*d/dt) (m/s)
iwave = 10; % Wave signal, 0: Off, 1: 2nd-order Stokes signal, 2: Wave form of Abreu et al. (2010), 3: Solitary wave, 
           % 4: N-wave, 5: Superimposed sech()^2 signals, 10: User-defined signal
if iwave == 1 % Input for second-order Stokes signal, with free stream velocity: U_1m*sin(omega*t) - U_2m*cos(2*omega*t)
  U_1m = 2.432; % Primary amplitude of free stream velocity (m/s)
  U_2m = 0.0; % Amplitude of second harmonic free stream velocity (m/s)
elseif iwave == 2 % Input for wave form of Abreu et al. (2010)
  U_1m = 2.432; % Total amplitude of free stream velocity signal
  r = 0.0; % Controls nonlinearity, recommended range 0<=r<=0.75
         % If r=0, then the signal will be sinusoidal
  phi = -pi/2; % Controls velocity asymmetry (phi=0) or skewness (phi=-pi/2)
elseif iwave == 3 | iwave == 4 % 3: Solitary wave or 4: N-wave
  T=25*60; U_1m = -2.432; % Maximum free stream velocity (m/s)
  t0 = 1770;
elseif iwave == 5 % Superimposed sech()^2 signals
  U_1m = 2.432; % Simply used as a scale e.g. in f_w calculation, etc.
  h0=14; h=14; Umvec = [-3.1 3.80 -1].*(h0/h)^0.25.*sqrt(g/h); % m/s
  Omegavec = [0.00420 0.00530 0.009];
  tnvec = [600 1000 1399];
  t0 = 1715;
elseif iwave == 10 % User defined signal (interpolated based solely on t_ud, U_ud, and Ut_ud variables)
  T = 14*60; omega = 2*pi/T; % Based on approximate tsunami period (s)
  load Data/MercatorVelocity.mat; % Loads t_ud (s) and U_ud (m/s)
  U_1m = max(U_ud) % Use max as velocity scale calculations
  %U_ud = U_1m*sech(omega.*(t_ud - T)).^2; % Free stream velocity (m/s), single wave
  %U_ud = U_1m*sech(omega.*(t_ud - T)).^2 - U_1m*sech(omega.*(t_ud-1.5*T)).^2; % N-wave (s)
  Ut_ud = gradient(U_ud,t_ud); % Compute free stream acceleration (m/s^2)
  subplot(2,1,1); plot(t_ud,U_ud); ylabel('u (m/s)');
  subplot(2,1,2); plot(t_ud,Ut_ud); ylabel('du/dt (m/s^2)'); xlabel('t (s)');
end

% Bed slope effects (applied if S~=0)
S = 0; % -dh/dx, modifies critical Shields parameter in up/down hill directions
iconv = 0; % Converging-diverging effects 0: Off, 1: On
depth = h_m; % Physical depth or tunnel half-depth

% Current input
Px = 0; % Constant added to 1/rho*dpdx (for no current, just set to zero) (m/s^2)
% For pure current, this will equal -Uf^2/h_m at steady state

% Time input
dt = 4; % Output time interval (s)
t1 = dt; % First time (after t=0) where the output should be saved (s)
t2 = 3300; % End time (s)
t_span = [0:dt:t2]; % Vector of time levels where output will be saved (s)

% Sediment input
d = 0.3/1000; % Grain size (m)
k_N = 2.5*d; % Nikuradse's equivalent sand grain roughness (m)
b = 2*d; % Reference height (m)
s = 2.65; % Relative sediment density
Shields_c0 = 0.045; % Critical Shields parameter (flat bed)
mu_d = 0.65; % Dynamic friction coefficient
phi_s = 32; % Angle of repose (deg)
w_s = 0.001; % Fall velocity (if 0, then calculated empirically) (m/s)

% Method for suspended sediment concentration boundary condition (only used is susp=1)
icb = 4; % 1: Zyserman-Fredsøe cb, 2: O'Donoghue-Wright cb, 3: Engelund-Fredsøe cb, 4: Einstein cb
% 11: Pickup function (van Rijn), 12: Pickup function based on Engelund-Fredsøe cb
beta_s = 2.0; % Sediment diffusivity (eps_s=beta_s*nu_T) (also, inverse of Prandtl-Schmidt number, if 0: Use van Rijn (1984) correction
iextrap = 1; % Use extrapolation to avoid overloading: 0: Off, 1: On

% Hindered settling of suspended sediment at large concentrations
Hind_set = 0; % Hindered settling [w_s = w_s0*(1-c)^nhs] 0: Off, 1: On
nhs = 0; % If nhs=0, then it is calculated according to Richardson and Zaki (1954)

% Turbulence supression
iturbsup = 0; % Suppression terms 0: Off, 1: On
sigma_p = 0.7; % Additional closure coefficient
n = 0.4; % Porosity (used to estimate maximum fluid-sediment density)	 
	 
% Method for computing bed-load sediment transport
BedLoad = 3; % 1: Meyer-Peter Muller, 2: Engelund-Fredsoe (7.59), 
             % 3: Engelund-Fredsoe (7.54) + (7.58), 4: Nielsen (1992)

% Wall boundary condition for turbulent kinetic energy k
ikbc = 0; % 0: k = 0, 1: dk/dy=0 wall boundary condition

% Closure coefficients k-omega model (disregard if turb=0)
alphaW = 13/25; beta_0 = 0.0708; betaW = beta_0; beta_star = 9/100;
sigmaW = 0.5; sigma_star = 3/5; sigma_do = 1/8;
if ikbc == 0
  C_lim = 0; % Turns off stress limiter for k=0 boundary condition
  coef1 = 50; % Coefficient in S_r for k_Np > 5 (Default = 100; Fuhrman et al. 2010 suggest 80) 
elseif ikbc == 1
  C_lim = 7/8; % Default = 7/8
  coef1 = 180; % As suggested by Fuhrman et al. (2010)
end
coef2 = 200; % Coefficient in S_r independent of k_Np (Default = 200)

% Closure coefficients for transitional k-w model (Wilcox 2006, p.206; disregard unless turb=2)
%R_beta = 8; R_k = 6; R_omega = 2.61; % Transitional Reynolds numbers (default values)
R_beta = 8; R_k = 3; R_omega = 2.61; % Re-calibrated values 
alpha0 = 1/9; alpha_star0 = 1/3*betaW; beta_star0 = 9/100; % Closure coefficients

% Filter coefficients, 0-0.25, 0: Off, 0.25: Kills Nyquist
filt_U = 0; % Filter on du/dt
filt_K = 0; % Filter on dk/dt
filt_W = 0; % Filter on d(omega)/dt
filt_C = 0; % Filter on dc/dt

% Move grid point nearest reference level to y=b
ib = find(abs(y-b) == min(abs(y-b))); % Index nearest y=b
if susp == 1; % Only do if a suspended sediment calculation is made
    y(ib) = b; % Replace nearest grid point value with reference level y=b
end
yc = y(ib:end); % y vector for c equation
nyc = length(yc); % Number of vertical points for c equation

% Initial conditions
C(1:nyc,1) = 0; % Suspended sediment concentration profile
U(1:ny,1) = 1e-6; % Velocity profile (m/s)
%U(1:ny,1) = -U_1m*exp(-y./delta1).*sin(-y./delta1); % Laminar wave boundary layer initial conditions
%W(1:ny,1) = 1e-8; % Omega (1/s)
%K(1:ny,1) = 0.00001*1e4.*W.*nu; % Turbulent kinetic energy density (m^2/s^2)
K(1:ny,1) = 0.0006.*U_1m.^2; 
W(1:ny,1) = 1e-6;1e8*U_1m.^2/nu; 

% Adjust bottom boundary conditions
U(1,1) = 0; % U boundary condition
if ikbc == 0 % k=0 wall boundary condition
  K(1,1) = 0; % K boundary condition
elseif ikbc == 1 % dk/dy=0 boundary condition
  K(1,1) = K(2,1);
end

% Settings for time integration (using ode15s)
RelTol = 1e-4; % Relative error tolerance
AbsTol = RelTol*1e-3;
MaxOrder = 5; % Maximum order for ode15s
MaxStep = 0.5; % Maximum time step
ibdf = 'off'; % Use backward differentiation: on | {off}
iNormControl = 'off'; % Use norm control:  on | {off}
istats = 'on'; % Display computational stats after simulation: on | {off}

% Output control during simulation
noutput = 200; % Displays time after every noutput RHS evaluations
ioutplot = 1; % Show plots, 0: No, 1: Yes

% Execute simulation script
main_MatRANS; % Location must be in Matlab path
