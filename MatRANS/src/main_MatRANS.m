%##############################################################################
% This is the main file for performing hydrodynamic and sediment transport
% calculations in Matlab for current and/or wave boundary layers.  The
% model solves Reynolds-Averaged Navier-Stokes (RANS) equations combined 
% with the k-omega turbulence closure of Wilcox (2006), coupled with a 
% gradient diffusion model for suspended sediment.  The equations are time 
% stepped using Matlab's ode15s solver, which utilizes a dynamic time step 
% control. Both k=0 and dk/dy=0 wall boundary conditions are implemented.
% Numerous options within bed load and suspended load computations are
% available, which should be set up in a separate MatRANS.m input file.
%
% Programmed by: 
% David R. Fuhrman
% Signe Schløer 
% Johanna Sterner 
% Ugur Caliskan
%##############################################################################

% Display to screen
disp('MatRANS: Matlab RANS + k-omega + sediment transport simulation ...'); disp(' '); tic;
pwd % Display current directory to screen

% SIMULATION --------------------------------------------------------------

% Determine additional wave parameters
omega = 2*pi/T; % Angular frequency (1/s)
k = omega/c; % Wavenumber (1/m)
L = 2*pi/k; % Wavelength (m)

% Calculate base sediment fall velocity for every grain size, only if w_s=0
for i=1:length(d)
    if w_s(i) == 0
        w_s0(i) = ( -108*nu/d(i) + sqrt((-108*nu/d(i))^2 + 4*4.2*(4*g*d(i)*(s-1))) ) / (2*4.2); % Explicit solution
        w_s(i) = w_s0(i);
    else
        w_s0(i) = w_s(i);
    end
end

% Hindered settling
if Hind_set == 1
    for i=1:length(d)
        if nhs(i) == 0 % Calculate based on Richardson and Zaki
            Rs(i) = w_s0(i)*d(i)/nu; % Calculate Reynolds number for settling grains
            if Rs(i) < 0.2
                error('Re for grains is to small for expression by Richardson and Zaki')
            elseif Rs(i) > 0.2 && Rs(i) < 1
                nhs(i) = 4.35*Rs(i)^-.03;
            elseif Rs(i) > 1 && Rs(i) < 500
                nhs(i) = 4.45*Rs(i)^-.1;
            else
                nhs(i) = 2.39;
            end
        end
    end
end

% Wave boundary layer variables
delta1 = (2*nu/omega)^(0.5); % Stoke's length (m)
A = abs(U_1m)/omega; % Amplitude of free stream orbital motion (m)
Re = A*abs(U_1m)/nu; % Reynolds number

% Estimate wave friction factor and wave boundary layer thickness, etc.
if iwave
  if (turb == 1) | (turb == 2 && (Re > 1.5e5)) % Turbulent wave boundary layer
    %fw_r = exp(5.5*(A/k_N)^-0.2-6.3); % Rough-turbulent
    fw_r = exp(5.5*(A/k_N)^-0.16-6.7); % Rough-turbulent
    delta_r = k_N*0.09*(A/k_N)^0.82; % Boundary layer thickness (m)
    %fw_s=0.035*Re^-0.16; % Smooth-turbulent
    fw_s=0.04*Re^-0.16; % Smooth-turbulent
    delta_s=A*0.086*Re^-0.11; % (2.56)
    fw = max([fw_r fw_s]); % Wave friction factor (taken as maximum of smooth and rough)
    delta = max([delta_r delta_s]); % Wave boundary layer thickness
  else
    fw=2/sqrt(Re); % Laminar, (5.22)
    delta = 3*pi/4*(2*nu/omega)^0.5; % Laminar boundary layer thickness
  end
Uf = sqrt(fw/2)*abs(U_1m); % Estimated max. friction velocity from wave (m/s)
end

% Determine modified critical Shields parameter due to beach slope
gamma = atan(S); % Longitudinal slope angle
Shields_cU = Shields_c0.*cos(gamma).*(1 + tan(gamma)./tand(phi_s)); % Uphill
Shields_cD = Shields_c0.*cos(gamma).*(1 - tan(gamma)./tand(phi_s)); % Downhill

% Calculate time shift required such that u(t=0)==0
if iwave == 1 % Stokes second order velocity signal
    t0 = fzero(@(tz) U_1m*sin(omega*tz)-U_2m*cos(2*omega*tz),0); % (s)
elseif iwave == 2 % Wave form of Abreau et al. (2010)
    t0 = fzero(@(tz) U_1m*omega*sqrt(1-r^2)*(sin(omega*tz) + r*sin(phi)/(1+sqrt(1-r^2))),0);
elseif iwave == 3 % Solitary wave
    %t0 = T;
end

% Display estimated dimensionless parameters to screen
if abs(U_1m) > 0
    disp(' ');
    disp('Estimated dimensionless wave parameters:')
    disp(['ak = U_1m/c = ' num2str(A*k)])
    disp(['Re = ' num2str(A*abs(U_1m)/nu, '%1.4g')])
    disp(['a/k_N = ' num2str(A/k_N)])
    disp(['dy+ = ' num2str(y(2)*Uf/nu)])
    disp(['dy/k_N = ' num2str(y(2)/k_N)])
    disp(['k_N^+ = ' num2str(k_N*Uf/nu)])
    disp(['fw = ' num2str(fw)])
    disp(['Shields = ' num2str(Uf^2./(g.*(s-1).*d), '%1.4g')])
    disp(['w_s/Uf = ' num2str(w_s/Uf)])
    disp(['w_s/(0.4*Uf) = ' num2str(w_s/(0.4*Uf))])
    disp(['a/h_m = ' num2str(A/h_m)])
    disp(['delta/h_m = ' num2str(delta/h_m)])
end

% Display estimated current parameters to screen
if abs(Px) > 0
    Ufc = sqrt(abs(Px)*h_m); % Estimated friction velocity for pure current
    disp(' ');
    disp('Estimated dimensionless steady current parameters:')
    disp(['dy+ = ' num2str(y(2)*Ufc/nu)])
    disp(['dy/k_N = ' num2str(y(2)/k_N)])    
    disp(['k_N^+ = ' num2str(k_N*Ufc/nu)])
    disp(['Shields = ' num2str(Ufc^2./(g.*(s-1).*d), '%1.4g')])
    disp(['w_s/Uf = ' num2str(w_s/Ufc)])
    disp(['w_s/(0.4*Uf) = ' num2str(w_s/(0.4*Ufc))])
end

% Create vector of initial conditions
nu_T = K./W; % Initial eddy viscosity
initial = [U; W; K; C];

% Add necessary variables in case a single grain size is used
if length(d) == 1
   d_50 = d; % Set d50 to the single grain size
   hiding_coef = 0; % Turn off hiding/exposure coefficient
   w_f = 1; % Set the single weight fraction to unity
   ref_conc_coef = 0; % Set the reference concentration modification coefficient to zero
end

% Perform time integration
nrhs = 0; % Initialize counter
disp(' '); disp('Running ...');
options = odeset('RelTol',RelTol,'AbsTol',AbsTol,'MaxOrder',MaxOrder,'MaxStep',MaxStep,'bdf',ibdf,'NormControl',iNormControl,'stats',istats);
[t,u] = ode15s(@(t,f) ddt_MatRANS(t,f),t_span,initial,options);
disp('Time integration finished.'); disp(' ');

% Extract variables
u = real(u);
U = u(:,1:ny);
W = u(:,1+ny:2*ny);
K = u(:,1+2*ny:3*ny);
C = u(:,1+3*ny:3*ny+nyc*length(d));
C = max(C,0);
C_total = zeros(length(t), nyc);

for i=1:length(d)        
    C_total = C_total + C(:,1+(i-1)*nyc:nyc*i);
end

% POST-PROCESSING -----------------------------------------------------------------------
disp('Post processing variables ...')

% Re-calculate eddy viscosity, shear stress, friction velocity, and Shields parameter
clear nu_T; % Clear eddy viscosity vector
dy = diff(y); % Vector of grid spacings
W_t = zeros(length(t),ny); % Initialize space for new matrices
nu_T = W_t; tau = W_t;
tau_b = zeros(length(t),1);
%U_f = tau_b; Shields = tau_b; Shields_rel = tau_b;
U_f = tau_b;
Shields = zeros(length(t), length(d));
Shields_rel = zeros(length(t), length(d));

% Eliminate any negative values
K = max(K,0); 
W = max(W,0);
C = max(C,0);

for i = 1:length(t)
    dudy = diff(U(i,:))./dy; dudy = [dudy 0]; % Velocity gradient vector (1/s)
    %if ikbc == 0 % k=0 boundary condition
    %  K(i,1) = 0;
    %elseif ikbc == 1 % dk/dy=0 boundary condition
    %  K(i,1) = K(i,2);
    %end

    if turb <= 1
        alpha_star = 1;
        W_t(i,:) = max(W(i,:), C_lim*abs(dudy)/sqrt(beta_star)); % Stress-limited Omega tilde
    elseif turb == 2 % Adjust for transitional turbulence model
        Re_T = K(i,:)./(W(i,:)*nu); % Turbulence Reynolds number
        Re_T(1) = 0; % Dissapears at the wall
        alpha_star = (alpha_star0 + Re_T/R_k)./(1 + Re_T/R_k);
        W_t(i,:) = max(W(i,:), C_lim*abs(dudy)./sqrt(beta_star0./alpha_star));
    end
    nu_T(i,:) = alpha_star.*K(i,:)./(W_t(i,:)+1e-16); % Eddy viscosity (m^2/s)
    tau(i,:) = rho.*(nu_T(i,:) + nu).*dudy; % Shear stress (kg/m/s^2)
    tau_b(i) = tau(i,1); % Bed shear stress (kg/m/s^2)
    U_f(i) = sqrt(abs(tau_b(i))./rho); % Friction velocity (m/s)

    % Shields parameters for every grain size
    for j=1:length(d)
        Shields(i,j) = (tau_b(i))./(rho.*g.*(s-1).*d(j));
        Shields(i,j) = (d(j)/d_50)^hiding_coef.*Shields(i,j);
    end
    
    % Determine omega boundary condition
    k_Np = max(U_f(i)*k_N/nu,1e-8);
    if turb >= 1 % Only do if turbulence model is on
        if k_Np > 5
            S_r = (coef1/k_Np) + ((coef2/k_Np)^2 - coef1/k_Np)*exp(5-k_Np); %-0.27*10^(-3);
        else
            S_r = (coef2/k_Np)^2;
        end
        W(i,1) = U_f(i)^2*S_r/nu;

        % Re-calculate bed shear stress quantities
        if turb == 1
            W_t(i,:) = max(W(i,:), C_lim*abs(dudy)/sqrt(beta_star)); % Stress-limited Omega tilde
	        alpha_star = 1;
        elseif turb == 2
	        W_t(i,1) = max(W(i,1), C_lim*abs(dudy(1))/sqrt(beta_star0/alpha_star(1))); % Stress-limited Omega tilde
        end
        nu_T(i,:) = alpha_star.*K(i,:)./(W_t(i,:)+1e-16); % Eddy viscosity (m^2/s)
        tau(i,:) = rho.*(nu_T(i,:) + nu).*dudy; % Shear stress (kg/m/s^2)
        tau_b(i) = tau(i,1); % Bed shear stress (kg/m/s^2)
        U_f(i) = sqrt(abs(tau_b(i))./rho); % Friction velocity (m/s)

        % Shields parameters for every grain size
        for j=1:length(d)
            Shields(i,j) = (tau_b(i))./(rho.*g.*(s-1).*d(j));
            Shields(i,j) = (d(j)/d_50)^hiding_coef.*Shields(i,j);
        end       
    end

    % Adjust critical Shields parameter
    for j=1:length(d) 
        if Shields(i,j) >= 0 % Positive
            Shields_c(i,j) = Shields_cU; % Uphill critcal Shields parameter
        else % Negative
            Shields_c(i,j) = Shields_cD; % Downhill critical Shields parameter
        end
        Shields_rel(i,j) = abs(Shields(i,j))-Shields_c(i,j); % Relative Shields parameter, (theta - theta_c)
        Shields_rel(i,j) = max(Shields_rel(i,j),0); % Set to zero if negative
    end   
end

% Re-calculate vertical velocity
V = zeros(length(t),ny); % Initialize
if streaming | iconv
    for i = 1:length(t)
        ddt_MatRANS(t(i),[U(i,:) W(i,:) K(i,:) C(i,:)]'); % Evaluate time step
        V(i,:) = v; % Assign vertical velocity
    end
end

% Calculate bed load transport
PhiB_total = 0;
qB_total = 0;
for j=1:length(d)
    if BedLoad == 1 % Meyer-Peter & Muller formula
        PhiB(:,j) = w_f(j).*8.*sign(Shields(:,j)).*Shields_rel(:,j).^1.5; % Non-dimensional
    elseif BedLoad == 2 % Engelund-Fredsoe (eq. (7.59) in Fredsoe & Deigaard 1992)
        PhiB(:,j) = w_f(j).*30/pi/mu_d*sign(Shields(:,j)).*Shields_rel(:,j).*...
	              (sqrt(abs(Shields(:,j))) - 0.7.*sqrt(Shields_c(:,j))); % (7.59)
    elseif BedLoad == 3 % Engelund-Fredsoe (eq (7.58) + (7.54)) 
        p(:,j) = w_f(j).*(1 + (pi/6*mu_d./Shields_rel(:,j)).^4).^(-1/4); % (7.58)
        PhiB(:,j) = 5.*p(:,j).*sign(Shields(:,j)).*(sqrt(abs(Shields(:,j))) - 0.7.*sqrt(Shields_c(j))); % (7.54)
        PhiB(:,j) = PhiB(:,j).*(abs(Shields(:,j))>=Shields_c(j));
    elseif BedLoad == 4 % Nielsen (1992, p. 112)
        PhiB(:,j) = w_f(j).*12.*sign(Shields(:,j)).*Shields_rel(:,j).*sqrt(abs(Shields(:,j))); % (2.3.12)
    end
    qB(:,j) = PhiB(:,j).*sqrt((s-1)*g*d(j)^3); % Dimensional (m^2/s)
    qB_total = qB_total + qB(:,j);
    PhiB_total = PhiB_total + PhiB(:,j);
end

% Calculate suspended sediment transport
if susp == 1
    for i = 1:length(t)      
        ddt_MatRANS(t(i), [U(i,:) W(i,:) K(i,:) C(i,:)]');
        C_total(i,1) = sum(cb);
        for j=1:length(d)
            C(i,1+nyc*(j-1))=cb(j);
            qS(i,j) = trapz(yc,U(i,ib:end).*C(i,1+nyc*(j-1):nyc*j)); % Integrate suspended sediment flux from y=b to h_m (m^2/s)
            PhiS(i,j) = qS(i,j)./sqrt((s-1)*g*d(j)^3); % Dimensionless
        end
        qS_total(i)=sum(qS(i,:));
        PhiS_total(i)=sum(PhiS(i,:));
    end
else
    qS = C; PhiS = C; % Set to zero
    qS_total = C; PhiS_total = 0;
end

% Save data to file
disp('Saving results ...')
save(OutFileName, 't','y','yc','ib','U','V','W','K','C','nu_T',...
    'nu','rho','g','tau','tau_b','U_f','Shields','qB','PhiB','qS','PhiS',...
    'S','w_s','d','k_N','s','U_1m','U_2m', 'C_total', 'qS_total',...
    'qB_total', 'PhiS_total', 'PhiB_total', 'w_f');

% Display elapsed time to screen
disp('Simulation complete.'); toc
disp(' ') % Leave a blank line when finished
