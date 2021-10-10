%##############################################################################
% Solves the leading order and higher order RANS-equation
% combined with the Wilcox (2006) k-omega turbulence
% model and gradient diffusion equation for suspended
% sediment.  The model uses various finite difference
% approximations for calculating vertical gradients.
% This function is intended to be used by any of Matlab's
% internal ODE solvers, e.g. ode15s.
%
% Release 2 also incorporates:
% 1) Transitional k-omega model (turb=2), see Williams & Fuhrman (2016)
% 2) Sediment mixtures, see Caliskan & Fuhrman (2017)
%
% Programmed by: 
% David R. Fuhrman
% Signe Schloer 
% Johanna Sterner 
% Ugur Caliskan
%##############################################################################

function [out]= ddt_MatRANS(t,f)
    GlobalVars; % Declare global variables

    % Uncomment to ramp up streaming effects at the beginning
    %nramp=8; streaming = tanh(omega.*t./(nramp*2*pi)*pi); % Ramps up over nramp periods
    %nramp=16; iconv = tanh(omega.*t./(nramp*2*pi)*pi); % Ramps up over nramp periods

    % Update RHS evaluation counter
    nrhs = nrhs + 1;

    % The unknowns f=f(u; W; k; C)
    u = real(f(1:ny,1));
    W = real(f(1+ny:2*ny,1));
    K = real(f(1+2*ny:3*ny,1));
    
    % Eliminate any negative values
    K = max(K,0); W = max(W,0); 

    % Reinforce bottom boundary conditions
    if ikbc == 0 % k=0 boundary condition
        K(1) = 0;
    elseif ikbc == 1 % dk/dy=0 boundary condition
        K(1) = K(2);
    elseif ikbc == 2 % Dirichlet condition
        K(1) = kbc;
    end
    u(1) = 0; % No slip condition
    %u(end) = u(end-1); % Forces du/dy=0 at top boundary

    % The pressure gradient term, 1/rho*dp/dx
    tp = t + t0; % Phase adjusted time, to ensure ufs=0 at t=0
    %t=t-500
    if iwave == 0 % No wave signal
        u_fs = 0; dudt_fs = 0;
    elseif iwave == 1 % Second-order Stokes acceleration
        u_fs = U_1m*sin(omega*tp) - U_2m*cos(2*omega*tp); % Desired free stream velocity
        dudt_fs = U_1m*omega*cos(omega*tp) + 2*U_2m*omega*sin(2*omega*tp); % Desired free stream acceleration
    elseif iwave == 2 % Wave form acceleration of Abreau et al. (2010)
        u_fs = U_1m*omega*sqrt(1-r^2)*(sin(omega*tp) + r*sin(phi)/(1+sqrt(1-r^2)));
        dudt_fs = U_1m*omega*sqrt(1-r^2)*(cos(omega*tp)-r*cos(phi)- ...
            r^2/(1+sqrt(1-r^2))*sin(phi)*sin(omega*tp+phi))/(1-r*cos(omega*tp+phi))^2;
    elseif iwave == 3 % Solitary wave acceleratino
        u_fs = U_1m*sech(omega*(t-t0)).^2; % t0=T set in main_MatRANS.m
        dudt_fs = -2*U_1m*omega*sech(omega*(t-t0)).^2*tanh(omega*(t-t0));
    elseif iwave == 4 % N-wave acceleration
        u_fs = U_1m*1.165.*( sech(omega*(t-t0)-3*pi/4).^2 - sech(omega*(t-t0-pi/(2*omega))-3*pi/4).^2 );
        dudt_fs = 2*U_1m*1.165*omega*( ...
                  sech(omega*(t-t0)-3*pi/4).^2*tanh(3*pi/4 - omega*(t-t0)) ...
                - sech(omega*(t-t0-pi/(2*omega))-3*pi/4).^2*tanh(3*pi/4 - omega*(t-t0-pi/(2*omega))) );
         %t, u_fs, dudt_fs
    elseif iwave == 5 % Superimposed sech^2 signals
        u_fs = 0; dudt_fs = 0; % Initialize
        for i = 1:length(Umvec)
          u_fs = u_fs + Umvec(i)*sech( Omegavec(i).*(t - t0 - tnvec(i)) )^2;
          dudt_fs = dudt_fs - 2*Umvec(i)*Omegavec(i)*sech(Omegavec(i)*(t-t0-tnvec(i))).^2*...
                              tanh(Omegavec(i)*(t-t0-tnvec(i)));
        end
        %t, u_fs, dudt_fs
    elseif iwave == 10 % User-defined time series
        u_fs = interp1(t_ud,U_ud,t); % Interpolate velocity values from user-given data
        dudt_fs = interp1(t_ud,Ut_ud,t); % Similarly, interpolate the acceleration
        %t, u_fs, dudt_fs
    end
    dpdx = (streaming*u(end)/c - 1)*dudt_fs; % Pressure gradient term, 1/rho*dp/dx
    dpdx = dpdx + Px; % Add steady contribution
    dpdx = dpdx - iconv*S*u(end)^2/depth; % Add contribution from converging-diverging effects

    % Determine velocity gradient
    dy = diff(y'); % Vector of grid spacings
    dudy = diff(u)./dy; % Velocity gradients
    dudy(ny) = 0; % Top boundary conditions (du/dy=0)

    % Transition modifications (Wilcox 2006, pp. 200--206)
    if turb == 2 % Transition model on
        Re_T = K./(W*nu); % Turbulence Reynolds number
      %if isinf(Re_T(1)); Re_T(1) = 0; end; % Force zero value at the wall
        if ikbc == 0 % Force zero value at the wall
            Re_T(1) = 0; 
        elseif ikbc == 1; % Zero gradient
            Re_T(1) = Re_T(2); 
        elseif ikbc == 2; % Dirichlet condition
            Re_T(1) = kbc./(W(1)*nu); 
            %if isnan(Re_T(1)); Re_T(1) = 0; end
            %Re_T(1) = 1;
        end 
        alpha_star = (alpha_star0 + Re_T/R_k)./(1 + Re_T/R_k); 
        alpha = 13/25.*(alpha0 + Re_T/R_omega)./(1 + Re_T/R_omega)./alpha_star; % Modified closure coefficients
        beta_star = 9/100.*(100*betaW/27 + (Re_T/R_beta).^4)./(1 + (Re_T/R_beta).^4);
        %alpha_star=1; beta_star0=9/100; beta_star=9/100; alpha=13/25; % Uncomment to force std. model values
    else
        alpha_star = 1; beta_star0 = beta_star;
    end
    
    % Calculate omega bottom wall boundary condition
    W_t = max(W, C_lim*abs(dudy)./sqrt(beta_star0./alpha_star)); % Stress-limited Omega tilde
    nu_T = alpha_star.*K./W_t; % Kinematic eddy viscosity
    if ikbc == 0 | ikbc == 2
        tau_b = rho.*(nu)*dudy(1); % Bed shear stress
    elseif ikbc == 1 % dk/dy=0 boundary condition or generalized Dirichlet condition
        tau_b = rho.*(nu_T(1)+nu)*dudy(1); % Bed shear stress
    end
    U_f = sqrt(abs(tau_b)./rho); % Friction velocity
    Shields=(tau_b)./(rho.*g.*(s-1).*d); % Shields parameter
    Shields = ((d./d_50).^hiding_coef).*Shields; % Hiding / exposure effects in mixtures (van Rijn 2007)

    k_Np = max(U_f*k_N/nu,1e-8); % Wall roughness
    if k_Np > 5 % Intermediate or hydraulically rough
        S_r = (coef1/k_Np) + ((coef2/k_Np)^2 - coef1/k_Np)*exp(5-k_Np);%-0.27*10^(-3);
    else % Hydraulically smooth
        S_r = (coef2/k_Np)^2;
    end
    
    W(1,1) = U_f^2*S_r/nu; % Omega at the wall boundary
    W_t = max(W, C_lim*abs(dudy)./sqrt(beta_star0./alpha_star)); % Stress-limited Omega tilde
    nu_T = alpha_star.*K./W_t; % Kinematic eddy viscosity

    % Re-calculate bed shear stress quantities
    if ikbc == 0 | ikbc == 2
        tau_b = rho.*(nu)*dudy(1); % Bed shear stress
    elseif ikbc == 1 % dk/dy=0 boundary condition or generalized Dirichlet condition
        tau_b = rho.*(nu_T(1)+nu)*dudy(1); % Bed shear stress
    end
    
    U_f = sqrt(abs(tau_b)./rho); % Friction velocity
    Shields=(tau_b)./(rho.*g.*(s-1).*d); % Shields parameter
    Shields = ((d./d_50).^hiding_coef).*Shields; % Hiding / exposure effects in mixtures (van Rijn 2007)
    k_Np = max(U_f*k_N/nu,1e-8); % Wall roughness

    % Leading order terms in Navier-Stokes equation
    dudt0(:,1) = -dpdx + nu*gradient(dudy,y') + (turb>0)*gradient(nu_T.*dudy,y');
    %dudt0 = dudt0./(1-streaming.*u./c); % Adds u*du/dx term here
    dudt0(1,1) = 0; % Boundary condition at the seabed (forces u=0 at the wall)

    % Approximate horizontal velocity gradient
    dudx = streaming.*(-dudt0(:,1)/c) + iconv.*(u.*S/depth);

    % Approximate vertical velocity from local continuity equation
    v = -cumtrapz(y,dudx); 
    %v = -cumsimps(y,dudx); % This also works

    % Derivatives in K and Omega equations
    dKdy = [diff(K)./dy; 0]; % Gradients
    if ikbc == 1 % dKdy=0 boundary condition
        dKdy(1) = 0;
    end
    dWdy = [diff(W)./dy; 0];
    sigma_d = sigma_do.*(dKdy.*dWdy >= 0); % Closure coefficient
    dKdy2 = gradient(dKdy,y'); % Second derivatives
    dWdy2 = gradient(dWdy,y');

    % Leading order K equation
    Kw = K./W;
    dKwdy = gradient(Kw,y); % Original
    Two_dudy = dudy.*dudy;
    dKdt0 = (nu_T.*Two_dudy - beta_star.*K.*W + nu*dKdy2 + ...
         gradient(sigma_star.*alpha_star.*Kw.*dKdy,y')).*(turb>0);

    % k wall boundary condition
    if ikbc == 0 | ikbc == 2 % k=0 boundary condition or generalized Dirichlet condition
        dKdt0(1) = 0;
    elseif ikbc == 1 % dk/dy=0 boundary condition
        dKdt0(1) = dKdt0(2);
    end
    dKdx = streaming.*(-dKdt0/c); % Approximate x derivative from leading-order time derivative
    dKdx = dKdx + iconv.*(K.*S/depth); % Also add slope contribution

    % Leading order Omega equation
    dWdt0 = (alphaW*W./W_t.*Two_dudy - betaW.*W.^2 + ...
         sigma_d./W.*dKdy.*dWdy + nu*dWdy2 + ...
         gradient(sigmaW.*alpha_star.*Kw.*dWdy,y')).*(turb>0);
    dWdx = streaming.*(-dWdt0/c); % Approximate x derivative from leading-order time derivative
    dWdx = dWdx + iconv.*(W.*S/depth); % Also add slope contribution

    % Final Navier-Stokes equation
    dudt = dudt0 - (u.*dudx + v.*gradient(u,y')) - 2/3*dKdx*(turb>0); % du/dt
    dudt(1,1) = 0; % Boundary condition at the seabed (keeps u=0)
    dudt(end,1) = dudt(end-1,1); % Forces du/dy=0 at top boundary

    % Final K equation
    dKdt = (-(u.*dKdx + v.*gradient(K,y')) + dKdt0)*(turb>0); % dK/dt
    if ikbc == 0 | ikbc == 2 % k=0 boundary condition or generalized Dirichlet condition
        dKdt(1) = 0; % Also forces k=0 at the bed, with k=0 initial condition
    elseif ikbc == 1
        dKdt(1) = dKdt(2); % dk/dy=0 boundary condition
    end

    % Final Omega equation
    dWdt = (-(u.*dWdx + v.*gradient(W,y')) + dWdt0)*(turb>0); % dW/dt
      

    % The C equation(s) (suspended sediment concentration)
    % Calculate the reference concentration
    if susp == 1 % Do if suspended sediment calculation is desired
        for i=1:length(d) % Loop over individual grain sizes  
            % Determine critical Shields parameter, modified for beach slope
            if Shields(i) >= 0 % Positive
                Shields_c(i) = Shields_cU; % Uphill critical Shields parameter
            else % Negative
                Shields_c(i) = Shields_cD; % Downhill critical Shields parameter
            end

            Shields_rel(i) = abs(Shields(i))-Shields_c(i); Shields_rel(i) = max(Shields_rel(i),0);

            if icb == 1 %  Zyserman & Fredsoe cb
                %cb(i) = 0.331*Shields_rel(i)^1.75/(1+(0.331/0.46)*Shields_rel(i)^1.75);
                cb(i) = w_f(i)*0.331*Shields_rel(i)^1.75/(1+(0.331/0.32)*Shields_rel(i)^1.75); % cb_max replaced with 0.32
            elseif icb == 2 % O'Donoghue & Wright cb
                cb(i) = w_f(i)*0.264;
            elseif icb == 3 || icb == 12 % Engelund & Fredsoe cb and Pickup function 2
                p_s(i) = (1+((pi*mu_d/6)/(Shields_rel(i)+1e-16))^4)^-0.25;
                if abs(Shields(i)) > Shields_c(i)+pi*p_s(i)*mu_d/6;
                    lambda(i) = 0.4*2*((Shields_rel(i)-pi*p_s(i)*mu_d/6)/(0.013*abs(Shields(i))*s))^0.5;
                    cb(i) = w_f(i)*0.6/(1+1/lambda(i))^3;
                else
                    lambda(i) = 0;
                    cb(i) = 0;
                end
            elseif icb == 4  % % Einstein's formula for the reference concentration
                p_s(i) = (1+((pi*mu_d/6)./(Shields_rel(i)+1e-16)).^4).^-0.25;
                cb(i) = w_f(i)*pi*p_s(i)/12;
            elseif icb == 5 % Brors
                cmax = 0.3; theta_crs = 0.25; theta_sheet = 0.75; % Brors values
                cb(i) = w_f(i).*cmax.*(abs(Shields(i))-theta_crs)./(theta_sheet-theta_crs);
                cb(i) = cb(i).*(abs(Shields(i))>theta_crs);
                cb(i) = min(cb(i),cmax);
            elseif icb == 6 % Tanh
                theta_susp = 0.25; cmax = 0.264; theta_sf = 2.0;
                cb(i) = w_f(i).*cmax.*tanh(pi./theta_sf.*(abs(Shields(i))-theta_susp));
                cb(i) = cb(i).*(abs(Shields(i))>theta_susp);
            elseif icb == 11 % Pickup function of van Rijn  
                if abs(Shields(i)) > Shields_c(i)
                    gammaS(i) = 0.00083*(abs(Shields(i))/Shields_c(i)-1)^1.5;
                    pickup(i) = w_s(i)*gammaS(i);
                else
                    pickup(i) = 0;
                end
            end
            cb(i) = cb(i)*(((2*d(i))/b)^(ref_conc_coef));
            if icb == 12 % Pickup function 2   
                if abs(Shields(i)) > Shields_c(i)
                    gammaS(i) = cb(i);
                    pickup(i) = w_s0(i)*gammaS(i);
                else
                    pickup(i) = 0;
                end
            end
        end
        
        C_total = zeros(nyc,1);
        for i=1:length(d)
            C = f(1+3*ny+nyc*(i-1):3*ny+nyc*i);
            C(1) = cb(i);
            C = max(C,0);
            C_total = C_total + C;
        end
    
        for i=1:length(d)
            C = f(1+3*ny+nyc*(i-1):3*ny+nyc*i);
            C(1) = cb(i);
            C = max(C,0);

            if iextrap && icb < 10
                cextrap = C(2) + (C(3)-C(2))/(yc(3)-yc(2))*(y(1)-y(2));
                C(1) = max(cb(i),cextrap);
            end        

            % Sediment diffusivity
            if beta_s == 0 % Use van Rijn (1984) correction
                eps_s = (1 + 2.0.*(w_s0/U_f)^2).*nu_T(end-nyc+1:end); 
            else
                eps_s = beta_s.*nu_T(end-nyc+1:end); % Sediment diffusivity
            end
            eps_s = eps_s + nu; % Add moledular diffusion
            dCdy = [diff(C)./diff(yc'); 0]; % dc/dy
            if icb > 10 % Gradient condition
                dCdy(1) = -pickup/eps_s(1);
            end

            % Settling term
            if Hind_set == 1 % Hindered settling 
                %w_s(1:nyc,i) = w_s0(:,i).*(1-C_total).^nhs(i); % w_s varies with concentration
                w_s = w_s0(i).*(1-C_total).^nhs(i); % w_s varies with concentration
                settling = gradient(w_s.*C,yc');
            else % Constant w_s
                settling = w_s0(i).*dCdy; % w_s*dc/dy
            end

            % Leading order C equation
            diffusion = gradient(eps_s.*dCdy,yc');
            dCdt0 = settling + diffusion;
            dCdx = streaming.*(-dCdt0./c);

            % Final C equation
            dCdt(1+(i-1)*nyc:i*nyc,1) = -(u(end-nyc+1:end).*dCdx + v(end-nyc+1:end).*gradient(C,yc')) + dCdt0; % Full dC/dt
        end
    else
        dCdt(1:length(d)*nyc,1) = zeros(nyc,1); % Otherwise, just set vector to zero
    end

    % Turbulence supression due to density gradients
    if turb && iturbsup && susp
        if iturbsup == 1
            rho_m = rho.*s.*C_total + rho.*(1-C_total); % Density of fluid-sediment mixture
            rhomax = n*rho + (1-n)*rho*s; % Maximum fluid-sediment density
            drdy_b = (rho_m(2)-rho_m(1))./(yc(2)-yc(1)); % d rho_m/dy @ y=b
            rho_m = [rho_m(1) + drdy_b.*(y(1:ib-1)' - y(ib)); rho_m]; % Extrapolate fluid-sediment density
            rho_m = min(rho_m,rhomax); rho_m = max(rho_m,rho); % Bound values
            drho_mdy = gradient(rho_m,y');
            N2 = -g./rho_m.*drho_mdy; % Square of Brunt-Vaisala frequency [1/s^2]
        elseif iturbsup == 2
            dCdy_b = (C_total(2)-C_total(1))./(yc(2)-yc(1)); % dCdy @ y=b
            C_m = [C_total(1) + dCdy_b.*(y(1:ib-1)' - y(ib)); C_total]; % Extrapolate concentration
            C_m = min(C_m,1-n); C_m = max(C_m,0); % Bound values
            dC_mdz = gradient(C_m,y'); % dC/dy
            N2 = -g.*(s-1).*dC_mdz; % Leading-order contribution
        end
        B = N2.*nu_T./sigma_p; % Buoyancy flux [m^2/s^3]
        dKdt = dKdt - B; % Modify k-equation
        c3eps = N2<0; % 1 for N^2<0, 0 for N^2>=0 (Reussink et al. 2009)
        dWdt = dWdt - c3eps.*N2; % Modify omega-equation

        % Enforce K boundary condition
        if ikbc == 0 | ikbc == 2 % k=0 boundary condition or generalized Dirichlet condition
            dKdt(1) = 0; % Also forces k=0 at the bed, with k=0 initial condition
        elseif ikbc == 1
            dKdt(1) = dKdt(2); % dk/dy=0 boundary condition
        end
    end

    % Filter du/dt
    if filt_U > 0
        filt = [filt_U 1-2.*filt_U filt_U]; % Filter
        dudt = filter(filt,1,dudt); % Initial filtering
        dudt = [dudt(2:end); dudt(end)]; % Eliminate phase shift
        dudt(1) = 0; % Enforce boundary condition
        dudt(end) = dudt(end-1); % du/dy=0 at top
    end

    % Filter dK/dt
    if filt_K > 0 && turb
        filt = [filt_K 1-2.*filt_K filt_K]; % Filter
        dKdt = filter(filt,1,dKdt); % Initial filtering
        dKdt = [dKdt(2:end); dKdt(end)]; % Eliminate phase shift
      % Enforce K boundary condition
        if ikbc == 0 | ikbc == 2 % k=0 boundary condition or generalized Dirichlet condition
            dKdt(1) = 0; % Also forces k=0 at the bed, with k=0 initial condition
        elseif ikbc == 1
            dKdt(1) = dKdt(2); % dk/dy=0 boundary condition
        end
    end

    % Filter d(omega)/dt
    if filt_W > 0 && turb
        filt = [filt_W 1-2.*filt_W filt_W]; % Filter
        dWdt = filter(filt,1,dWdt); % Initial filtering
        dWdt = [dWdt(2:end); dWdt(end)]; % Eliminate phase shift
    end

    % Filter dc/dt
    if filt_C > 0 && susp
        filt = [filt_C 1-2.*filt_C filt_C]; % Filter
        dCdt = filter(filt,1,dCdt); % Initial filtering
        dCdt = [dCdt(2:end); dCdt(end)]; % Eliminate phase shift
        dCdt(end) = 0; % Forces C=0 at top boundary condition
    end

    % Create output vector containing all time derivatives
    out = [dudt; dWdt; dKdt; dCdt];

    % Display time to screen every 500th stage evaluation
    if mod(nrhs,noutput) == 0
        disp([ 't = ' num2str(t) ' s, k_N^+ = ' num2str(k_Np) ...
         ', Shields = ' num2str(Shields)...
         ', CPU time = ' num2str(toc) ' s, cb = ' num2str(cb)]) % Display time to screen

        if ioutplot % Also show plots, if desired (to turn off set to 'if 0')
            figure(100); 
            ym = h_m;

            % Plot physical variables
            subplot(3,4,1); plot(u,y); xlabel('u (m/s)'); ylabel('y (m)'); ylim([0 ym]);
            subplot(3,4,2); plot(K,y); xlabel('k (m^2/s^2)'); ylim([0 ym]);
            subplot(3,4,3); plot(W,y); xlabel('\omega (1/s)'); ylim([0 ym]);
            if susp; subplot(3,4,4); plot(C,yc); xlabel('C'); ylim([0 ym]); end

            % Plot time derivatives
            subplot(3,4,5); plot(dudt,y); xlabel('du/dt (m/s^2)'); ylabel('y (m)'); ylim([0 ym]);
            subplot(3,4,6); plot(dKdt,y); xlabel('dk/dt (m^2/s^3)'); ylim([0 ym]);
            subplot(3,4,7); plot(dWdt,y); xlabel('d\omega/dt (1/s^2)'); ylim([0 ym]);
            if susp; subplot(3,4,8); plot(dCdt,yc); xlabel('dC/dt (1/s)'); ylim([0 ym]); end; 

            subplot(3,1,3); hold on; plot(t,u(end),'b.'); hold off;
            %hold on; plot(t,u_fs,'r.'); hold off; % Add analytic value
            xlabel('t (s)'); ylabel('u_0(t) (m/s)'); box on; grid on;

            if turb == 2; % Transitional model
                disp(['max(Re_T) = ' num2str(max(Re_T))]); 
                %figure(200); plot(alpha_star,y,'b-o'); xlim([0 1]); 
                %xlabel('\alpha^*'); ylabel('y (m)');
            end;

            % Update plot
            drawnow; 
        end
    end
