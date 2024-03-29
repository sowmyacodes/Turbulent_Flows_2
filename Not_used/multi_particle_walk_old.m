function [P, stats] = multi_particle_walk(MatRANS, Np, beta_star0, ...
    max_Tit, max_T, yplus_bottom, Tobj, random_selection, ...
    plot_trajectories, plot_dispersant_cloud)

    % Function to solve the random_walk of particles in the turbulent flow
    % and perform statistical analysis of the results. As a general
    % notation, the capital letters indicate non-dimensional quantities,
    % whereas the lower case is used for dimensional
    % 
    % Authors:
    %  · Sowmya Srinivasan Iyer
    %  · Videep Goverdhan Kamath
    %  · Yann Birnie Scott 
    %  · Carlos Perez Moreno
    % European Wind Energy Master students · Rotor Design Aerodynamics
    %
    % Denmark Technical University ---- 41129 Turbulent Flows
    % Assignment 2
    % Date: 13/10/2021
    % Version 1 ---- 
    % *********************************************************************
    %%% INPUTS ____________________________________________________________
    %
    % · MatRANS    ---> Output from the RANS flow model
    % · Np         ---> Number of particles
    % · beta_star0 ---> Length scale coefficient
    % · max_Tit    ---> Maximum number of time iterations
    % · max_T      ---> Maximum time for trajectories
    % · random_selection ---> Flag to indicate wether random or linearly 
    %                         spaced initial points are used
    % · plot_trajectories ---> Flag for plotting the trajectories
    % · plot_dispersant_cloud ---> Flag to plot the cloud of points at
    %                              max_T
    % *********************************************************************
    %%% OUTPUTS ___________________________________________________________
    % · P          ---> Stored data of the particles' trajectories
    % · stats      ---> Statistical analysis of the trajectories
    % *********************************************************************
    %%% AUXILIARY _________________________________________________________
    %
    % *********************************************************************

    % If the frictional velocity Uf is not included in the output calculate it
    % Probably unnecesary for our case but nevermind
    if ~isfield(MatRANS, 'Uf')
        MatRANS.Uf = sqrt(MatRANS.tau0/MatRANS.rho);
    end
    
    %% Load the MatRANS Output Data
    
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
    kp_last     = k_last ./ (Uf^2);     % Non-dimensionalise with Uf
    
    % Turbulent dissipation rate
    omega            = MatRANS.omega; 
    omega_last       = omega(end,:);
    % Kinematic viscosity of the fluid
    nu   = MatRANS.nu; 
    
    P = struct('Xp',[],'Yp',[],'Tp',[]); % structure to save particle tracks
    
    y_bottom = yplus_bottom *  nu/ Uf(end); % Dimensional bottom height
    y_top = h;                              % Dimensional top height
    
    yp_bottom = y_bottom/h;                 % Non-dimensional bottom height
    yp_top    = y_top/h;                    % Non-dimensional top height
    
    % Calculate the length scale 
    l = beta_star0^(-0.25) * sqrt(k_last) ./ omega_last;
    
    Lp = l ./ h; % Non-dimensionalise the length scale
    
    rms_Vprime = sqrt(1 / 3 .* kp_last);  % Calculate the v' for each height
    
    % Calculate the time increment vector (the time step size depends on the
    % vertical location)
    Delta_t = Lp ./ rms_Vprime; % Non-dimensional
    
    % Plot the time step size against the vertical dimension
%     plot_me_(Delta_t,y,1, 'Step Size for each vertical position', ...
%              '$\Delta t$', '$y/h$', 'k-o');
    
    %% Preallocate arrays to store the data in each loop iteration
    % Vertical and horizontal positions
    Yp  = zeros(max_Tit,1); % Non-dimensional
    Xp  = zeros(max_Tit,1);
    % Time step size and time values
    Dt  = zeros(max_Tit,1); % Non-dimensional
    Tp  = zeros(max_Tit,1);
    % Streamwise velocity
    Upy = zeros(max_Tit,1); % Non-dimensional
    % Vertical and horizontal displacements
    Dy  = zeros(max_Tit-1,1); % Non-dimensional
    Dx  = zeros(max_Tit-1,1);
    
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
        Y_initial = linspace(yp_bottom, yp_top, Np);
    end
    
    % loop over the number of particles
    for jj = 1 : Np
    
        % Time counter
        ii = 1; 
    
        if random_selection
            % initialize particle start position with a random distribution
            ran_num   = rand(1);    % Generate a random number to get the starting 
                                    % point of the particle
            yp0 = yp_bottom + (yp_top-yp_bottom)*ran_num; % Scale the random number 
                                                          % to the flume 
        else
            yp0 = Y_initial(jj);
        end
                                                      
        % Store the first value
        Yp(1) = yp0;
        % Initialize time variable 
        time = 0;
    
        % loop to calculate one particle track based on the initial position
        while ii <= max_Tit && time <= max_T 
                                    % Limit the number of iterations and 
                                             % the maximum time
            
            % Interpolate to get the time step
            Dt(ii)   = interp1(Y, Delta_t, Yp(ii));
            Tp(ii+1) = time + Dt(ii);
            time     = Tp(ii+1);
    
            % Interpolate to get the non-dimensional streamwise velocity at y
            if ii == 1      % Only for the first time step
                Upy(ii) = interp1(Y, U, Yp(ii));
            else            % Otherwise just assign the previous value
                Upy(ii) = Up_dy;
            end
    
            % Generate random number between -1 and 1 using a normal
            % distribution with mean 0 and stddev 1
            a_r = randn;

            % Interpolate to get the V' rms
            Vpri_rms = interp1(Y,rms_Vprime, Yp(ii));
            % Calculate the vertical displacement and position
            Dy(ii)   = a_r * Vpri_rms * Dt(ii);
            Yp(ii+1) = Yp(ii) + Dy(ii); 
    
            % Check that the vertical position is not out of bounds
            while Yp(ii+1) > y_top/h || Yp(ii+1) < y_bottom/h
                if Yp(ii+1) > y_top/h
                    auxY = Yp(ii+1) - y_top/h;
                    Yp(ii+1) = y_top/h - auxY;
                    Dy(ii) = Yp(ii+1) - Yp(ii);
                elseif Yp(ii+1) < y_bottom/h
                    auxY = y_bottom/h - Yp(ii+1);
                    Yp(ii+1) = y_bottom/h + auxY;
                    Dy(ii) = Yp(ii+1) - Yp(ii);
                end
            end
    
            % Obtain the streamwise velocity at the next vertical position
            Up_dy    = interp1(Y, U, Yp(ii+1));
            % Calculate the streamwise displacement and position
            Dx(ii)   = 0.5 * Dt(ii)*(Up_dy + Upy(ii));
            Xp(ii+1) = Xp(ii) + Dx(ii);
            
            ii = ii + 1; % update time step at the end of the while loop
            
        end
    
        % Plot particle track
        if plot_trajectories
            figure(50)
            plot(Xp(1:ii),Yp(1:ii),'-')
        end
    
        if plot_dispersant_cloud
            figure(51)
            plot(Xp(ii),Yp(ii), 'o')
        end
    
        % Store particle track in variable P using e.g. P(jj).Xp = ...; 
        P(jj).Yp  = Yp(1:ii); 
        P(jj).Xp  = Xp(1:ii);
        P(jj).Tp  = Tp(1:ii);
        P(jj).Dx  = Dx(1:ii);
        P(jj).Dy  = Dy(1:ii);
        P(jj).Dt  = Dt(1:ii);
        P(jj).upy = Upy(1:ii);
    
        
    end
    
    
    %% Analyse results for a series of time targets
    
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

    stats.Xmeans = XP_means;
    stats.Xvars  = XP_vars;
    stats.Ymeans = YP_means;
    stats.Yvars  = YP_vars;
    if plot_trajectories
        figure(50)
        grid on
        xlabel('$X$ [-]','Interpreter','latex')
        ylabel('$Y$ [-]', 'Interpreter','latex')
        title(['Particles path for $N_{P} = $ ', num2str(Np)], ...
            'Interpreter','latex')
        figure_name_path = join('particles_path_',num2str(Np));
        set(gcf, 'PaperPosition', [-2.5 0 30 10]); 
        set(gcf, 'PaperSize', [25 10]); 
        print(figure_name_path, '-dpdf')
    end

    if plot_dispersant_cloud
        figure(51)
        grid on
        xlabel('$X$ [-]','Interpreter','latex')
        ylabel('$Y$ [-]', 'Interpreter','latex')
        title(['Particles dispersant cloud for $N_{P} = $ ', num2str(Np)],...
            'Interpreter','latex')
        figure_name_cloud = join('particles_cloud_',num2str(Np));
        print(figure_name_cloud, '-dpdf')
        set(gcf, 'PaperPosition', [0 0 15 20]); 
        set(gcf, 'PaperSize', [15 20]); 
        print(figure_name_path, '-dpdf')
    end
end
