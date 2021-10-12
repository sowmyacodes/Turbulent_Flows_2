function [Xmean, Ymean, Xvar, Yvar] = get_particle_stats(P,Tobj)
    % _____________________________________________________________________
    %%% INPUTS
    % P       ---> Structure with the particle trajectories and time
    % Tobj    ---> Time for which the stats are calculated
    % _____________________________________________________________________
    %%% OUTPUTS
    % Xmean   ---> Mean value of the x-position
    % Ymean   ---> Mean value of the y-position
    % Xvar    ---> Variance of the x-position
    % Yvar    ---> Variance of the y-position
    % _____________________________________________________________________
    %%% AUXILIARY
    % Np      ---> Number of particles
    % Xp_vals ---> Auxiliary array for the X values at the Tobj
    % Yp_vals ---> Auxiliary array for the Y values at the Tobj
    % ii      ---> Loop variable

    
    Np = size(P,2); % Number of particles to analyse from the structure

    % Preallocate for auxiliary calculations
    XP_vals = zeros(Np,1);
    YP_vals = zeros(Np,1);
    
    % Interpolate for each particle at the target time
    for ii = 1:Np
    
        XP_vals(ii) = interp1(P(ii).Tp, P(ii).Xp, Tobj);
        YP_vals(ii) = interp1(P(ii).Tp, P(ii).Yp, Tobj);
    
    end
    
    % Calculate the mean
    Xmean = mean(XP_vals);
    Ymean = mean(YP_vals);
    
    % Calculate the variance
    Xvar  = var(XP_vals);
    Yvar  = var(YP_vals);

end