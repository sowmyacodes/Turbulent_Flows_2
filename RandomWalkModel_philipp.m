function P = RandomWalkModel_philipp(Tmax,Np,PLOTPATH,PLOTEND)
    %__________________________________________________________________________
    %
    %==============   R A N D O M     W A L K     M O D E L   =================
    %__________________________________________________________________________


    % % % Tmax        = 25;                           % Maximum Time
    % % % Np          = 2500;                         % Number of Particles
    % % % 
    % % % PLOTPATH    = 1;                            % Plot ? -->  0: No  //  1: Yes
    % % % PLOTEND     = 1;                            % Plot ? -->  0: No  //  1: Yes

    % % %     outFileName = ('outRndWalkModel.mat');  % Uncomment last Lines if want to save

    %__________________________________________________________________________
    %
    %==========================================================================
    %__________________________________________________________________________


    %% Inputs =================================================================
    load('out_MatRANS.mat');
    if ~isfield(MatRANS, 'Uf')
        MatRANS.Uf = sqrt(MatRANS.tau0/MatRANS.rho);
    end

    n_t     = MatRANS.n_t;
    u       = MatRANS.u(n_t,:);
    k       = MatRANS.k(n_t,:);
    omega   = MatRANS.omega(n_t,:);
    Uf      = MatRANS.Uf(n_t);                          % Friction Speed
    h       = MatRANS.h_m;
    nu      = MatRANS.nu;
    y       = MatRANS.y;

    % Parameters
    beta    = 9/100;
    vp2     = 1/3 * k;
    vrms    = sqrt(vp2);
    l       = beta^(-1/4) * sqrt(k) ./ omega;
    dt      = l ./ vrms;

    P       = struct('Xp',[],'Yp',[],'Tp',[]);                                  % structure to save particle tracks

    % Boundary Limits
    ybot    = 70 * nu / Uf;
    ytop    = h;

    % ar Random Number Parameters
    RND_mean    = 0;
    RND_std_dev = 1;

    %% Loops ==================================================================
    for jj = 1:Np

        clear yp xp tp dy dx dtp;

        yp(1) = ybot + (ytop-ybot) * rand(1);
        xp(1) = 0;
        tp(1) = 0;

        ii    = 1;

        while 1

            ar          = normrnd(RND_mean,RND_std_dev);

            dtp(ii)     = interp1(y, dt, yp(ii));                               % Time Steps
            tp(ii+1)    = tp(ii) + dtp(ii);                                     % Append particle times

            vrmsp(ii)   = interp1(y, vrms, yp(ii));                             % interpolate particle vrms
            dy(ii)      = ar * vrmsp(ii) * dtp(ii);                             % Change in y


            if (yp(ii) + dy(ii)) > ytop                                         % TOP BOUNDARY REFLECTION
               yp(ii+1) = ytop - (dy(ii) - (ytop - yp(ii)) );

            elseif (yp(ii) + dy(ii)) < ybot                                     % BOT BOUNDARY REFLECTION
               yp(ii+1) = ybot + (-dy(ii) - (yp(ii) - ybot) );

            else
               yp(ii+1) = yp(ii) + dy(ii);                                      % Between Boundaries: Standard case
            end

            u1          = interp1(y, u, yp(ii));                                 % Speed at y
            u2          = interp1(y, u, yp(ii+1));                               % Speed at y + dy

            dx(ii)      = dtp(ii) * (u1+u2)/2;                                  % Change in x

            xp(ii+1)    = xp(ii) + dx(ii);

            if ( tp(ii+1) * Uf / h ) > Tmax
                break;                                                          % Break Statement (wrt. Maximum Non-Dim Time)
            end

            ii = ii+1;                                                          % Increment Time step number
        end

        %% Plot Path ==========================================================
        if PLOTPATH == 1
            figure(1)
                plot(xp(:)./h,yp(:)./h,'.-')
                hold on
        end

        %% Store Results ======================================================
        P(jj).Xp    = xp./h;
        P(jj).Yp    = yp./h;
        P(jj).Tp    = tp*Uf/h;

        P(jj).dXp   = dx./h;
        P(jj).dYp   = dy./h;
        P(jj).dTp   = dtp*Uf/h;
    end

    %% Plot End ===============================================================
    if PLOTEND == 1
        for jj = 1:size(P,2)
            figure(2)
                plot(P(jj).Xp(end), P(jj).Yp(end),'o')
                hold on
        end
    end

% % %     %% Save Output File =======================================================
% % % 
% % %     delete(outFileName);
% % % 
% % %     P.Tmax      = Tmax;
% % %     P.Np        = Np;
% % % 
% % %     save(outFileName,'P');
end