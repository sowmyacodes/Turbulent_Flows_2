function handle = plot_me_(x, y, nplots, st_title, ...
    st_xlabel, st_ylabel, st_linespec)

    % Function to plot results with a standarise style
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
    % · x           ---> Values for the x-axis
    % · y           ---> Values for the y-axis
    % · nplots      ---> Number of plots
    % · st_title    ---> Title of the plot
    % · st_xlabel   ---> X label string
    % · st_ylabel   ---> Y label string
    % · st_linespec ---> Line style specification
    % *********************************************************************
    %%% OUTPUTS ___________________________________________________________
    % · handle      ---> Figure handle 
    % *********************************************************************
    
    if nplots < 1
        disp('Error calling plot_me_, not enough number of data');
        handle = - 1;
        return
    elseif nplots == 1
        handle = figure('Name', st_title);
        plot(x,y, st_linespec);
        grid on
        xlabel(st_xlabel, 'Interpreter','latex');
        ylabel(st_ylabel, 'Interpreter','latex');
        title(st_title, 'Interpreter','latex');
    else
        handle = figure("Name",st_title);
        plot(x(1,:), y(1,:), st_linespec(1));
        hold on
        for ii = 2:size(x,1)
            plot(x(ii,:), y(ii,:), st_linespec(ii));
        end
        grid on
        xlabel(st_xlabel, 'Interpreter','latex');
        ylabel(st_ylabel, 'Interpreter','latex');
        title(st_title, 'Interpreter','latex');
    end

end