function handle = plot_me_(x, y, nplots, st_title, ...
    st_xlabel, st_ylabel, st_linespec)
    
    if nplots < 1
        disp('Error calling plot_me_, not enough number of data')
        handle = - 1;
        return
    elseif nplots == 1
        figure('Name', st_title)
        plot(x,y, st_linespec)
        grid on
        xlabel(st_xlabel, 'Interpreter','latex')
        ylabel(st_ylabel, 'Interpreter','latex')
        title(st_title, 'Interpreter','latex')
    else
        figure("Name",st_title)
        plot(x(1,:), y(1,:), st_linespec(1))
        hold on
        for ii = 2:size(x,1)
            plot(x(ii,:), y(ii,:), st_linespec(ii))
        end
        grid on
        xlabel(st_xlabel, 'Interpreter','latex')
        ylabel(st_ylabel, 'Interpreter','latex')
        title(st_title, 'Interpreter','latex')
    end

    handle = 1;

end