function [yh_output,uUf_output,kUf_output] = Assignment1(on_off_plot)
nu = 1e-6; %[m2/s]
rho = 1e3; %[Kg/m3]
mu = nu*rho;

depth = 0.07; %[m]
width = 0.3; %[m]
Ufree = 0.3; %[m/s]

load('..\Ex_1_channel_flow\Exercise1.mat')
for i=1:23
    u = Channel(i).u;
    v = Channel(i).v;
    Dt = Channel(i).tt;
    y(i) = Channel(i).y;
    
    u_mean(i) = sum( u .* Dt ) / sum(Dt);
    u_rmsq(i) = sqrt( sum( (u-u_mean(i) ).^2 .* Dt)/sum(Dt));
    v_mean(i) = sum( v .* Dt ) / sum(Dt);
    v_rmsq(i) = sqrt( sum( (v-v_mean(i) ).^2 .* Dt)/sum(Dt));
    uv(i) = sum( (v-v_mean(i) ).*(u-u_mean(i)) .* Dt)/sum(Dt);
end

u_ext = u_mean; %U extended
y_ext = y;      %y extended
u_ext(end+1) = 0;
y_ext(end+1) = 0;
u_ext(end+1) =  Ufree;
y_ext(end+1) = depth;
y_ext = sort(y_ext);
u_ext = sort(u_ext);

if on_off_plot == 1
    figure('Name','item 1');
    plot(y_ext,u_ext,'k-o');
    xlabel('y [m]','Interpreter','latex');
    ylabel('$u$ [m/s]','Interpreter','latex');
    grid on
end
%% item 2 - depth average velocity
V=1/depth *trapz(y_ext,u_ext);
%% item 3 - Bed friction velocity
Re = depth*width/(2*depth+width) * V/nu;
f = 0.0557/Re^0.25;
Uf = sqrt(f/2)*V; %Darcy-Weisbach
fprintf('From Darcy weissbach  Uf = %.4f m/s\n',Uf);
%% item 4
y_lower_log = 0.003;
y_upper_log = depth*0.3;
sel = y>y_lower_log & y_upper_log > y; % I select the logaritmic zone
myfit = fittype('a + b*log(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b'});
fit_var = fit(y(sel)',u_mean(sel)',myfit);
func =@(x) fit_var.a+ fit_var.b * log(x);

Uf_plot = fit_var.b/2.5; %taken from formula mean(u)=B.Uf ln(y)+A

yh_output=y./depth;
uUf_output = u_mean ./Uf;


if on_off_plot == 1
    figure('Name','item 4')
    semilogx(y_ext,func(y_ext),'k--','DisplayName','fit')
    hold on
    semilogx(y,u_mean,'ko','DisplayName','$\bar{u}$');
    xlabel('y [m]','Interpreter','latex');
    ylabel('$\overline{u}$ [m/s]','Interpreter','latex');
    legend show
    legend('Interpreter','latex','Location','SouthEast');
    grid on
    xline(y_lower_log,'HandleVisibility','off');
    xline(y_upper_log,'HandleVisibility','off');
    text(y_lower_log*1.1,0.1,{'log layer'; sprintf('$ %.0f < y^{+} < %.0f $',...
        y_lower_log*Uf_plot/nu,y_upper_log*Uf_plot/nu) },'interpreter','latex');
    fprintf('From the figure Uf = %.4f m/s \n',Uf_plot);
end
yPlus_lower_log = y_lower_log * Uf_plot/nu;
yPlus_upper_log = y_upper_log * Uf_plot/nu;

fprintf('Y^{+} of the lower boundary of the log-layer is %.4f \n',...
    yPlus_lower_log);


fprintf('Lower bound of the log-Layer is:\n y= %.4f \n y^+ = %.4f \n\n',...
        y_lower_log,yPlus_lower_log);
%% Item 5 - Item 6

uplus = u_mean./Uf_plot;
yplus = y.*Uf_plot/nu;

%Van Driest velocity distribution
k=0.4;
Ad=25;
ydplus = linspace(0,depth*Uf/nu,1000);
VD = 2 * ...
 cumtrapz(ydplus, 1./(1+sqrt(1+4*k^2*ydplus.^2.*(1-exp(-ydplus./Ad)).^2)));

if on_off_plot==1
    figure('Name','Item 5 and 6');
    semilogx(yplus,uplus,'ko','DisplayName','data');
    xlabel('$y^{+}$','interpreter','latex')
    ylabel('$\bar{u}/U_{f}$','interpreter','latex')
    text(1.1,18,'Viscous layer')
    xline(5,'HandleVisibility','off');
    text(6,18,'Buffer layer')
    xline(yPlus_lower_log,'HandleVisibility','off');
    text(yPlus_lower_log*1.1,18,'log layer')
    xline(yPlus_upper_log,'HandleVisibility','off');
    text(yPlus_upper_log*1.1,18,'\tau \neq cte')
    hold on 
    %Van Driest velocity distribution PLOT
    semilogx(ydplus,VD,'k','DisplayName','Van Driest')
    legend('Interpreter','latex','Location','SouthEast');
    legend show
    grid on
end
%% Item 7

% Reynolds stress
tau = rho.*Uf^2*(1-y_ext./depth);
Re_stress = tau(1:end-1) - mu*diff(u_ext)./diff(y_ext);
% turbulent kinetic energy
k = 0.5*(u_rmsq.^2+v_rmsq.^2 .*(1 + 1.8)); %+1.8 for the missing w component

kUf_output = k./Uf^2;

if on_off_plot ==1
    figure('Name','Item 7');
    sgtitle('Turbulence data')
    subplot(1,3,1);
    semilogx(yplus,u_rmsq./Uf,'ks-','DisplayName',...
        "$\sqrt{\overline{u'^{2}}}/U_{f}$");
    hold on
    semilogx(yplus,v_rmsq./Uf,'kd-','DisplayName',...
        "$\sqrt{\overline{v'^{2}}}/U_{f}$");
    semilogx(yplus,-uv./Uf^2,'ko-','DisplayName',...
        "$-\overline{u'v'}/U_{f}^{2}$");

    [hleg1, hobj1] = legend('Interpreter','latex','Location','NorthEast');
    textobj = findobj(hobj1, 'type', 'text');
    set(textobj, 'Interpreter', 'latex', 'fontsize', 9);
    xlabel('$y^{+}$','Interpreter', 'latex', 'fontsize', 12);
    grid on
    text(1.1,1.3,'Viscous layer')
    xline(5,'HandleVisibility','off');
    text(6,1.3,'Buffer layer')
    xline(yPlus_lower_log,'HandleVisibility','off');
    text(yPlus_lower_log*1.1,1.3,'log layer')
    legend show

    subplot(1,3,2);
    semilogx(yplus./(Uf/nu)/depth,u_rmsq./Uf,'ks-','DisplayName',...
        "$\sqrt{\overline{u'^{2}}}/U_{f}$");
    hold on
    semilogx(yplus./(Uf/nu)/depth,v_rmsq./Uf,'kd-','DisplayName',...
        "$\sqrt{\overline{v'^{2}}}/U_{f}$");
    semilogx(yplus./(Uf/nu)/depth,-uv./Uf^2,'ko-','DisplayName',...
        "$-\overline{u'v'}/U_{f}^{2}$");

    %Re stress plot
    semilogx(y_ext(1:end-1)./depth , Re_stress./rho/Uf^2,'k--','DisplayName',...
        "$\overline{u'v'}/U_{f}^{2}$ by the Re tensor");

    text(5*nu/Uf/depth-0.004,1.3,'Viscous layer')
    xline(5*nu/Uf/depth,'HandleVisibility','off');
    text(5*nu/Uf/depth *1.1,1.3,'Buffer layer')
    xline(y_lower_log/depth,'HandleVisibility','off');
    text(y_lower_log/depth*1.1,1.3,'log layer')
    [hleg1, hobj1] = legend('Interpreter','latex','Location','NorthEast');
    textobj = findobj(hobj1, 'type', 'text');
    set(textobj, 'Interpreter', 'latex', 'fontsize', 9);
    xlabel('y/h','Interpreter', 'latex', 'fontsize', 12);
    grid on
    legend show

    %k subplot
    subplot(1,3,3);

    semilogx(yplus./(Uf/nu)/depth , k./Uf^2,'k')
    xlabel('y/h','Interpreter', 'latex', 'fontsize', 12);
    ylabel('$K_{t} / U_{f}$ [J/kg]','Interpreter', 'latex', 'fontsize', 14);
    grid on
    text(5*nu/Uf/depth-0.004,2.5,'Viscous layer')
    xline(5*nu/Uf/depth,'HandleVisibility','off');
    text(5*nu/Uf/depth *1.1,2.5,'Buffer layer')
    xline(y_lower_log/depth,'HandleVisibility','off');
    text(y_lower_log/depth*1.1,2.5,'log layer')
end
%% Item 11

turb_prod = - rho.* uv .* diff(u_ext(1:end-1))./diff(y_ext(1:end-1));

if on_off_plot==1
    figure
    semilogx(y_ext(1:end-2)./depth, turb_prod,'k');
    xlabel('$y/h$', 'interpreter','latex', 'fontsize', 12);
    ylabel("$-\rho.\overline{u'.v'} ~ \partial \overline{u}/\partial y$",...
        'interpreter','latex', 'fontsize', 12);

    text(5*nu/Uf/depth -0.004,7.5,'Viscous layer')
    xline(5*nu/Uf/depth,'HandleVisibility','off');
    text(5*nu/Uf/depth *1.1,7.5,'Buffer layer')
    xline(y_lower_log/depth,'HandleVisibility','off');
    text(y_lower_log/depth*1.1,7.5,'log layer')
    grid on
end


