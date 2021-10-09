clear; clc; close all;

%% initialize variables

load('Exercise1.mat');

h = 0.07; %m flow depth
u = 0.30; %m/s flow velocity at surface where y = h
nu = 1e-6; %m2/s kinematic viscosity 
rho = 1e3; %kg/m3 density
b = 0.3 %width
%% import data from the channels and run loop to find U mean and rms u

for i = 1:23
    u = Channel(i).u;
    v = Channel(i).v;
    dt = Channel(i).tt;
    y(i) = Channel(i).y;
    h(i) = Channel(i).h;
   
    u_mean(i) = sum(u.*dt)/sum(dt);
    rms_u(i) = sqrt(sum(((u - u_mean(i)).^2).*dt)/sum(dt));
    
    v_mean(i) = sum(v.*dt)/sum(dt);
    rms_v(i) = sqrt(sum(((v - v_mean(i)).^2).*dt)/sum(dt));
    uv_rms(i) = sqrt(sum(abs(( v - v_mean(i) ).*(u - u_mean(i)) .* dt))/sum(dt));
end

%% add in the boundary conditions for u, and y
u_mean_complete = [0, u_mean, 0.3]
v_mean_complete = [0, v_mean, 0]
rms_u_complete = [0.0281, rms_u, 0.0130]
rms_v_complete = [0.0023, rms_v, 0.008]
uv_rms_complete = [0.005, rms_v, 0.008]


h_complete = [0.07,h,0.07]
y_complete = [0,y,0.07]

%% 1 plot u mean against y
figure 
hold on
grid on
plot(y_complete,u_mean_complete)
title('U mean against Y')
xlabel('Y in m')
ylabel('U mean in m/s')
axis([0 0.07 0 0.3])

%% 2 depth avg velocity
%V = (1./h_complete).*trapz(u_mean_complete,h_complete);
V = (1./h_complete).*trapz(y_complete, u_mean_complete);

P = 2.*h + b %perimeter
A = h.*b %area
Rh = A/P %hydraulic radius
RE = Rh.*V./nu
friction = 0.0557./(RE.^0.25) % friction coefficitnt
%% 3 bed friction velocity
U_f = sqrt(friction/2).*V

%% 4 log u mean against y
% % % % % figure
% % % % % hold on
% % % % % grid on
% % % % % semilogx(y_complete,u_mean_complete)
% % % % % title('Logarithmic U mean against Y')
% % % % % xlabel('Y in m')
% % % % % ylabel('U mean in m/s')
% % % % % axis([0 0.05 0 0.3])

sel = y>0.01; % I select the logaritmic zone
myfit = fittype('a + b*log10(x)',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b'});
f = fit(y_complete(sel)',u_mean_complete(sel)',myfit);
func =@(x) f.a+ f.b * log10(x);
figure
semilogx(y_complete,func(y_complete),'DisplayName','fit')
hold on
grid on
semilogx(y_complete,u_mean_complete,'DisplayName','U mean');
title('Logarithmic U mean against Y')
xlabel('y [m]');
ylabel('u [m/s]');
axis([0 0.05 0 0.3])
legend show
fprintf('from the figure Uf = %.4f m/s \n',f.a);

%% 5 u+ against y+
dimensionless_u = u_mean_complete./U_f
dimensionless_y = y_complete.*U_f./nu


semilogx(dimensionless_y,dimensionless_u,'DisplayName','dimensionless U mean');
xlabel('y+');
ylabel('u+');
title('Logarithmic Dimensionless U against Y')
fprintf('from the figure Uf = %.4f m/s \n',f.a);
xlabel('y^{+}')
ylabel('Uf^{+}')
text(1.1,3,'Viscous layer')
xline(6);
text(7,3,'Buffer layer')
xline(47);
text(50,3,'Log layer')
xline(500);
%% 6 Van Driest velocity distribution
k=0.4;
Ad=25;
VD = 2* cumtrapz(dimensionless_y, 1./(1+sqrt(1+4*k^2*dimensionless_y.^2.*(1-exp(-dimensionless_y./Ad)).^2)));

figure
hold on 
semilogx(dimensionless_y,VD,'DisplayName','Vam Driest')
title('van Driest velocity distribution')
xlabel('y+')
ylabel('van Driest velocity')


%% 7 turbulence data in wall units
U_f_square = U_f.^2;

figure


subplot(1,3,1);
loglog(dimensionless_y,rms_u_complete./U_f);
hold on
loglog(dimensionless_y,rms_v_complete./U_f);
loglog(dimensionless_y,(uv_rms_complete./U_f_square));
legend('uprime^2/Uf','vprime^2/Uf','uvprime^2/Uf^2');
xlabel('y+')
title('Turbulence data in wall units')

% 8 turbulence data in outer-flow parameters
y_h = y_complete./h_complete;

subplot(1,3,2);
loglog(y_h,rms_u_complete./U_f);
hold on
loglog(y_h,rms_v_complete./U_f);
loglog(y_h,(uv_rms_complete./U_f_square));
legend('uprime^2/Uf','vprime^2/Uf','uvprime^2/Uf^2');
xlabel('y/h')
title('Turbulence data in outer flow parameters')

% 9 turbulent KE
subplot(1,3,3);
k = 0.5*(rms_u_complete.^2+rms_v_complete.^2 .*(1 + 1.8)); % w'^2 = 1.8v'^2
semilogx(y_complete./(U_f./nu)./h_complete , k./U_f_square)
xlabel('y/h')
title('Turbulent Kinetic Energy per unit mass')

%% 10 Reynolds Stress

tau = rho.*U_f_square.*(1 - y_h);
rms_uv = -tau./rho + nu./rho.*(gradient(u_mean_complete, y_complete));
re_stress = -rms_uv./U_f_square;
figure
plot(y_h, re_stress)
xlabel('y/h')
ylabel('non-dimensional Reynolds Stress')
title('Reynolds stress')


