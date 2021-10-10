% Get data from plot

% FA5010, Figure 9
%fname = 'Fig9a.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9b.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9c.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9d.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9e.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9f.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9g.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig9h.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;

% MA5010, Figure 10
%fname = 'Fig10a.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10b.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10c.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10d.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10e.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10f.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10g.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig10h.jpg'; xvec = [-200:0.1:200]./1000; yvec = [-5:0.1:30]./1000;

% Period averaged sediment flux, Figure 12
%fname = 'Fig12a.jpg'; xvec = [-30:0.1:30]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig12b.jpg'; xvec = [-30:0.1:30]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig12c.jpg'; xvec = [-30:0.1:30]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig12d.jpg'; xvec = [-30:0.1:30]./1000; yvec = [-5:0.1:30]./1000;
%fname = 'Fig12e.jpg'; xvec = [-30:0.1:30]./1000; yvec = [-5:0.1:30]./1000;
fname = 'Fig12f.jpg'; xvec = [-30:0.1:30]./1000; yvec = [-5:0.1:30]./1000;


% Read and plot the image
im = imread(fname);
figure(1), clf
im = im(end:-1:1,:,:);
image(xvec,yvec,im)
axis xy; grid on;

%set(gca,'xscale','log'); % Uncomment for SajjadiWaywell1997 data

% Collect data
[x,y] = ginput;

% Plot collected data
hold on; plot(x,y,'r-o'); hold off;

% Then save data as:
% save filename.mat x y