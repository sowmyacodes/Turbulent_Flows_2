% Get data from plot

% Input
%fname = 'Fig6.jpg'; xvec = [-150:0.1:150]; yvec = [-0.8:0.01:1.5].*1.25+0.126; xt = [-150:50:150]; yt = [-2:0.5:2];
%fname = 'Fig12a.jpg'; xvec = [-170.5:0.1:200.5]-10; yvec = [-1.0:0.01:1.8].*1.09-0.15; xt = [-100:100:100]; yt = [-2:1:2];
%fname = 'Fig12b.jpg'; xvec = [-170.5:0.1:200.5]-10; yvec = [-1.0:0.01:1.8].*1.09-0.18; xt = [-100:100:100]; yt = [-2:1:2];
%fname = 'Fig12c.jpg'; xvec = [-170.5:0.1:200.5]-10; yvec = [-1.0:0.01:1.8].*1.09-0.18; xt = [-100:100:100]; yt = [-2:1:2];
%fname = 'Fig12d.jpg'; xvec = [-170.5:0.1:200.5]-10; yvec = [-1.0:0.01:1.8].*1.09-0.21; xt = [-100:100:100]; yt = [-2:1:2];
%fname = 'Fig12e.jpg'; xvec = [-170.5:0.1:200.5]-10; yvec = [-1.0:0.01:1.8].*1.09-0.22; xt = [-100:100:100]; yt = [-2:1:2];
fname = 'Fig12f.jpg'; xvec = [-170.5:0.1:200.5]-10; yvec = [-1.0:0.01:1.8].*1.09-0.22; xt = [-100:100:100]; yt = [-2:1:2];

% Read and plot the image
im = imread(fname);
figure(1), clf
im = im(end:-1:1,:,:);
image(xvec,yvec,im)
axis xy; grid on;

% Set ticks
set(gca,'xtick',xt,'ytick',yt); % Match tick marks


% Collect data
[x,y] = ginput; % Click on points, and then hit ENTER to finish

% Plot collected data
hold on; plot(x,y,'r-o'); hold off;

% Then save data as:
% save Filename.mat x y
