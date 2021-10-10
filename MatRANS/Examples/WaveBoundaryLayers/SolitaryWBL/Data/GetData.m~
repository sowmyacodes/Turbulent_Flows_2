% Get data from plot

% FA5010, Figure 9
fname = 'MercatorTimeSeries.jpg'; xvec = [0:1:5460]; yvec = [22:-0.01:-11.7]+0.5;

% Read and plot the image
im = imread(fname);
figure(1), clf
im = im(end:-1:1,:,:);
image(xvec,yvec,im)
axis xy; grid on;

% Set ticks
%set(gca,'xtick',[0:30*60:10000]-593,'ytick',[0:5:20]);
set(gca,'xtick',[0 1207:30*60:4807],'ytick',[0:5:20]); % Match tick marks


% Collect data
[x,y] = ginput; % Click on points, and then hit ENTER to finish

% Plot collected data
hold on; plot(x,y,'r-o'); hold off;

% Then save data as:
% save Mercator.mat x y
