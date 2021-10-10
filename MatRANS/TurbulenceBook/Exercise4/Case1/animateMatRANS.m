% Example file animating computed MatRANS results.
clear all; close all;

% Input
load('out_MatRANS.mat'); % Load the output workspace
nrepeat = 2; % No. of times to repeat animation
ymax = MatRANS.y(end)/2; % Upper limit to plot for y
ivec = 1:length(MatRANS.t); % Animate all simulated time steps
%ivec = 4*72 + 1 + [0:72]; % Animate only 5th period               

% Animate
figure, clf % Create a new figure
for n = 1:nrepeat % Repeat animation
    for i = ivec % Loop over time
    
        % Plot u(y)
        subplot(1,2,1);
        plot(MatRANS.u(i,:),MatRANS.y,'b-')
        xlabel('u (m/s)'); ylabel('y (m)');
        xlim([-1 1].*1.2*MatRANS.U0m); ylim([0 ymax])
        title(['t = ' num2str(MatRANS.t(i),4) ' s']);
        
        % Plot k(y)
        subplot(1,2,2);
        plot(MatRANS.k(i,:),MatRANS.y,'b-')
        xlabel('k (m^2/s^2)'); ylabel('y (m)');
        xlim([0 1].*0.012.*MatRANS.U0m^2); ylim([0 ymax])
    
        % Update plot
        drawnow
    end
end