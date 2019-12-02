% OptoComparisonScript.m: This script generates a comparison across genotypes
% to visualize the efficacy of the optogenetics method.

%% Load the three datasets and compute time series for limb hits in each
% Bristle
load(bristle_datapath, 'newData');
[bri_meanVel, bri_clVel, bri_cuVel, bri_n] = computeOptoTimeSeries(newData);
clear newData;
% Gr5a 
load(gr5a_datapath, 'newData');
[gr5a_meanVel, gr5a_clVel, gr5a_cuVel, gr5a_n] = computeOptoTimeSeries(newData);
clear newData;
% Control
load(ctrl_datapath, 'newData');
[ctrl_meanVel, ctrl_clVel, ctrl_cuVel, ctrl_n] = computeOptoTimeSeries(newData);
clear newData;

%% Plot the comparison for midlimb v_r under optogenetic manipulation

% NOTE: meanVel has the following structure: timepoints x hit type [fore, mid, hind] x velocity type [V_r, V_par, V_perp]

% Define a legend and y-axis labels for the plots below
legendStr = {'control hit','bristle hit', 'gr5a hit'};
seriesLabels = {'v_{r} (\circ/s)'}; %y-axis labels

% Define the time variable that will be used for plotting
winlen = 36;
t = (-winlen:winlen)';
t_ms = t./(0.15);

% Define the line width of the dashed lines
dashWidth = .5;

% Define the y-limits of each of the plots
yLimArray = [-200,200; -.25, 1.75; -7, 3];

% Define the colors for the plots
corder = [.5,.5,.5; 0,0,1; 1,0,0]; % (gray, blue, red);

% Group the plot data 
meanVel = [ctrl_meanVel(:,2,1), bri_meanVel(:,2,1), gr5a_meanVel(:,2,1)];
clVel = [ctrl_clVel(:,2,1), bri_clVel(:,2,1), gr5a_clVel(:,2,1)];
cuVel = [ctrl_cuVel(:,2,1), bri_cuVel(:,2,1), gr5a_cuVel(:,2,1)];

% Plot the figure
MakeFigure;
hold on;
PlotAsymmetricErrorPatch(t_ms, meanVel, clVel, cuVel, corder);

% Define the y-limits for the plot
ylim(yLimArray(1,:));

% Plot a line to indicate no change
plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');

% Plot a line to indicate onset of activation
plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
hold off;

% Format the plot and add labels
xlabel('time (ms)');
ylabel(seriesLabels);
legend(legendStr);
ConfAxis('fontSize', 16);
title('Midlimb Hits');
axis('square');