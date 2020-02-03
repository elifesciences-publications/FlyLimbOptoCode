%% OptoMethodPaperFigureScript.m
% A script to regenerate the panels (and more) from Optogenetic Methods
% Paper

% Set data paths
SetTurningDataPaths;

%% Load data from file

% Load the desired data
% NOTE: Select the desired genotype for creating each set of figures
load(bristle_datapath, 'newData');
% load(gr5a_datapath, 'newData');
% load(ctrl_datapath, 'newData');

%% Set options for plotting

% Font size
fontSize = 12;

% Set color scheme
corder = [
    0.105882352941176,0.619607843137255,0.466666666666667;
    0.850980392156863,0.372549019607843,0.007843137254902;
    0.458823529411765,0.439215686274510,0.701960784313725;
    0.905882352941176,0.160784313725490,0.541176470588235;
    0.400000000000000,0.650980392156863,0.117647058823529;
    0.901960784313726,0.670588235294118,0.007843137254902;
    0.650980392156863,0.462745098039216,0.113725490196078;
    0.400000000000000,0.400000000000000,0.400000000000000
    ];

% Set color scheme
centroidCorder = corder(1:3,:);
limbCorder = [
    0.346666666666667,0.536000000000000,0.690666666666667;
    0.915294117647059,0.281568627450980,0.287843137254902;
    0.441568627450980,0.749019607843137,0.432156862745098
    ];

%% Prepare data for analysis 

% Compute some additional variables required for subsequent analyses
vfMin = 0.5;
[ newData ] = PrepareTurningData(newData, vfMin);

%% Generate body hit panels

% Suppress addtional plots
suppressPlots = true;

% Run the analysis function
tic;
BodyActivationMeanKinematics(newData, suppressPlots, fontSize, corder);
fprintf('\nCompleted analysis of body activations in %f seconds.\n', toc);

%% Generate limb hit panels

% Set parameters of more conservative body exclusion region than was used
% in initial hit detection
semiMajorAxis = 1.9; % mm
semiMinorAxis = 0.5; % mm
xShift = 0;
yShift = -0.2; % mm

% Suppress addtional plots
suppressPlots = false;

% Run the analysis function
tic;
LimbActivationMeanKinematics(newData,...
    suppressPlots, fontSize,centroidCorder, limbCorder,...
    semiMajorAxis, semiMinorAxis, xShift, yShift)
fprintf('\nCompleted analysis of limb activations in %f seconds.\n', toc);

%% Generate the comparison plot

% Plot comparison of midlimb hits
OptoComparisonScript;
