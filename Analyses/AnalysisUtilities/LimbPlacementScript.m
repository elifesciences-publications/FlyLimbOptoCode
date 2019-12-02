% LimbPlacementScript.m: This script looks at the distribution of next step
% locations following optogenetic manipulation

% Set data paths
SetTurningDataPaths;

% Decide whether you're running for limb or body hits
% hit_case = 'Limb';
hit_case = 'Body';

%% Load data from file and prepare data

% Load the desired data
% NOTE: Select the desired genotype for creating each set of figures
load(bristle_datapath, 'newData');
% load(gr5a_datapath, 'newData');
% load(hs_datapath, 'newData');
% load(ctrl_datapath, 'newData');

% Compute some additional variables required for subsequent analyses
vfMin = 0.5;
[ newData ] = PrepareTurningData(newData, vfMin);

%% Convert the dataset newData into a distribution of next limb placements

% Define the number of bins for the analysis
nBins = 75;

% Get distribution of touchdown x, y positions for all the limbs
switch hit_case
    case 'Limb'
        [ stepDensity, binCenters ] = calcNextStepDensity_LimbHits( newData, nBins );
        hit_strings = {'forelimb hit', 'midlimb hit', 'hindlimb hit'};
    case 'Body'
        [ stepDensity, binCenters ] = calcNextStepDensity_BodyHits( newData, nBins );
        hit_strings = {'head hit', 'thorax hit', 'abdomen hit'};
    otherwise 
        disp('Error: Invalid hit_case.');
end

%% Define some formatting for the plots

% Define some colormaps
[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningPaperColormaps();

% Width of dashed lines
dashWidth = .5;
dataWidth = 1;
numLevels = 20;

%% Plot the joint distributions for the desired dataset

% Plot the raw distributions
for i = 1:3

    % Plot the data
    MakeFigure;
    hold on;
    imagesc(binCenters, binCenters, stepDensity(:,:,i)');
%     contour(binCenters, binCenters, stepDensity(:,:,i)', numLevels, 'linewidth', dataWidth, 'LineColor', 'k');
%     contour(binCenters, binCenters, stepDensity(:,:,i)', numLevels, 'linewidth', dataWidth);
    
    % Format the plot
    hold on;
    plot([0,0],[-3,3], '--k', 'linewidth', dashWidth);
    plot([-3,3],[0,0], '--k', 'linewidth', dashWidth);
    cbar = colorbar;
    caxis([0,1.5]); % NOTE: Check plots before setting
    colormap(cmpRed);
    xlim([-3,3]);
    ylim([-3,3]);
    title(hit_strings{i});
    xlabel('location_{\perp} (mm)');
    ylabel('location_{||} (mm)');
    ylabel(cbar, '1/mm^2');
    ConfAxis('fontSize', 16);
    axis('xy','equal','tight');

end

%% Compute the control for comparison

clear newData;

% Load the control data
load(ctrl_datapath, 'newData');
newData_ctrl = newData;
clear newData;

% Compute some additional variables required for subsequent analyses
vfMin = 0.5;
[ newData_ctrl ] = PrepareTurningData(newData_ctrl, vfMin);

% Get distribution of touchdown x, y positions for all the limbs
switch hit_case
    case 'Limb'
        [ stepDensity_ctrl, binCenters ] = calcNextStepDensity_LimbHits( newData_ctrl, nBins );
    case 'Body'
        [ stepDensity_ctrl, binCenters ] = calcNextStepDensity_BodyHits( newData_ctrl, nBins );
    otherwise 
        disp('Error: Invalid hit_case.');
end

% Compute the difference between 
stepDensity_diff = stepDensity - stepDensity_ctrl;

%% Plot the difference from the control dataset

% Plot the raw distributions
for i = 1:3

    % Plot the data
    MakeFigure;
    hold on;
    imagesc(binCenters, binCenters, stepDensity_diff(:,:,i)');
%     contour(binCenters, binCenters, stepDensity(:,:,i)', numLevels, 'linewidth', dataWidth, 'LineColor', 'k');
%     contour(binCenters, binCenters, stepDensity_diff(:,:,i)', numLevels*2, 'linewidth', dataWidth);
    
    % Format the plot
    hold on;
    plot([0,0],[-3,3], '--k', 'linewidth', dashWidth);
    plot([-3,3],[0,0], '--k', 'linewidth', dashWidth);
    cbar = colorbar;
    caxis([-1.5,1.5]); % NOTE: Check plots before setting
    colormap(cmpBlueRed);
    xlim([-3,3]);
    ylim([-3,3]);
    xlabel('location_{\perp} (mm)');
    ylabel('location_{||} (mm)');
    ylabel(cbar, '1/mm^2');
    title({'Difference From Control',hit_strings{i}});
    ConfAxis('fontSize', 16);
    axis('xy','equal','tight');

end
