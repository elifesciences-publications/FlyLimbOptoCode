function BodyActivationMeanKinematics(newData, suppressPlots, fontSize, corder)

%% Set some parameters for the analysis

% Window length
winlen = 36;

% Whether to use smoothed and random control trigger
smoothing = false;
useRandomTrigger = false;

% Various analysis parameters
nboot = 1000;
vfMin = 3;
postStimTime = 100;
bootAlpha = 0.01;

% Define desired line widths
dashWidth = .5;
dataWidth = 1;

% Define the y limits for various plots
yLimArray = [-200,200; -.25, 1.75; -7, 3];

%% Set plotting options

% Extract quadrant colorings from main color order
symQuadrantCorder = corder(1:2,:);
quadrantCorder = corder(3:6,:);
segmentCorder = corder(3:8,:);
symSegmentCorder = corder(4:6,:);

% Define the colormaps
ptSz = 100;
[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningPaperColormaps();

%% Set variable lists

% Define the list of variables
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the associated hit variables
hitVarList = {'bodyHitOnset', 'bodyHit','bodyHitX_mm','bodyHitY_mm','bodyEllipseDist_mm', 'videoID'};

% Define all limb related variables for analysis
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');
phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
downVarList = strcat(limbList, '_down_cam');

% Define the body related variables for the analysis
comVarList = {'Orient_Rad_FullRotation', 'xCOM','yCOM'};
hitLocVarListX = strcat(limbList, '_xHitLoc');
hitLocVarListY = strcat(limbList, '_yHitLoc');

% Define the body velocity variables
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
if smoothing
    velVarList = strcat('smooth_', velVarList);
    phaseVarList = strcat('smooth_',phaseVarList);
end

% Get the full variable list for inclusion in the analysis
varList = [velVarList, hitVarList, limbVarListX, limbVarListY,...
    phaseVarList, comVarList, hitLocVarListX, hitLocVarListY, downVarList];


%% Make the trigger

% Extract all the body hits from the dataset
hitData = newData.bodyHitOnset;
trigger = any(hitData, 2);

% Decide if you are using a random trigger
if useRandomTrigger
    trigger = trigger(randperm(length(trigger)));
end

% Define all the time variables
t = (-winlen:winlen)';
t_ms = t./(0.15);

%% Extract needed data

% Get the variables associated with individual time series trials
S = getSlicesSimple(newData, varList, winlen, trigger);
S = cat(3, S{:});

% Partition the data into separate variable types
V = S(:, 1:3, :); % Body velocity data
H = S(:, 4:9, :); % Body hit related variables
L = S(:, 10:21, :); % Limb related variables
Phi = S(:,22:27,:); % Phase related variables
CM = S(:,28:30,:); % Body related position data
hitLocXcam = S(:,31:36,:); % Camera frame limb hit locations x variables
hitLocYcam = S(:,37:42,:); % Camera frame limb hit locations y variables
down = S(:,43:48,:); % Swing/stance variables

% Invert the limb position variables (Changes the coordinate system)
L = -L;

% Convert the rotational velocities from radians to degrees
V(:,1,:) = rad2deg(V(:,1,:));

% Calculate the distances between contralateral limbs
D = hypot(L(:, 1:3, :) - L(:, 4:6,:), L(:, 7:9,:) - L(:, 10:12, :));

% Extract the body hit data
if useRandomTrigger
    hit = rand(size(H,3),size(H,2)) < 0.3;
else
    hit = squeeze(H(t==0,:,:))';
    hitRaw = squeeze(H(t==0,:,:))';
end

%% Group hits by location

% Extract hit locations
hitLocX = squeeze(H(t==0, 3, :));
hitLocY = squeeze(H(t==0, 4, :));

% Group hits by left & right fore & hind quadrants
quadrantIdx = [(hitLocX > 0) & (hitLocY > 0),...
    (hitLocX > 0) & (hitLocY < 0),...
    (hitLocX < 0) & (hitLocY > 0),...
    (hitLocX < 0) & (hitLocY < 0)];

% Group hits by right and left head, thorax, abdomen...
segmentIdx = [
    (hitLocX > 0) & (hitLocY > 1),...
    (hitLocX > 0) & (hitLocY > 0) & (hitLocY < 1),...
    (hitLocX > 0) & (hitLocY < 0),...
    (hitLocX < 0) & (hitLocY > 1),...
    (hitLocX < 0) & (hitLocY > 0) & (hitLocY < 1),...
    (hitLocX < 0) & (hitLocY < 0)];

%% Symmetrize left and right hits

% Define indicators
leftHit = squeeze(any(hit(:,1:3),2));
rightHit = squeeze(any(hit(:,4:6),2));

% Define symmetric version of each of the groups of variables
hitSym = squeeze((hit(:,1:3)>0) | (hit(:,4:6)>0));
Vsym = V; % Body velocity
Vsym(:, [1,3], leftHit) = -Vsym(:, [1,3], leftHit);
Lsym = L; % Limb positions
Lsym(:, :, leftHit) = L(:, [4,5,6,1,2,3,10,11,12,7,8,9], leftHit);
Lsym(:,1:6,leftHit) = -Lsym(:,1:6, leftHit);
Phisym = Phi; % Limb phase
Phisym(:,:,leftHit) = Phisym(:,[4,5,6,1,2,3],leftHit);
downSym = down; % Swing/stance
downSym(:,:,leftHit) = down(:,[4,5,6,1,2,3],leftHit);

%% Compute averages that will be used repeatedly

% Define an indicator for restricting examples based on minimum forward walking speed
vfBefore = squeeze(nanmean(V(t<0,2,:),1));
vfIdx = vfBefore > vfMin;

% Define an indicator for the post stimulus window
postStimWindow = t_ms>0 & t_ms<postStimTime;

% Compute the average post stimulus velocities for each trial
vAfter = squeeze(nanmean(V(postStimWindow,:,:),1))';

%% Kruskal-Wallis one-way analysis of variance

% Compute post stimulus average velocity components
Vsel = squeeze(nanmean(Vsym(postStimWindow,:,:),1))';
Vsel(:,2) = Vsel(:,2) ./ vfBefore; % Fold change for walking speed

% Initialize some variables for storing test statistics
pAnova = nan(3,1);
tAnova = cell(3,1);

% Compute test statistics
for ind = 1:3
    idx = any(hitSym,2) & vfIdx;
    g = sum(hitSym(idx,:) .* [1,2,3],2); % Define a grouping variable
    [pAnova(ind),tAnova{ind}] = kruskalwallis( Vsel(idx,ind), g,'off' );
end

%% Scatter plots of centroid kinematics by activation location

% Plot the forward velocity scatter
MakeFigure;
hold on;
plot([0,0],[-2,2], '--k', 'linewidth', dashWidth);
plot([-2,2],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocX(vfIdx), hitLocY(vfIdx), ptSz, vAfter(vfIdx,2) ./ vfBefore(vfIdx), 'filled');
cbar = colorbar;
caxis([-1,3]);
colormap(cmpBlueRed);
axis('equal');
xlim([-1,1]);
ylim([-2,2]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{||} (fold change)');
ConfAxis(fontSize);

% Plot the yaw velocity scatter
MakeFigure;
hold on;
plot([0,0],[-2,2], '--k', 'linewidth', dashWidth);
plot([-2,2],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocX(vfIdx), hitLocY(vfIdx), ptSz, vAfter(vfIdx,1), 'filled');
cbar = colorbar;
caxis([-500,500]);
colormap(cmpBlueRed);
axis('equal');
xlim([-1,1]);
ylim([-2,2]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{r} (\circ/s)');
ConfAxis(fontSize);

% Plot the lateral velocity scatter
MakeFigure;
hold on;
plot([0,0],[-2,2], '--k', 'linewidth', dashWidth);
plot([-2,2],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocX(vfIdx), hitLocY(vfIdx), ptSz, vAfter(vfIdx,3), 'filled');
cbar = colorbar;
caxis([-15,15]);
colormap(cmpBlueRed);
axis('equal');
xlim([-1,1]);
ylim([-2,2]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{\perp} (mm/s)');
ConfAxis(fontSize);

%% Mean timeseries, split by quadrant

if ~suppressPlots
    
    % Initialize variables for storing values
    meanVel = nan(length(t), 4, 3); % timepoints x quadrants x vel components
    clVel = nan(length(t), 4, 3);
    cuVel = nan(length(t), 4, 3);
    n = nan(4,1);
    
    % Extract out the velocity values
    Vsel = V(:,:,vfIdx);
    Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(vfIdx), [3 2 1]);
    
    % Loop through each hit quadrant
    for ind = 1:4
        % Get the number of trials
        n(ind) = nnz(quadrantIdx(vfIdx,ind));
        % Compute the mean time series for the quadrant
        meanVel(:,ind,:) = nanmean(Vsel(:,:,quadrantIdx(vfIdx,ind)),3);
        
        % Compute confidence interval errorbars
        % NOTE: Default is 95% confidence intervals
        ci = bootci(nboot, {@(x) nanmean(x,1), permute(Vsel(:,:,quadrantIdx(vfIdx,ind)), [3,1,2])});
        clVel(:,ind,:) = ci(1,:,:,:);
        cuVel(:,ind,:) = ci(2,:,:,:);
    end
    
    % Create fome labels
    legendStr = {sprintf('right fore hit (N=%d)', n(1)),sprintf('right hind hit (N=%d)', n(2)),...
        sprintf('left fore hit (N=%d)', n(3)),sprintf('left hind hit (N=%d)', n(4))};
    seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};
    
    % Loop through each of the velocity components and plot
    for ind = 1:3
        MakeFigure;
        % Plot the data
        PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), quadrantCorder);
        
        % Format the plot area
        ylim(yLimArray(ind,:));
        hold on;
        if ind == 1 || ind == 3
            plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        elseif ind==2
            plot(t_ms, ones(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        hold off;
        xlim([min(t_ms), max(t_ms)]);
        xlabel('time (ms)');
        ylabel(seriesLabels{ind});
        legend(legendStr, 'location','northwest');
        ConfAxis(fontSize);
        axis('square');
    end
    
end

%% Mean kinematics, split by quadrant

if ~suppressPlots
    
    % Initialize variables for storing the data
    meanVel = nan(4, 3);
    clVel = nan(4, 3);
    cuVel = nan(4, 3);
    n = nan(4,1);
    
    % Loop through each of the four quadrants
    for ind = 1:4
        
        % Indicator for the current quadrant's trials
        idx = quadrantIdx(:,ind) & vfIdx;
        
        % Number of trials in condition
        n(ind) = nnz(idx);
        
        % Compute the mean post stimulus velocity components per trial
        Vsel = squeeze(mean(V(postStimWindow,:,idx),1));
        
        % Convert forward speed to fold change
        Vsel(2,:) = Vsel(2,:) ./ vfBefore(idx)';
        
        % Compute the group means for each velocity component
        meanVel(ind,:) = nanmean(Vsel,2);
        
        % Compute the error bars
        ci = bootci(nboot, {@(x) nanmean(x,1), permute(Vsel, [2,1])}, 'Alpha', bootAlpha);
        clVel(ind,:) = ci(1,:,:,:);
        cuVel(ind,:) = ci(2,:,:,:);
    end
    
    % Define some legend labels
    legendStr = {'right fore','right hind','left fore', 'left hind'};
    seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};
    
    % Make a separate subplot for each velocity component
    MakeFigure;
    for ind = 1:3
        subplot(1,3,ind);
        % Plot the data
        bar(1:4, meanVel(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', quadrantCorder);
        hold on;
        % Add error bars
        errorbar(1:4, meanVel(:,ind), meanVel(:,ind)-clVel(:,ind), cuVel(:,ind)-meanVel(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0);
        if ind==2
            plot([0 5], [1,1], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        % Format the plot
        title(sprintf('average over %d ms, \\alpha = %0.1e', postStimTime, bootAlpha));
        ylim(yLimArray(ind,:));
        xticks(1:4);
        xticklabels(legendStr);
        xtickangle(20);
        ylabel(seriesLabels{ind});
        ConfAxis(fontSize);
        axis('square');
    end
    
end

%% Mean timeseries, split fore / hind and left / right symmetrized

if ~suppressPlots
    
    % Initialize some variables for storing values
    meanVel = nan(length(t), 2, 3); % Mean values
    clVel = nan(length(t), 2, 3); % Confidence interval lower bound
    cuVel = nan(length(t), 2, 3); % Confidence interval upper bound
    n = nan(4,1); % Number of trials
    
    % Select trajectories
    Vsel = V(:,:,vfIdx);
    
    % Convert forward speed to fold change
    Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(vfIdx), [3 2 1]);
    
    % Loop through fore and hind hit conditions (left-right symmetrized)
    for ind = 1:2
        
        % Get the number of trials
        n(ind) = nnz(quadrantIdx(vfIdx,ind) | quadrantIdx(vfIdx,ind+2));
        
        % Get the velocities (in symmetrized form)
        v = cat(3,Vsel(:,:,quadrantIdx(vfIdx,ind)),[-1,1,-1].*Vsel(:,:,quadrantIdx(vfIdx,ind+2)));
        
        % Compute the mean time series
        meanVel(:,ind,:) = nanmean(v,3);
        
        % Compute the confidence intervals
        ci = bootci(nboot, {@(x) nanmean(x,1), permute(v, [3,1,2])});
        clVel(:,ind,:) = ci(1,:,:,:);
        cuVel(:,ind,:) = ci(2,:,:,:);
    end
    
    % Create labels for plots
    legendStr = {sprintf('fore hit (N=%d)', n(1)),sprintf('hind hit (N=%d)', n(2))};
    seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};
    
    % Loop through each velocity component and plot
    for ind = 1:3
        MakeFigure;
        
        % Plot the data
        PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), symQuadrantCorder);
        
        % Add additional lines indicating prestimulus behavior and stimulus onset
        ylim(yLimArray(ind,:));
        hold on;
        if ind == 1 || ind == 3
            plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        elseif ind==2
            plot(t_ms, ones(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        hold off;
        
        % Format the plot
        xlim([min(t_ms), max(t_ms)]);
        xlabel('time (ms)');
        ylabel(seriesLabels{ind});
        legend(legendStr, 'location','northwest');
        ConfAxis(fontSize);
        axis('square');
    end
    
end

%% Mean kinematics, split fore / hind and left / right symmetrized

if ~suppressPlots
    
    % Initialize some variables for storing values
    meanVel = nan(2, 3);
    clVel = nan(2, 3);
    cuVel = nan(2, 3);
    n = nan(2,1);
    
    % Select the velocity data for included trials
    Vsel = V(:,:,vfIdx);
    
    % Convert forward speed to fold change
    Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(vfIdx), [3 2 1]);
    
    % Loop through fore and hind conditions
    for ind = 1:2
        
        % Count the number of trials
        n(ind) = nnz(quadrantIdx(vfIdx,ind) | quadrantIdx(vfIdx,ind+2));
        
        % Left-right symmetrize the velocity data
        v = cat(3,Vsel(:,:,quadrantIdx(vfIdx,ind)),[-1,1,-1].*Vsel(:,:,quadrantIdx(vfIdx,ind+2)));
        
        % Compute the average post stimulus velocities per trial
        v = squeeze(nanmean(v(postStimWindow,:,:),1));
        
        % Compute the mean post stimulus velocity
        meanVel(ind,:) = nanmean(v,2);
        
        % Compute error bards
        ci = bootci(nboot, {@(x) nanmean(x,1), permute(v, [2,1])}, 'Alpha', bootAlpha);
        clVel(ind,:) = ci(1,:,:,:);
        cuVel(ind,:) = ci(2,:,:,:);
    end
    
    % Create some plot lables
    legendStr = {'fore hit', 'hind hit'};
    seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};
    
    % Loop through each velocity type and plot in subplots
    MakeFigure;
    for ind = 1:3
        subplot(1,3,ind);
        
        % Plot the mean post stimulus velocites
        bar(1:2, meanVel(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', symQuadrantCorder);
        hold on;
        
        % Add error bars
        errorbar(1:2, meanVel(:,ind), meanVel(:,ind)-clVel(:,ind), cuVel(:,ind)-meanVel(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0);
        
        % Format the plot
        if ind==2
            plot([0 3], [1,1], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        title(sprintf('average over %d ms, \\alpha = %0.1e', postStimTime, bootAlpha));
        ylim(yLimArray(ind,:));
        xticks(1:2);
        xticklabels(legendStr);
        %     xtickangle(20);
        ylabel(seriesLabels{ind});
        ConfAxis(fontSize);
        axis('square');
    end
    
end

%% Mean timeseries, split by segment

% Check if this plot should be suppressed
if ~suppressPlots
    
    % Initialize some variables
    meanVel = nan(length(t), 6, 3); % Mean velocity trace
    clVel = nan(length(t), 6, 3); % Confidence interval lower bound
    cuVel = nan(length(t), 6, 3); % Confidence interval upper bound
    n = nan(6,1); % Number of trials
    
    % Select the trials for inclusion
    Vsel = V(:,:,vfIdx);
    
    % Convert forward speed to fold change
    Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(vfIdx), [3 2 1]);
    
    % Loop through each of the segment conditions
    for ind = 1:6
        
        % Get the number of trials
        n(ind) = nnz(segmentIdx(vfIdx,ind));
        
        % Compute mean time series
        meanVel(:,ind,:) = nanmean(Vsel(:,:,segmentIdx(vfIdx,ind)),3);
        
        % Compute the error bars
        ci = bootci(nboot, {@(x) nanmean(x,1), permute(Vsel(:,:,segmentIdx(vfIdx,ind)), [3,1,2])});
        clVel(:,ind,:) = ci(1,:,:,:);
        cuVel(:,ind,:) = ci(2,:,:,:);
    end
    
    % Create some plot labels
    legendStr = {sprintf('right head hit (N=%d)', n(1)),sprintf('right thorax hit (N=%d)', n(2)),sprintf('right abdomen hit (N=%d)', n(3)),...
        sprintf('left head hit (N=%d)', n(4)),sprintf('left thorax hit (N=%d)', n(5)),sprintf('left abdomen hit (N=%d)', n(6))};
    seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};
    
    % Loop through each velocity type and plot
    for ind = 1:3
        MakeFigure;
        % Plot the data
        PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), segmentCorder);
        % Format the plot
        ylim(yLimArray(ind,:));
        hold on;
        if ind == 1 || ind == 3
            plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        elseif ind==2
            plot(t_ms, ones(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        hold off;
        xlim([min(t_ms), max(t_ms)]);
        xlabel('time (ms)');
        ylabel(seriesLabels{ind});
        legend(legendStr, 'location','northwest');
        ConfAxis(fontSize);
        axis('square');
    end
end

%% Mean kinematics, split by segment

% Check if this plot should be suppressed
if ~suppressPlots
    
    % Initialize some variables
    meanVel = nan(6, 3); % Mean data
    clVel = nan(6, 3); % Confidence interval lower bound
    cuVel = nan(6, 3); % Confidence interval upper bound
    n = nan(6,1); % Number of trials
    
    % Loop through each hit type (body segment, left/right)
    for ind = 1:6
        
        % Get the trials associated with the current body region
        idx = segmentIdx(:,ind) & vfIdx;
        
        % Store the number of trials
        n(ind) = nnz(idx);
        
        % Compute the post stimulus average of the velocity components
        Vsel = squeeze(mean(V(postStimWindow,:,idx),1));
        
        % Convert forward speed to fold change
        Vsel(2,:) = Vsel(2,:) ./ vfBefore(idx)';
        
        % Compute the mean across trials
        meanVel(ind,:) = nanmean(Vsel,2);
        
        % Compute the confidence intervals
        ci = bootci(nboot, {@(x) nanmean(x,1), permute(Vsel, [2,1])}, 'Alpha', bootAlpha);
        clVel(ind,:) = ci(1,:,:,:);
        cuVel(ind,:) = ci(2,:,:,:);
    end
    
    % Generate some plot labels
    legendStr = {'right head','right thorax','right abdomen','left head', 'left thorax','left abdomen'};
    seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};
    
    % Plot each of the velocity components in separate subplots
    MakeFigure;
    for ind = 1:3
        subplot(1,3,ind);
        
        % Plot the data
        bar(1:6, meanVel(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', segmentCorder);
        
        % Format the plot
        hold on;
        errorbar(1:6, meanVel(:,ind), meanVel(:,ind)-clVel(:,ind), cuVel(:,ind)-meanVel(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0);
        if ind==2
            plot([0 7], [1,1], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        title(sprintf('average over %d ms, \\alpha = %0.1e', postStimTime, bootAlpha));
        ylim(yLimArray(ind,:));
        xticks(1:6);
        xticklabels(legendStr);
        xtickangle(20);
        ylabel(seriesLabels{ind});
        ConfAxis(fontSize);
        axis('square');
        xlim([0 7]);
    end
end

%% Mean timeseries, split by body segment and left / right symmetrized

% Initialize variables for storing the data
meanVel = nan(length(t), 3, 3); % Mean time series
clVel = nan(length(t), 3, 3); % Confidence interval lower bound
cuVel = nan(length(t), 3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of trials

% Select the velocity data that satisfies the exclusion
Vsel = V(:,:,vfIdx);

% Convert forward velocity to fold change
Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(vfIdx), [3 2 1]);

% Loop through each body segment
for ind = 1:3
    
    % Get the number of trials
    n(ind) = nnz(segmentIdx(vfIdx,ind) | segmentIdx(vfIdx,ind+3));
    
    % Left right symmetrize the data
    v = cat(3,Vsel(:,:,segmentIdx(vfIdx,ind)),[-1,1,-1].*Vsel(:,:,segmentIdx(vfIdx,ind+3)));
    
    % Compute the mean time series
    meanVel(:,ind,:) = nanmean(v,3);
    
    % Compute the error bars
    ci = bootci(nboot, {@(x) nanmean(x,1), permute(v, [3,1,2])});
    clVel(:,ind,:) = ci(1,:,:,:);
    cuVel(:,ind,:) = ci(2,:,:,:);
end

% Generate some plot labels
legendStr = {sprintf('head hit (N=%d)', n(1)),sprintf('thorax hit (N=%d)', n(2)),sprintf('abdomen hit (N=%d)', n(3))};
seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};

% Loop through each velocity component and plot
for ind = 1:3
    MakeFigure;
    
    % Plot the data
    PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), symSegmentCorder);
    
    % Format the plot
    ylim(yLimArray(ind,:));
    hold on;
    if ind == 1 || ind == 3
        plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    elseif ind==2
        plot(t_ms, ones(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    end
    plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    hold off;
    xlim([min(t_ms), max(t_ms)]);
    xlabel('time (ms)');
    ylabel(seriesLabels{ind});
    legend(legendStr, 'location','northwest');
    ConfAxis(fontSize);
    axis('square');
end

%% Mean kinematics, split by segment and left / right symmetrized

% Initialize some variables for storing data
meanVel = nan(3, 3); % Mean data
clVel = nan(3, 3); % Confidence interval lower bound
cuVel = nan(3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of trials

% Select the trials that have above threshold prestimulus forward speed
Vsel = V(:,:,vfIdx);

% Convert forward speed to fold change
Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(vfIdx), [3 2 1]);

% Loop through the body segments
for ind = 1:3
    
    % Count the number of trials
    n(ind) = nnz(segmentIdx(vfIdx,ind) | segmentIdx(vfIdx,ind+3));
    
    % Left-right symmetrize the data
    v = cat(3,Vsel(:,:,segmentIdx(vfIdx,ind)),[-1,1,-1].*Vsel(:,:,segmentIdx(vfIdx,ind+3)));
    % Compute the mean post stimulus response per trial
    v = squeeze(nanmean(v(postStimWindow,:,:),1));
    
    % Compute the mean across trials
    meanVel(ind,:) = nanmean(v,2);
    
    % Compute the error bars
    ci = bootci(nboot, {@(x) nanmean(x,1), permute(v, [2,1])}, 'Alpha', bootAlpha);
    clVel(ind,:) = ci(1,:,:,:);
    cuVel(ind,:) = ci(2,:,:,:);
end

% Create some plot labels
legendStr = {'head hit', 'thorax hit', 'abdomen hit'};
seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};

% Loop through each velocity component and plot in separate subplots
MakeFigure;
for ind = 1:3
    subplot(1,3,ind);
    
    % Plot the data
    bar(1:3, meanVel(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', symSegmentCorder);
    hold on;
    
    % Add error bars
    errorbar(1:3, meanVel(:,ind), meanVel(:,ind)-clVel(:,ind), cuVel(:,ind)-meanVel(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0);
    
    % Format the plot
    if ind==2
        plot([0 4], [1,1], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    end
    title(sprintf('average over %d ms, \\alpha = %0.1e', postStimTime, bootAlpha));
    ylim(yLimArray(ind,:));
    xticks(1:3);
    xticklabels(legendStr);
    xtickangle(20);
    ylabel(seriesLabels{ind});
    ConfAxis(fontSize);
    axis('square');
end

% Define some labels for print statements
hit_names = {'Head Hit', 'Thorax Hit', 'Abdomen Hit'};
vel_names = {'V_r', 'V_par', 'V_perp'};

% Print the CI bounds and whether each is significant
for c = 1:length(vel_names)
    for r = 1:length(hit_names)
        
        % Check if the range intersects with no effect
        if c == 2
            if (clVel(r,c) > 1 && cuVel(r,c) > 1) || (clVel(r,c) < 1 && cuVel(r,c) < 1)
                sig = 'Significant';
            else
                sig = 'Not Significant';
            end
        else
            if (clVel(r,c) > 0 && cuVel(r,c) > 0) || (clVel(r,c) < 0 && cuVel(r,c) < 0)
                sig = 'Significant';
            else
                sig = 'Not Significant';
            end
        end
        
        % Print the result
        fprintf('%s, %s: %.2f to %.2f (%s) \n', hit_names{r}, vel_names{c}, clVel(r,c), cuVel(r,c), sig);
    end
end

%% Timeseries of mean absolute deviation, split by body segment and left / right symmetrized

madLimArray = [0,400;0,20;0,10];

% Initialize variables for storing the data
meanVel = nan(length(t), 3, 3); % Mean time series
clVel = nan(length(t), 3, 3); % Confidence interval lower bound
cuVel = nan(length(t), 3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of trials

% Select the velocity data that satisfies the exclusion
Vsel = V(:,:,vfIdx);

% Loop through each body segment
for ind = 1:3
    
    % Get the number of trials
    n(ind) = nnz(segmentIdx(vfIdx,ind) | segmentIdx(vfIdx,ind+3));
    
    % Left right symmetrize the data
    v = cat(3,Vsel(:,:,segmentIdx(vfIdx,ind)),[-1,1,-1].*Vsel(:,:,segmentIdx(vfIdx,ind+3)));
    
    % Compute the mean time series
    meanVel(:,ind,:) = mad(v,0,3);
    
    % Compute the error bars
    ci = bootci(nboot, {@(x) mad(x,0,1), permute(v, [3,1,2])});
    clVel(:,ind,:) = ci(1,:,:,:);
    cuVel(:,ind,:) = ci(2,:,:,:);
end

% Generate some plot labels
legendStr = {sprintf('head hit (N=%d)', n(1)),sprintf('thorax hit (N=%d)', n(2)),sprintf('abdomen hit (N=%d)', n(3))};
seriesLabels = {'v_{r} (\circ/s)','v_{||} (mm/s)','v_{\perp} (mm/s)'};

% Loop through each velocity component and plot
for ind = 1:3
    MakeFigure;
    
    % Plot the data
    PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), symSegmentCorder);
    
    % Format the plot
    ylim(madLimArray(ind,:));
    hold on;
    plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    hold off;
    xlim([min(t_ms), max(t_ms)]);
    xlabel('time (ms)');
    ylabel([{'mean absolute deviation of '}, seriesLabels{ind}]);
    legend(legendStr, 'location','northwest');
    ConfAxis(fontSize);
    axis('square');
end

%% Bar plots of change in stance duration after activation

end