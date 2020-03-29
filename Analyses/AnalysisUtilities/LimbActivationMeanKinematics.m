function LimbActivationMeanKinematics(newData, fontSize, centroidCorder, limbCorder, semiMajorAxis, semiMinorAxis, xShift, yShift)

%% Set some parameters for the analysis

% Window length
winlen = 36;

% Whether to use smoothed or raw variables
smoothCentroid = false;
smoothPhase = false;

% Whether to remove instances where multiple limbs are hit
removeMultipleHits = true;

% Various analysis parameters
useRandomTrigger = false;

nboot = 1000;
vfMin = 3;
postStimTime = 100;
bootAlpha = 0.01; % Defines the confidence interval for the plots

% Distance converstion from pixels to mm
mmPerPix = 0.043;

% Define desired line widths
dashWidth = 0.5;
dataWidth = 1;

% Define the y-limits of each of the plots
yLimArray = [-200,200; -.25, 1.75; -7, 3];

if useRandomTrigger
    limbHitsOnly = false;
else
    limbHitsOnly = true;
end

% Define the colormaps
ptSz = 100;
[ cmpBlueRed, cmpRed, cmpBlue ] = MakeTurningPaperColormaps();

%% Extract hit data

[t, t_ms,  hit, hitLocX, hitLocY,V, hitSym, Vsym, Lsym, Phisym, D, downSym]...
    = ExtractLimbHits(newData, winlen, useRandomTrigger, limbHitsOnly,...
    smoothCentroid, smoothPhase, removeMultipleHits, mmPerPix,...
    semiMajorAxis, semiMinorAxis, xShift, yShift);

%% Precompute some averages

% Compute the average trigger forward velocity for each time series
vfBefore = squeeze(nanmean(Vsym(t<0,2,:),1));

% Identify time series where forward velocity is above threshold and there
% is a hit (that satisfies all the above criteria for a hit)
vfIdx = (vfBefore > vfMin) & any(hitSym,2);

% Define the post stimulus window for averaging values
postStimWindow = t_ms>0 & t_ms < postStimTime;

% Compute the post stimulus average velocities
vAfter = squeeze(nanmean(Vsym(postStimWindow,:,:),1))';
vAfterNotSym = squeeze(nanmean(V(postStimWindow,:,:),1))';

% Normalize the poststimulus average forward walking speed by the
% prestimulus average forward walking sped
vAfter(:,2) = vAfter(:,2) ./ vfBefore;
vAfterNotSym(:,2) = vAfterNotSym(:,2) ./ vfBefore;

%% Spatial map of mean kinematic responses

% Create full versions of all the indices and values to treat each limb hit
% separately. Flatten the hitLoc data to match
vfIdxFull = (repmat(vfBefore,6,1) > vfMin) & hit(:);
vAfterFull = repmat(vAfterNotSym,6,1);

hitLocXFlat = hitLocX(:);
hitLocYFlat = hitLocY(:);

% Forward velocity: Plot the locations of each activation colored by forward velocity change
MakeFigure;
hold on;
plot([0,0],[-3,3], '--k', 'linewidth', dashWidth);
plot([-3,3],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocXFlat(vfIdxFull), hitLocYFlat(vfIdxFull), ptSz, vAfterFull(vfIdxFull,2), 'filled');
cbar = colorbar;
caxis([-1,3]);
colormap(cmpBlueRed);
axis('equal');
xlim([-3,3]);
ylim([-3,3]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{||} (fold change)');
ConfAxis(fontSize);

% Yaw velocity: Plot the locations of each activation colored by yaw velocity change
MakeFigure;
hold on;
plot([0,0],[-3,3], '--k', 'linewidth', dashWidth);
plot([-3,3],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocXFlat(vfIdxFull), hitLocYFlat(vfIdxFull), ptSz, vAfterFull(vfIdxFull,1), 'filled');
cbar = colorbar;
caxis([-500,500]);
colormap(cmpBlueRed);
axis('equal');
xlim([-3,3]);
ylim([-3,3]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{r} (\circ/s)');
ConfAxis(fontSize);

% Lateral velocity: Plot the locations of each activation colored by lateral velocity change
MakeFigure;
hold on;
plot([0,0],[-3,3], '--k', 'linewidth', dashWidth);
plot([-3,3],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocXFlat(vfIdxFull), hitLocYFlat(vfIdxFull), ptSz, vAfterFull(vfIdxFull,3), 'filled');
cbar = colorbar;
caxis([-15,15]);
colormap(cmpBlueRed);
axis('equal');
xlim([-3,3]);
ylim([-3,3]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{\perp} (mm/s)');
ConfAxis(fontSize);


%% Distribution of phases of hit limb at activation onset

% Define values for sampling the distribution
xq = (0:0.01:1)';
fPDF = nan(length(xq),3);
n = nan(3,1);

% Compute the pdf for fore, mid, and hindlimb hits
for indH = 1:3
    idx = hitSym(:,indH) & vfIdx; % Checks for limb hit AND above min speed
    n(indH) = nnz(idx); % Number of hits
    fPDF(:,indH) = bqksdensity(squeeze(mod(Phisym(t==0,3+indH,idx),2*pi)), 2*pi*xq); % Compute the pdf
end

legendStr = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};

% Plot Distribution of phases of hit limb at activation onset
MakeFigure;
hold on;
set(gca, 'colororder', centroidCorder);
plot(xq, fPDF, 'linewidth',dataWidth);
legend(legendStr);
axis('square');
xlabel('Phase of Hit Limb at Activation Onset (cycles)');
ylabel('pdf (1/cycles)');
ylim([0 0.3]);
ConfAxis(fontSize);

%% Kruskal-Wallis one-way analysis of variance

% Compute the test statistic
pAnova = nan(3,1);
tAnova = cell(3,1);
for ind = 1:3
    idx = any(hitSym,2) & vfIdx;
    g = sum(hitSym(idx,:) .* [1,2,3],2); % Define a grouping variable for the different hit types
    [pAnova(ind),tAnova{ind}] = kruskalwallis(vAfter(idx,ind), g, 'off');
end

%% Timeseries of average centroid kinematics

% Initialize some empty variables
meanVel = nan(length(t), 3, 3); % Mean velocity
clVel = nan(length(t), 3, 3); % Lower confidence interval
cuVel = nan(length(t), 3, 3); % Upper confidence interval
n = nan(3,1); % Number of examples

% Loop through each of the limb hit types
for ind = 1:3
    
    % Identify the instances where the current limb type was hit
    % NOTE: Trigger selection already deals with case of multiple limb hits
    idx = hitSym(:,ind) & vfIdx;
    
    % Count the number of entries
    n(ind) = nnz(idx);
    
    % Select the velocity data associated with trial for the given hit type
    % NOTE: Vsym has dimensions timepoints, velocity components, trials
    Vsel = Vsym(:,:,idx);
    
    % Compute the forward velocity fold change
    % NOTE: This call to permute gets the sizes to match for element-wise division
    Vsel(:,2,:) = Vsel(:,2,:) ./ permute(vfBefore(idx), [3 2 1]);
    
    % Compute the mean time series for each of the three body velocity components
    meanVel(:,ind,:) = nanmean(Vsel,3);
    
    % Bootstrap to get errorbars
    ci = bootci(nboot, {@(x) nanmean(x,1), permute(Vsel, [3,1,2])});
    clVel(:,ind,:) = ci(1,:,:,:);
    cuVel(:,ind,:) = ci(2,:,:,:);
end

% Define a legend and y-axis labels for the plots below
legendStr = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};
seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'}; %y-axis labels

% Plot the mean time series for each of the body components
for ind = 1:3 % Loop through each body velocity component
    MakeFigure;
    PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), centroidCorder);
    
    % Define the y-limits for the plot
    ylim(yLimArray(ind,:));
    hold on;
    
    % Plot a line to indicate no change (unity for v_{||}, else zero)
    if ind == 2
        plot(t_ms, ones(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    else
        plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    end
    
    % Plot a line to indicate onset of activation
    plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    hold off;
    
    % Format the plot and add labels
    xlim([min(t_ms), max(t_ms)]);
    xlabel('time (ms)');
    ylabel(seriesLabels{ind});
    legend(legendStr);
    ConfAxis(fontSize);
    axis('square');
end

%% Averaged centroid kinematics

% Initialize some variables for storing the plot data
meanVel = nan(3, 3); % Mean velocity
clVel = nan(3, 3); % Confidence interval lower bound
cuVel = nan(3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of examples

% Loop through each limb hit type (fore, mid, hind)
for ind = 1:3
    
    % Select the trials associated with the given limb type AND above threshold walking speed
    idx = hitSym(:,ind) & vfIdx;
    n(ind) = nnz(idx); % Count the number of examples
    
    % Select the mean post stimulus velocities
    Vsel = vAfter(idx,:);
    
    % Compute the mean post-stimulus velocities
    meanVel(ind,:) = nanmean(Vsel,1);
    
    % Bootstrap to generate confidence intervals
    ci = bootci(nboot, {@(x) nanmean(x,1), Vsel}, 'Alpha', bootAlpha);
    clVel(ind,:) = ci(1,:,:,:);
    cuVel(ind,:) = ci(2,:,:,:);
end

% Define a legend and y-axis labels for the plots below
legendStr = {'forelimb hit','midlimb hit', 'hindlimb hit'};
seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'};

% Subplots for each of the velocity components (Post-stimulus mean)
MakeFigure;
for ind = 1:3
    subplot(1,3,ind);
    % Plot a bar for the mean of the velocity component
    bar(1:3, meanVel(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', centroidCorder);
    hold on;
    
    % Plot errorbars based on the confidence intervals
    errorbar(1:3, meanVel(:,ind), meanVel(:,ind)-clVel(:,ind), cuVel(:,ind)-meanVel(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0, 'HandleVisibility','off');
    
    % Plot a line at unity for fold change in forward walking
    if ind==2
        plot([0 4], [1,1], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    end
    
    % Format the plot and add titles
    % NOTE: bootstrapped confidence interval is 100*(1-bootAlpha)
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
hit_names = {'Forelimb Hit', 'Midlimb Hit', 'Hindlimb Hit'};
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

%% Average timeseries of limb separation

% Initialize some variables for storing data
meanLimbDist = nan(length(t), 3, 3); % Mean limb distance
clLimbDist = nan(length(t), 3, 3); % Confidence interval lower bound
cuLimbDist = nan(length(t), 3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of trials

% Loop through each of the limb hit types (fore, mid, hind)
for ind = 1:3
    idx = hitSym(:,ind) & vfIdx;
    n(ind) = nnz(idx);
    
    % Select the corresponding set of contralateral limb distances
    Dsel = D(:,:,idx);
    
    % Compute the mean prestimulus distance between limbs
    DselBefore = nanmean(Dsel(t<0,:,:),1);
    
    % Compute the change in contralateral limb distance
    Dsel = Dsel - DselBefore;
    
    % Compute the time series average
    meanLimbDist(:,:,ind) = nanmean(Dsel,3);
    
    % Bootstrap error bars
    ci = bootci(nboot, {@(x) nanmean(x,1), permute(Dsel, [3,1,2])});
    clLimbDist(:,:,ind) = ci(1,:,:,:);
    cuLimbDist(:,:,ind) = ci(2,:,:,:);
    
end

% Create some labels for the plots
legendStr = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};
seriesLabels = {'forelimb','midlimb','hindlimb'};

% Loop through each hit type and create a separate plot for each
for ind = 1:3
    MakeFigure;
    % Plot the time series data
    PlotAsymmetricErrorPatch(t_ms, meanLimbDist(:,:,ind), clLimbDist(:,:,ind), cuLimbDist(:,:,ind), limbCorder);
    
    % Add lines to indicate mean behavior and stimulus onset
    hold on;
    plot([-250,250],[0 0], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    plot([0 0],[-0.2,0.2], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    hold off;
    
    % Format the plot
    xlim([min(t_ms), max(t_ms)]);
    xlabel('time (ms)');
    ylabel('change in limb separation (mm)');
    legend(seriesLabels);
    title(legendStr{ind});
    ylim([-0.2,0.2]);
    yticks(-0.2:0.05:0.2);
    ConfAxis(fontSize);
    axis('square');
end

%% Average limb separation

% Initialize some variables
meanLimbDist = nan(3, 3); % Mean limb distance
clLimbDist = nan(3, 3); % Confidence interval lower bound
cuLimbDist = nan(3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of trials

% Loop through each limb hit type (fore, mid, hind)
for ind = 1:3
    
    % Create an index for selecting the appropriate limb hit types
    idx = hitSym(:,ind) & vfIdx;
    
    % Count the number of trials
    n(ind) = nnz(idx);
    
    % Select the contralateral limb distances
    Dsel = D(:,:,idx);
    
    % Compute the prestimulus mean limb distance for each example
    DselBefore = nanmean(Dsel(t<0,:,:),1);
    
    % Compute the mean change in contralateral limb distance
    Dsel = nanmean(Dsel(postStimWindow,:,:) - DselBefore, 1);
    
    % Compute the average
    meanLimbDist(:,ind) = nanmean(Dsel,3);
    
    % Generate bootstrapped confidence intervals
    ci = bootci(nboot, {@(x) nanmean(x,1), permute(Dsel, [3,1,2])}, 'Alpha', bootAlpha);
    clLimbDist(:,ind) = ci(1,:,:,:);
    cuLimbDist(:,ind) = ci(2,:,:,:);
    
end

% Create some plot labels
legendStr = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};
seriesLabels = {'forelimb','midlimb','hindlimb'};

% Plot each limb hit type as a subplot
MakeFigure;
for ind = 1:3
    subplot(1,3,ind);
    
    % Generate the bar graphs
    bar(1:3, meanLimbDist(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', limbCorder);
    hold on;
    
    % Add errorbars
    errorbar(1:3, meanLimbDist(:,ind), meanLimbDist(:,ind)-clLimbDist(:,ind), cuLimbDist(:,ind)-meanLimbDist(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0, 'HandleVisibility','off');
    
    % Format the plot
    xticks(1:3);
    xticklabels(seriesLabels);
    ylabel('change in limb separation (mm)');
    title(legendStr{ind});
    ylim([-0.2,0.2]);
    yticks(-0.2:0.05:0.2);
    ConfAxis(fontSize);
    axis('square');
end

%% Timeseries of mean absolute deviation of centroid kinematics

madLimArray = [0,400;0,20;0,10];

% Initialize some empty variables
meanVel = nan(length(t), 3, 3); % Mean velocity
clVel = nan(length(t), 3, 3); % Lower confidence interval
cuVel = nan(length(t), 3, 3); % Upper confidence interval
n = nan(3,1); % Number of examples

% Loop through each of the limb hit types
for ind = 1:3
    
    % Identify the instances where the current limb type was hit
    % NOTE: Trigger selection already deals with case of multiple limb hits
    idx = hitSym(:,ind) & vfIdx;
    
    % Count the number of entries
    n(ind) = nnz(idx);
    
    % Select the velocity data associated with trial for the given hit type
    % NOTE: Vsym has dimensions timepoints, velocity components, trials
    Vsel = Vsym(:,:,idx);
    
    % Compute the mean time series for each of the three body velocity components
    meanVel(:,ind,:) = mad(Vsel,0,3);
    
    % Bootstrap to get errorbars
    ci = bootci(nboot, {@(x) mad(x,0,1), permute(Vsel, [3,1,2])});
    clVel(:,ind,:) = ci(1,:,:,:);
    cuVel(:,ind,:) = ci(2,:,:,:);
end

% Define a legend and y-axis labels for the plots below
legendStr = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};
seriesLabels = {'v_{r} (\circ/s)','v_{||} (mm/s)','v_{\perp} (mm/s)'}; %y-axis labels

% Plot the mean time series for each of the body components
for ind = 1:3 % Loop through each body velocity component
    MakeFigure;
    PlotAsymmetricErrorPatch(t_ms, meanVel(:,:,ind), clVel(:,:,ind), cuVel(:,:,ind), centroidCorder);
    
    % Define the y-limits for the plot
    ylim(madLimArray(ind,:));
    hold on;
    
    % Plot a line to indicate onset of activation
    plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    hold off;
    
    % Format the plot and add labels
    xlim([min(t_ms), max(t_ms)]);
    xlabel('time (ms)');
    ylabel([{'mean absolute deviation of '}, seriesLabels{ind}]);
    legend(legendStr,'location','northwest');
    ConfAxis(fontSize);
    axis('square');
end

%% Swing/stance split bar plots

meanVelHit = nan(3,3,2);
clVelHit = nan(3,3,2);
cuVelHit = nan(3,3,2);
nHit = nan(3,2);

% Iterate over hit limbs
for indH = 1:3
    % Form various indexing vectors
    idxHit = hitSym(:,indH) & vfIdx;
    downHit = squeeze(downSym(t==0,3+indH,:))==1;
    nHit(indH,1) = nnz(idxHit & ~downHit);
    nHit(indH,2) = nnz(idxHit & downHit);
    
    % Compute the mean of each velocity component
    meanVelHit(indH,:,1) = nanmean(vAfter(idxHit & ~downHit,:),1);
    meanVelHit(indH,:,2) = nanmean(vAfter(idxHit & downHit,:),1);
    
    % Bootstrap to generate confidence intervals
    ci = bootci(nboot, {@(x) nanmean(x,1), vAfter(idxHit & ~downHit,:)}, 'Alpha', bootAlpha);
    clVelHit(indH,:,1) = ci(1,:,:,:);
    cuVelHit(indH,:,1) = ci(2,:,:,:);
    
    ci = bootci(nboot, {@(x) nanmean(x,1), vAfter(idxHit & downHit,:)}, 'Alpha', bootAlpha);
    clVelHit(indH,:,2) = ci(1,:,:,:);
    cuVelHit(indH,:,2) = ci(2,:,:,:);
end

seriesLabels = {'v_{r} (\circ/s)','v_{||} (fold change)','v_{\perp} (mm/s)'}; %y-axis labels
legendStr = {'forelimb hit','midlimb hit','hindlimb hit'};

% One figure for each hit type
for ind = 1:3
    MakeFigure;
    
    % One subplot for each of the velocity components (Post-stimulus mean)
    for indV = 1:3
        subplot(1,3,indV);

        hold on;
        bar(1:2, squeeze(meanVelHit(ind,indV,:)), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', centroidCorder(ind,:));
        errorbar(1:2, squeeze(meanVelHit(ind,indV,:)),...
            squeeze(meanVelHit(ind,indV,:)-clVelHit(ind,indV,:)), squeeze(cuVelHit(ind,indV,:)-meanVelHit(ind,indV,:)),...
            ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0, 'HandleVisibility','off');
        
        % Plot a line at y=1 for fold change in forward walking
        if indV==2
            plot([0 3], [1,1], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        end
        
        % Format the plot and add titles
        ylim(yLimArray(indV,:));
        ylabel(seriesLabels{indV});
        title(legendStr{ind});
        xlim([0 3]);
        xticks(1:2);
        xticklabels({sprintf('swing (%d)',nHit(ind,1)),sprintf('stance (%d)',nHit(ind,2))});
        xtickangle(20);
        ConfAxis(fontSize);
        axis('square');
    end
end

end
