%% Set some parameters for the analysis

% Window length
winlen = 36;

% Whether to use smoothed or raw variables
smoothCentroid = false;
smoothPhase = false;

% Whether to remove instances where multiple limbs are hit
removeMultipleHits = true;

% Suppress addtional plots
supressPlots = true;

% Various analysis parameters
useRandomTrigger = false;
nboot = 1000;
vfMin = 3;
postStimTime = 100;
bootAlpha = 0.01; % Defines the confidence interval for the plots

% Distance converstion from pixels to mm
mmPerPix = 0.043;

% Define desired line widths
dashWidth = .5;
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

%% Set variable lists

% Define the list of limb prefixes
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the list of hit variables
hitVarList = strcat(limbList, '_hit');

% Define lists for the limb position and limb phase related variables
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');
phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
freqVarList = strcat('InstantaneousFrequency_', limbList, 'y');

% Define the body variables of interest
comVarList = {'Orient_Rad_FullRotation', 'xCOM','yCOM'};

% Define the list of limb hit location variables
hitLocVarListX = strcat(limbList, '_xHitLoc');
hitLocVarListY = strcat(limbList, '_yHitLoc');

% Define the the list of body velocity related variables
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
if smoothCentroid
    velVarList = strcat('smooth_', velVarList);
end

% Update the limb phase and frequency variables, if smoothing
if smoothPhase
    phaseVarList = strcat('smooth_',phaseVarList);
    freqVarList = strcat('smooth_', freqVarList);
end

% Append all the variable lists together
varList = [velVarList, hitVarList, limbVarListX, limbVarListY,...
    phaseVarList, comVarList, hitLocVarListX, hitLocVarListY, freqVarList];


%% Make the trigger

% Gather all the hit variables into one location
hitData = newData{:, hitVarList};

% Define the trigger variable for selecting trajectories
trigger = any(hitData, 2);

% Remove cases from the trigger where multiple limbs are hit simultaneously
if removeMultipleHits
    trigger = trigger & (sum(hitData,2) == 1);
end

% If using a random trigger, generate the random trigger instead
if useRandomTrigger
    trigger = trigger(randperm(length(trigger)));
end

% Define the time variable that will be used for plotting
t = (-winlen:winlen)';
t_ms = t./(0.15);

%% Extract needed data

% Convert the trigger variable into an array of time series centered on the trigger
S = getSlicesSimple(newData, varList, winlen, trigger);
S = cat(3, S{:});

% Extract out each of the groups of variables from the slice data structure
V = S(:, 1:3, :); % Body velocity
H = S(:, 4:9, :); % Hit logicals
L = S(:, 10:21, :); % Limb positions in mm
Phi = S(:,22:27,:); % Limb phase variables
CM = S(:,28:30,:); % Body position variables
hitLocXcam = S(:,31:36,:); % Hit x_location (Centers of each activation)
hitLocYcam = S(:,37:42,:); % Hit y_location (Centers of each activation)
Phidot = S(:,43:48,:); % Limb frequency variables

L = -L; % Reverse both limb position axes

V(:,1,:) = rad2deg(V(:,1,:)); % Convert yaw rate to degrees
Phidot = Phidot / (2*pi) * 150; % Convert limb frequency units of cycles per second

% Compute the distance between contralateral limbs
D = hypot(L(:, 1:3, :) - L(:, 4:6,:), L(:, 7:9,:) - L(:, 10:12, :));

if useRandomTrigger
    hit = rand(size(H,3),size(H,2)) < 0.2;
else
    hit = squeeze(H(t==0,:,:))';
end

%% Transform hit locations into the fly frame

% Get the body and hit locations (in camera frame) at activation time
cm0 = squeeze(CM(t==0,:,:))';
hXcam0 = squeeze(hitLocXcam(t==0,:,:))';
hYcam0 = squeeze(hitLocYcam(t==0,:,:))';

% Transform the hit coordinates
% Yes, this seems backwards, but it's the correct transformation between
% right-handed and left-handed coordinates
th = cm0(:,1)-pi/2;
hitLocX = (hXcam0~=0).*(hXcam0-cm0(:,2)).*sin(th)+(hYcam0~=0).*(hYcam0-cm0(:,3)).*cos(th);
hitLocY = (hXcam0~=0).*(hXcam0-cm0(:,2)).*cos(th)-(hYcam0~=0).*(hYcam0-cm0(:,3)).*sin(th);

% Convert hit coodinates from pixels to mm
hitLocX = hitLocX .* mmPerPix;
hitLocY = hitLocY .* mmPerPix;

% Define whether a hit is within the body region
bodyEllipseDist = hypot((hitLocX-xShift)./semiMinorAxis, (hitLocY-yShift)./semiMajorAxis);

% Decide if excluding limb hits within the body ellipse
if limbHitsOnly
    outsideBodyEllipse = bodyEllipseDist > 1;
else
    outsideBodyEllipse = bodyEllipseDist < 1;
end
hit = hit & outsideBodyEllipse;

%% Symmetrize left and right hits

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
Phidotsym = Phidot; % Limb frequency
Phidotsym(:,:,leftHit) = Phidotsym(:,[4,5,6,1,2,3],leftHit);

%% Spatial map of mean kinematic responses

% Compute the average trigger forward velocity for each time series
vfBefore = squeeze(nanmean(V(t<0,2,:),1));

% Identify time series where forward velocity is above threshold and there
% is a hit (that satisfies all the above criteria for a hit)
vfIdx = (vfBefore > vfMin) & any(hit,2);

% Define the post stimulus window for averaging values
postStimWindow = t_ms>0 & t_ms<postStimTime;

% Compute the post stimulus average velocities (For the averaging window)
% Note: This includes values for all three components of the body velocity, 
% unlike vfBefore
vAfter = squeeze(nanmean(V(postStimWindow,:,:),1))';

% Create full versions of all the indices and values to treat each limb hit 
% separately. Flatten the hitLoc data to match
vfIdxFull = (repmat(vfBefore,6,1) > vfMin) & hit(:);
vfBeforeFull = repmat(vfBefore,6,1);
vAfterFull = repmat(vAfter,6,1);

hitLocXFlat = hitLocX(:);
hitLocYFlat = hitLocY(:);

% Forward velocity: Plot the locations of each activation colored by forward velocity change 
MakeFigure;
hold on;
plot([0,0],[-3,3], '--k', 'linewidth', dashWidth);
plot([-3,3],[0,0], '--k', 'linewidth', dashWidth);
scatter(hitLocXFlat(vfIdxFull), hitLocYFlat(vfIdxFull), ptSz, vAfterFull(vfIdxFull,2) ./ vfBeforeFull(vfIdxFull), 'filled');
cbar = colorbar;
caxis([-1,3]);
colormap(cmpBlueRed);
axis('equal');
xlim([-3,3]);
ylim([-3,3]);
xlabel('hit_{\perp} (mm)');
ylabel('hit_{||} (mm)');
ylabel(cbar, 'v_{||} (fold change)');
ConfAxis('fontSize', 16);

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
ConfAxis('fontSize', 16);

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
ConfAxis('fontSize', 16);


%% Distribution of phases of hit limb at activation onset

% Check if this plot should be suppressed
if ~supressPlots
    
    % Define values for sampling the distribution
    xq = (0:0.01:1)';
    fPDF = nan(length(xq),3);
    n = nan(3,1);

    % % NOTE: These are both calculated above and can be deleted
    % % Get the average prestimulus forward velocity for each time series
    % vfBefore = squeeze(nanmean(V(t<0,2,:),1));
    % % Get the indices of time series that have a prestimulus forward velocity above threshold
    % vfIdx = vfBefore > vfMin;

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
    ConfAxis('fontSize', 16);
end

%% Kruskal-Wallis one-way analysis of variance

% Define the post stimulus window
postStimWindow = t_ms>0 & t_ms<postStimTime;

% % NOTE: These are defined above
% vfBefore = squeeze(nanmean(V(t<0,2,:),1));
% vfIdx = vfBefore > vfMin;

% Calculate the average post stimulus symmetric velocities
Vsel = squeeze(nanmean(Vsym(postStimWindow,:,:),1))';

% Calculate the ratio of prestimulus and poststimulus forward velocity
Vsel(:,2) = Vsel(:,2) ./ vfBefore;

% Compute the test statistic
pAnova = nan(3,1);
tAnova = cell(3,1);
for ind = 1:3
    idx = any(hitSym,2) & vfIdx;
    g = sum(hitSym(idx,:) .* [1,2,3],2); % Define a grouping variable for the different hit types
    [pAnova(ind),tAnova{ind}] = kruskalwallis( Vsel(idx,ind), g,'off' );
end

%% Timeseries of average centroid kinematics

% Initialize some empty variables
meanVel = nan(length(t), 3, 3); % Mean velocity
clVel = nan(length(t), 3, 3); % Lower confidence interval
cuVel = nan(length(t), 3, 3); % Upper confidence interval
n = nan(3,1); % Number of examples

% % NOTE: These have been computed above
% vfBefore = squeeze(nanmean(V(t<0,2,:),1));
% vfIdx = vfBefore > vfMin;

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
    
    % Plot a line to indicate no change
    % Rotational and Lateral velocities
    if ind == 1 || ind == 3
        plot(t_ms, zeros(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    % Forward velocity
    elseif ind==2
        plot(t_ms, ones(length(t_ms),1), '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    end
    
    % Plot a line to indicate onset of activation
    plot([0 0], ylim, '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
    hold off;
    
    % Format the plot and add labels
    xlabel('time (ms)');
    ylabel(seriesLabels{ind});
    legend(legendStr);
    ConfAxis('fontSize', 16);
    axis('square');
end

%% Averaged centroid kinematics

% Initialize some variables for storing the plot data
meanVel = nan(3, 3); % Mean velocity
clVel = nan(3, 3); % Confidence interval lower bound
cuVel = nan(3, 3); % Confidence interval upper bound
n = nan(3,1); % Number of examples

% % NOTE: These are defined above and can be commented out/ removed.
% % Define the post-stimulus window
% postStimWindow = t_ms>0 & t_ms<postStimTime;
% % Define the index for selecting trials above minimum forward speed
% vfBefore = squeeze(nanmean(V(t<0,2,:),1));
% vfIdx = vfBefore > vfMin;

% Loop through each limb hit type (fore, mid, hind)
for ind = 1:3
    
    % Select the trials associated with the given limb type AND above threshold walking speed
    idx = hitSym(:,ind) & vfIdx;
    n(ind) = nnz(idx); % Count the number of examples
    
    % Select the post stimulus velocities
    Vsel = squeeze(mean(Vsym(postStimWindow,:,idx),1));
    
    % Convert forward velocities to fold-change
    Vsel(2,:) = Vsel(2,:) ./ vfBefore(idx)';
    
    % Compute the mean post-stimulus velocities
    meanVel(ind,:) = nanmean(Vsel,2);
    
    % Bootstrap to generate confidence intervals
    ci = bootci(nboot, {@(x) nanmean(x,1), permute(Vsel, [2,1])}, 'Alpha', bootAlpha);
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
    % Plot errorbars based on the confidene interval
    errorbar(1:3, meanVel(:,ind), meanVel(:,ind)-clVel(:,ind), cuVel(:,ind)-meanVel(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0, 'HandleVisibility','off');
    
    % Plot a line at y=1 for fold change in forward walking
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
    ConfAxis('fontSize', 16);
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
%% Timeseries of mean relative phase

% Check if this plot should be suppressed
if ~supressPlots
    
    % Compute relative phases
    phiRel = mod(Phisym(:,[1,2,3],:) - Phisym(:,[4,5,6],:),2*pi);

    % Initialize variables for storing computed mean relative phases
    meanRelPhase = nan(length(t), 3, 3);
    clRelPhase = nan(length(t), 3, 3);
    cuRelPhase = nan(length(t), 3, 3);
    n = nan(3,1);

    % % NOTE: These are defined above and can be commented out/ removed.
    % % Define the indexing to remove trials where the forward speed is below a threshold
    % vfBefore = squeeze(nanmean(V(t<0,2,:),1));
    % vfIdx = vfBefore > vfMin;

    tic;
    % Loop through each of the limb hit types (fore, mid, hind)
    for ind = 1:3
        % Select the trials associated with the given limb hit type
        idx = hitSym(:,ind) & vfIdx & ~squeeze(any(any(isnan(phiRel),2),1));
        % Get the number of trials
        n(ind) = nnz(idx);
        % Select the limb relative phases associated with the current condition
        phiSel = phiRel(:,:,idx);

        % Compute the mean prestimulus relative phase
        phiBef = mod(circmean(phiSel(t<0,:,:),1),2*pi);

        % Compute the difference between the current timepoint relative phase 
        % and mean prestimulus relative phase
        % NOTE: Adding pi moves this away from the branch point. It is
        % subtracted later.
        phiSel = mod(phiSel - phiBef + pi,2*pi);

        % Get the average difference between prestimulus relative phase and
        % current phase at all timepoints
        meanRelPhase(:,:,ind) = mod(circmean(phiSel,3), 2*pi);

        % Compute confidence interval error bars
        ci = bootci(nboot, {@(x) mod(circmean(x,1),2*pi), permute(phiSel, [3,1,2])});
        clRelPhase(:,:,ind) = ci(1,:,:,:);
        cuRelPhase(:,:,ind) = ci(2,:,:,:);
    end

    % Convert from radians to cycle fractions
    meanRelPhase = (meanRelPhase - pi) / (2*pi);
    clRelPhase = (clRelPhase - pi) / (2*pi);
    cuRelPhase = (cuRelPhase - pi) / (2*pi);
    toc;

    % Define titles and legend labels
    seriesLabels = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};
    legendStr = {'C1-I1','C2-I2','C3-I3'};

    % Loop through each hit type and create a separate plot for each
    for ind = 1:3
        MakeFigure;
        % Plot the mean relative phase
        PlotAsymmetricErrorPatch(t_ms, meanRelPhase(:,:,ind), clRelPhase(:,:,ind), cuRelPhase(:,:,ind), limbCorder);

        % Format the plot
        ylim([-0.15, 0.15]);
        yticks(-0.15:0.05:0.15);
        hold on;
        plot([-250,250], [0,0], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        plot([0,0],[-0.15, 0.15], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        xlabel('time (ms)');
        ylabel('relative phase advance (cycles)');
        legend(legendStr);
        title(seriesLabels{ind});
        ConfAxis('fontSize', 16);
        axis('square');
    end
    
end
%% Mean relative phase advance

% Check if this plot should be suppressed
if ~supressPlots

    % Compute relative phases
    phiRel = mod(Phisym(:,[1,2,3],:) - Phisym(:,[4,5,6],:),2*pi);

    % Compute mean relative phases
    meanRelPhase = nan(3, 3);
    clRelPhase = nan(3, 3);
    cuRelPhase = nan(3, 3);
    n = nan(3,1);

    % % NOTE: This is identically calculated above
    % % Compute an index for restricting based on prestimulus speed
    % vfBefore = squeeze(nanmean(V(t<0,2,:),1));
    % vfIdx = vfBefore > vfMin;

    tic;
    % Loop through each of the limb hit types (fore, mid, hind)
    for ind = 1:3
        % Define an index for selecting hits for this category
        idx = hitSym(:,ind) & vfIdx & ~squeeze(any(any(isnan(phiRel),2),1));
        % Get the number of hits
        n(ind) = nnz(idx);
        % Select the relevant phase data
        phiSel = phiRel(:,:,idx);

        % Compute the average mean phase
        phiBef = mod(circmean(phiSel(t<0,:,:),1),2*pi);
        % Compute the difference between the average prestimulus phase and the
        % poststimulus mean
        phiSel = mod(circmean(mod(phiSel(postStimWindow,:,:) - phiBef + pi,2*pi), 1),2*pi);
        % Average over examples to get the mean phase differnce
        meanRelPhase(:,ind) = mod(circmean(phiSel,3), 2*pi);

        % Compute some error bars
        ci = bootci(nboot, {@(x) mod(circmean(x,1),2*pi), permute(phiSel, [3,1,2])}, 'Alpha', bootAlpha);
        clRelPhase(:,ind) = ci(1,:,:,:);
        cuRelPhase(:,ind) = ci(2,:,:,:);
    end

    % Convert from radians to cycle fractions
    meanRelPhase = (meanRelPhase - pi) / (2*pi);
    clRelPhase = (clRelPhase - pi) / (2*pi);
    cuRelPhase = (cuRelPhase - pi) / (2*pi);
    toc;

    % Define labels for the plot below
    seriesLabels = {sprintf('forelimb hit (N=%d)', n(1)),sprintf('midlimb hit (N=%d)', n(2)), sprintf('hindlimb hit (N=%d)', n(3))};
    legendStr = {'C1-I1','C2-I2','C3-I3'};

    % Make a single plot for each of the mean phase difference barplots
    MakeFigure;

    % Loop through each hit type as a subplot
    for ind = 1:3
        subplot(1,3,ind);
        % Plot the data
        bar(1:3, meanRelPhase(:,ind), 'EdgeColor','none', 'FaceColor', 'flat', 'CData', limbCorder);
        hold on;
        % Add confidence interval error bars
        errorbar(1:3, meanRelPhase(:,ind), meanRelPhase(:,ind)-clRelPhase(:,ind), cuRelPhase(:,ind)-meanRelPhase(:,ind), ' k', 'LineStyle','none','LineWidth', dataWidth, 'CapSize', 0, 'HandleVisibility','off');
        % Format the plot
        ylim([-0.15, 0.15]);
        yticks(-0.15:0.05:0.15);
        ylabel('Relative Phase Advance (Cycles)');
        xticks(1:3);
        xticklabels(legendStr);
        title(seriesLabels{ind});
        ConfAxis('fontSize', 16);
        axis('square');
    end
    
end
%% Average timeseries of limb separation

% Check if this plot should be suppressed
if ~supressPlots
    
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

        % Compute the fractional change in contralateral limb distance
        Dsel = (Dsel - DselBefore) ./ DselBefore;

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
        plot([0 0],[-0.15,0.15], '--k', 'linewidth', dashWidth, 'HandleVisibility','off');
        hold off;

        % Format the plot
        xlabel('time (ms)');
        ylabel('limb separation (fractional change)');
        legend(seriesLabels);
        title(legendStr{ind});
        ylim([-0.15,0.15]);
        yticks(-0.15:0.05:0.15);
        ConfAxis('fontSize', 16);
        axis('square');
    end
end
%% Average limb separation

% Check if this plot should be suppressed
if ~supressPlots

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
        % Compute the fractional change in limb distances
        Dsel = nanmean((Dsel(postStimWindow,:,:) - DselBefore) ./ DselBefore,1);

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
        ylabel('Limb Separation (Fractional Change)');
        title(legendStr{ind});
        ylim([-0.15,0.15]);
        yticks(-0.15:0.05:0.15);
        ConfAxis('fontSize', 16);
        axis('square');
    end
end
