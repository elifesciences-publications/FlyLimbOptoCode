function [meanVel, clVel, cuVel, n] = computeOptoTimeSeries(newData)
% This function takes as input the raw datasets and returns time series
% variables for plotting

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

%% Set variable lists

% Define the list of limb prefixes
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the list of hit variables
hitVarList = strcat(limbList, '_hit');

% Define lists for the limb position and limb phase related variables
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');

% Define the the list of body velocity related variables
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
if smoothCentroid
    velVarList = strcat('smooth_', velVarList);
end

% Append all the variable lists together
varList = [velVarList, hitVarList, limbVarListX, limbVarListY];

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

%% Extract needed data

% Convert the trigger variable into an array of time series centered on the trigger
S = getSlicesSimple(newData, varList, winlen, trigger);
S = cat(3, S{:});

% Extract out each of the groups of variables from the slice data structure
V = S(:, 1:3, :); % Body velocity
H = S(:, 4:9, :); % Hit logicals
L = S(:, 10:21, :); % Limb positions in mm

L = -L; % Reverse both limb position axes

V(:,1,:) = rad2deg(V(:,1,:)); % Convert yaw rate to degrees

% Compute the distance between contralateral limbs
D = hypot(L(:, 1:3, :) - L(:, 4:6,:), L(:, 7:9,:) - L(:, 10:12, :));

% Define the time variable that will be used for plotting
t = (-winlen:winlen)';
t_ms = t./(0.15);

if useRandomTrigger
    hit = rand(size(H,3),size(H,2)) < 0.2;
else
    hit = squeeze(H(t==0,:,:))';
end

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

%% Timeseries of average centroid kinematics

% Initialize some empty variables
meanVel = nan(length(t), 3, 3); % Mean velocity
clVel = nan(length(t), 3, 3); % Lower confidence interval
cuVel = nan(length(t), 3, 3); % Upper confidence interval
n = nan(3,1); % Number of examples

% Compute an indicator for restricting based on initial forward speed
vfBefore = squeeze(nanmean(V(t<0,2,:),1));
vfIdx = vfBefore > vfMin;

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

end

