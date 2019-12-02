function [ stepDensity, binCenters ] = calcNextStepDensity_LimbHits( newData, nBins )
% This function takes as input raw datasets from optogenetic experiments
% and returns a spatial distribution of next steps

%% Set some parameters for the analysis

% Window length
winlen = 36;

% Whether to use smoothed or raw variables
smoothCentroid = false;

% Whether to remove instances where multiple limbs are hit
removeMultipleHits = true;

% Various analysis parameters
useRandomTrigger = false;
vfMin = 3;

%% Set variable lists

% Define the list of limb prefixes
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the list of hit variables
hitVarList = strcat(limbList, '_hit');

% Define lists for the limb position variables
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');

% Define the lists for the limb position up/down variables
limbUpDownList = strcat(limbList, '_down_cam');

% Define the the list of body velocity related variables
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};
if smoothCentroid
    velVarList = strcat('smooth_', velVarList);
end

% Append all the variable lists together
varList = [velVarList, hitVarList, limbVarListX, limbVarListY, limbUpDownList];

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
U = S(:, 22:27, :);% Limb up/down logicals 

L = -L; % Reverse both limb position axes

V(:,1,:) = rad2deg(V(:,1,:)); % Convert yaw rate to degrees

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

% Symmetric version of each of the UP/DOWN time-series 
% NOTE: Perform same operation as in limb positions (excl. sign inversion)
% NOTE: This converts everything to look like right side hits
Usym = U;
Usym(:, :, leftHit) = U(:, [4,5,6,1,2,3], leftHit);

% Define the post stimulus index range
postStimWindow = t_ms>0;

%% Compute distributions of next step locations

% Define the binning of the joint distribution
binEdges = -3:(6/nBins):3;
% Get the bin centers
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

% Initialize some empty variables
stepDensity = nan(nBins,nBins,3); % image height x image width x # hit types
xPos = cell(3,1);
yPos = cell(3,1);
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
    
    % Select the limb position and limb up/down data associated with the given hit type
    % NOTE: dimensions --> timepoints, variables, trials
    Lsel = Lsym(:,:,idx);
    Usel = Usym(:,:,idx);
    
    % Split the limb position data into X and Y variables 
    Xsel = Lsel(:,1:6,:);
    Ysel = Lsel(:,7:12,:);
   
    % Get the limb positions data at first touchdown post stimulation for each limb
    [ xPos{ind}, yPos{ind} ] = findLocsTD(Usel, Xsel, Ysel, postStimWindow);

    % Compute the joint distribution for the limb
    stepDensity(:,:,ind) = histcounts2(xPos{ind}, yPos{ind}, binEdges, binEdges, 'normalization','pdf');

end

end

