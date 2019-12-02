function [ stepDensity, binCenters ] = calcNextStepDensity_BodyHits( newData, nBins )
% This function takes as input raw datasets from optogenetic experiments
% and returns a spatial distribution of next steps (for body hits)

%% Set some parameters for the analysis

% Window length
winlen = 36;

% Various analysis parameters
vfMin = 3;

%% Set variable lists

% Define the list of variables
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the associated hit variables
hitVarList = {'bodyHitOnset', 'bodyHit','bodyHitX_mm','bodyHitY_mm','bodyEllipseDist_mm', 'videoID'};

% Define all limb related variables for analysis
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');

% Define the body related variables for the analysis
comVarList = {'Orient_Rad_FullRotation', 'xCOM','yCOM'};

% Define the lists for the limb position up/down variables
limbUpDownList = strcat(limbList, '_down_cam');

% Define the body velocity variables
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec','translationalSpeed_mmPerSec'};

% Get the full variable list for inclusion in the analysis
varList = [velVarList, hitVarList, limbVarListX, limbVarListY,...
           comVarList, limbUpDownList];

%% Make the trigger

% Extract all the body hits from the dataset
hitData = newData.bodyHitOnset;
trigger = any(hitData, 2);

% Define all the time variables
t = (-winlen:winlen)';
t_ms = t./(0.15);

% Define the post stimulus index range
postStimWindow = t_ms>0;

%% Extract needed data

% Get the variables associated with individual time series trials
S = getSlicesSimple(newData, varList, winlen, trigger);
S = cat(3, S{:});

% Partition the data into separate variable types
V = S(:, 1:3, :); % Body velocity data
H = S(:, 4:9, :); % Body hit related variables
L = S(:, 10:21, :); % Limb related variables
CM = S(:,22:24,:); % Body related position data
U = S(:, 25:30, :);% Limb up/down logicals 

% Invert the limb position variables (Changes the coordinate system)
L = -L;

% Convert the rotational velocities from radians to degrees
V(:,1,:) = rad2deg(V(:,1,:)); 

%% Group hits by location

% Extract hit locations
hitLocX = squeeze(H(t==0, 3, :));
hitLocY = squeeze(H(t==0, 4, :));

% Group hits by right and left head, thorax, abdomen...
segmentIdx = [
    (hitLocX > 0) & (hitLocY > 1),...
    (hitLocX > 0) & (hitLocY > 0) & (hitLocY < 1),...
    (hitLocX > 0) & (hitLocY < 0),...
    (hitLocX < 0) & (hitLocY > 1),...
    (hitLocX < 0) & (hitLocY > 0) & (hitLocY < 1),...
    (hitLocX < 0) & (hitLocY < 0)];

%% Symmetrize left and right hits

% Create an indicator for left side hits
leftHit = any(segmentIdx(:,4:6),2);

% Symmetrize the segment hit data
hitSym = segmentIdx(:,1:3);
hitSym(leftHit,1:3) = segmentIdx(leftHit,4:6);

% Symmetrize the corresponding body velocity data by reversing v_lateral
% and v_rotational for cases of left hits
Vsym = V;
Vsym(:, [1,3], leftHit) = -Vsym(:, [1,3], leftHit);

% Swap the limb position data for left hits to symmetrize
Lsym = L;
Lsym(:, :, leftHit) = L(:, [4,5,6,1,2,3,10,11,12,7,8,9], leftHit);
Lsym(:,1:6,leftHit) = -Lsym(:,1:6, leftHit);

% Symmetric version of each of the UP/DOWN time-series 
% NOTE: Perform same operation as in limb positions (excl. sign inversion)
% NOTE: This converts everything to look like right side hits
Usym = U;
Usym(:, :, leftHit) = U(:, [4,5,6,1,2,3], leftHit);

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