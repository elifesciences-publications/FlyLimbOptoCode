function [t, t_ms,  hit, hitLocX, hitLocY, V, hitSym, Vsym, Lsym, Phisym, D, downSym] = ExtractLimbHits(newData, winlen, useRandomTrigger, limbHitsOnly, smoothCentroid, smoothPhase, removeMultipleHits, mmPerPix, semiMajorAxis, semiMinorAxis, xShift, yShift)
% Utility function to extract limb hit data

%% Set variable lists

% Define the list of limb prefixes
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the list of hit variables
hitVarList = strcat(limbList, '_hit');

% Define lists for the limb position and limb phase related variables
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');
phaseVarList = strcat('InstantaneousPhase_', limbList, 'y');
downVarList = strcat(limbList, '_down_cam');

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

% Update the limb phase variables, if smoothing
if smoothPhase
    phaseVarList = strcat('smooth_',phaseVarList);
end

% Append all the variable lists together
varList = [velVarList, hitVarList, limbVarListX, limbVarListY,...
    phaseVarList, comVarList, hitLocVarListX, hitLocVarListY, downVarList];


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
down = S(:,43:48,:); % Swing/stance variables

L = -L; % Reverse both limb position axes

V(:,1,:) = rad2deg(V(:,1,:)); % Convert yaw rate to degrees

% Compute the distance between contralateral limbs
D = hypot(L(:, 1:3, :) - L(:, 4:6,:), L(:, 7:9,:) - L(:, 10:12, :));

if useRandomTrigger
    hit = rand(size(H,3),size(H,2)) < 0.2;
else
    hit = squeeze(H(t==0,:,:))' > 0;
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

% Remove hits that do not satisfy this criterion
includeIdx = any(outsideBodyEllipse,2);
if ~all(any(outsideBodyEllipse,2),1)
    hit = hit(includeIdx,:);
    hitLocX = hitLocX(includeIdx,:);
    hitLocY = hitLocY(includeIdx,:);
    V = V(:,:,includeIdx);
    L = L(:,:,includeIdx);
    Phi = Phi(:,:,includeIdx);
    D = D(:,:,includeIdx);
    down = down(:,:,includeIdx);
end

%% Symmetrize left and right hits

% Define indicators
leftHit = squeeze(any(hit(:,1:3),2));
rightHit = squeeze(any(hit(:,4:6),2));

% Define symmetric version of each of the groups of variables
hitSym = squeeze((hit(:,4:6) & rightHit) | (hit(:,1:3) & leftHit));
Vsym = V; % Body velocity
Vsym(:, [1,3], leftHit) = -Vsym(:, [1,3], leftHit);
Lsym = L; % Limb positions
Lsym(:, :, leftHit) = L(:, [4,5,6,1,2,3,10,11,12,7,8,9], leftHit);
Lsym(:,1:6,leftHit) = -Lsym(:,1:6, leftHit);
Phisym = Phi; % Limb phase
Phisym(:,:,leftHit) = Phisym(:,[4,5,6,1,2,3],leftHit);
downSym = down; % Swing/stance
downSym(:,:,leftHit) = down(:,[4,5,6,1,2,3],leftHit);

end

