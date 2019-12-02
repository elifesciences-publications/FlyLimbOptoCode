% This script tests the temporal alignment of flashes in the projector and
% camera

%% Check that the flashes are happening when the camera says that they are happening

% Load the camera data
load('cameraData.mat');

% Plot the mean pixel intensity vs. the gpioPixel values
flashOn = zeros(size(gpioPixels,1),1);
flashOn(gpioPixels(:,1)==240) = 1;
makeFigure;
hold on;
yyaxis left
plot(timeStamps, meanPixels);
yyaxis right
plot(timeStamps, flashOn);
legend({'Pixel Intensity', 'Camera Detection'});


%% Load the camera clock and the projector clock data

% Load the projector data
load('projectorData.mat','mapIdx', 'timestamp');
projTime = timestamp;
clear timestamp

% Load the camera data
load('cameraData.mat','gpioPixels', 'timeStamps', 'meanPixels');
camTime = timeStamps;
clear timeStamps

%% Check that the camera clock and the projector clock are advancing at the same rate
% Unwind the camera timestamps
% Define what the max time value is for the cycling camera timestamps 
max_time_value = 128; % Note: This is in seconds

% Calculate the differences in time 
diffTime = diff(camTime);

% Find all the locations in diffTime where the camTime decrement
% NOTE: Could use 0 as the bound but this is more robust in case there is some weirdness with the values
decrease_locs = find(diffTime < -1*(max_time_value-20)) + 1;

% Loop through the decreased locations and add the max values to all subsequent camTime
for n = 1:length(decrease_locs) 
    loc = decrease_locs(n);
    camTime(loc:end) = camTime(loc:end) + max_time_value;
end

% Find the alignment pulses in the projector time
projAlignTimes = projTime(mapIdx ==-1);

% Find the alignment pulses in camera time
temp = camTime(gpioPixels==240);
camAlignTimes = [temp(1); temp(end)];

projTimeDur = projAlignTimes(2) - projAlignTimes(1);
camTimeDur = camAlignTimes(2) - camAlignTimes(1);

% % % Rescaling to remove the effect of one clock running faster than the other % % %
projTime = projTime * (camTimeDur/projTimeDur);
projAlignTimes = projTime(mapIdx ==-1); % Recalculate the new alignment times
% % % Rescaling to remove the effect of one clock running faster than the other % % %

% Convert the projector time to camera time by subtracting off the offset at the first alignment pulse
offset = projAlignTimes(1) - camAlignTimes(1);
projTime_aligned = projTime - offset;

% Define a logical for the flash to be on
mapLogical = (mapIdx ~= 0);

% Convert the gpioPixels in to a logical
gpioLogical = (gpioPixels(:,1) == 240);

% Convert meanPixels into a logical
meanLogical = (meanPixels >= 140);

%% Plot the gpioLogical and mapLogical as a function of cameraTime
makeFigure;
hold on;
hold on;
yyaxis left
plot(projTime_aligned, mapLogical);
yyaxis right
% plot(camTime, gpioLogical);
plot(camTime, meanLogical);
legend({'Projector', 'Camera'});

%% Find all the paired occurrences of flashes in the camera and projector

% Calculate the absolute differences in time between the projector and the camera timestamps 
active_proj_time = projTime_aligned(mapLogical);
active_cam_time = camTime(meanLogical);
proj_array = repmat(active_proj_time',length(active_cam_time),1);
cam_array = repmat(active_cam_time,1,length(active_proj_time) );
diff_array = abs(proj_array - cam_array);

% Find the closest camera timestamp for every projector timestamp
min_diff_vec = min(diff_array,[],1);

% Remove the examples where the timestamp differs by more the 800 ms as this means that the camera missed the detection of this pulse
min_diff_vec(min_diff_vec >.8) = [];

% Plot the result
% makeFigure;
hold on;
plot(min_diff_vec);
ylabel('\Delta Time (Secs)');
xlabel('Projector Pulse Index');
title('Differences in Projector and Camera Timing'); 
ConfAxis;