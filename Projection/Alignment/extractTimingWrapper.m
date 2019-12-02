function [ timeStamps, gpioPixels, meanPixels ] = extractTimingWrapper()
% This function extracts the timing information from a set of videos and
% aligns it to the corresponding camera stimulus

% Get the list of videos corresponding to this particular set of flies
startDir = 'C:\Users\ClarkLab\Desktop\SingleLegActivation';
videoDir = uigetdir(startDir, 'Select a folder of videos.');

% Extract the list of avi files from the video directory
fileList = dir(strcat(videoDir, '\*.avi'));

% % Remove the non-file elements from the fileList structure
% fileList = fileList(~[fileList.isdir]);

timeStamps = [];
gpioPixels = [];
meanPixels = [];

% Loop through all the videos in the directory and run the pixel extraction
for n = 1:length(fileList)
    
    % Select the current file in the directory
    curFile = fileList(n).name;
    
    % Extract the current set of time and gpio states
    [ curTime, curPixels, curMeanPix ] = extractCameraTiming( fullfile(videoDir,curFile));
    
    % Append the current set of values to the full list
    timeStamps = [timeStamps;curTime];
    gpioPixels = [gpioPixels;curPixels];
    meanPixels = [meanPixels;curMeanPix]; % For testing
    
    % Print a status update
    fprintf('Extracted timing from video %i of %i. \n', n, length(fileList));
    
end

% % Unwind the camera timestamps so that the timestamps only increment

% Define what the max time value is for the cycling camera timestamps 
max_time_value = 128; % Note: This is in seconds

% Calculate the differences in time 
diffTime = diff(timeStamps);

% Find all the locations in diffTime where the timestamps decrement
% NOTE: Could use 0 as the bound but this is more robust in case there is some weirdness with the values
decrease_locs = find(diffTime < -1*(max_time_value-20)) + 1;

% Loop through the decreased locations and add the max values to all subsequent timestamps
for n = 1:length(decrease_locs)
    val = decrease_locs(n);
    timeStamps(val:end) = timeStamps(val:end) + max_time_value;
end

end