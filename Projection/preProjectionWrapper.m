function [ outMap, mapIdx, framesOn, frameOff, projBitMap, cameraBitMap, alignPulseOutMap, framesPerUp ] = preProjectionWrapper( varargin )
% This function is a wrapper for the code used to generate an outMap as
% input into the projectionWrapper. Prior to running this function need to 
% calculate tform using calculateCoordinateTransformationMatrix.m

%% 
%CHECK THAT YOU'VE REMOVED THE CONSTANT STIMULUS THING BELOW

%% 
% NOTE: A second function will take the output of this function and run the projector

% Define some default values

% % Define the experiment duration
% % NOTE: The camera will capture 4000s of data. We want the projector stimulation to run to completion before the camera stops so that we get both alignment pulses in the video.
% NOTE: With the desired buffer of 5 seconds on each side on an alignement pulse we have slightly less than 60 seconds less projection than camera time. 
% [60 sec minus the 2*16.6 msec (the two alignment pulses) and 2*.05 sec (Pauses at each alignment pulse)] than the total amount of video recording using 3920.
expDuration = 3920; % seconds 

% % % % Parameters for Positive Control % % %
% expDuration = 246; % seconds
% expDuration = 113; 

% Define an alignment pulse flag
alignmentOn = true;

% Define a buffer added to the front and end of the experiment startEndAlignmentBuffer
startEndAlignmentBuffer = 5; % seconds

% Evaluate any name-value argument pairs
for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

% Select the directory of the desired alignment file
folder = uigetdir('F:\Data\_FeatureExtractionFiles\ProjectorAlignmentFiles\Alignment_Files\','Select the folder for the desired alignment file.');
% folder = uigetdir('C:\Users\labuser\Desktop\Brian_SL_Tests\Alignment_Files','Select the folder for the desired alignment file.');
tformPath = fullfile(folder,'tform.mat');

% Generate the desired camera bitMap
% [ cameraBitMap, framesOn, frameOff, timesOn, timeOff, framesPerUp, shiftsPerDim ] = generateSquareLatticeBitmaps_SingleFrames(varargin{:});
[ cameraBitMap, framesOn, frameOff, timesOn, timeOff, framesPerUp, shiftsPerDim ] = generatePoissonBitmap_SingleFrames(varargin{:});
% [ cameraBitMap, framesOn, frameOff, timesOn, timeOff ] = generateSquareLatticeBitmaps(varargin{:});

% % Optional: Add a full-field activation to the cameraBitMap
% cameraBitMap = cat(3, cameraBitMap, ones(size(cameraBitMap,1),size(cameraBitMap,2)));
% shiftsPerDim = size(cameraBitMap,3)^(1/2);

% Load the transformation matrix that converts the desired camera bitMap to the projector bitMap
% NOTE: Need to generate tform and save in a location in the repo
load(tformPath);

% Transform camera's bitMap into the necessary projector bitMap
[ projBitMap, cameraBitMap ] = generateProjectorBitmap( cameraBitMap, tform );


% % TEMPORARY METHOD (skips the conversion from projector to camera coordinates)
% cameraBitMap = [];
% [ projBitMap, framesOn, frameOff, timesOn, timeOff, framesPerUp, shiftsPerDim ] = generateSquareLatticeBitmaps_SingleFrames(varargin{:});

% Take the projector bitMap and generate the necessary outMap(s) for projection
[ outMap ] = generateOutmapsFromSingleFrameBitmaps( projBitMap, framesPerUp, framesOn, frameOff );
% % [ outMap ] = generateOutmaps( projBitMap, framesPerUp );

% Calculate the number of projector flips in the experiment
% NOTE: rounded down in case the expDuration given is not an integer
numFlips = floor(expDuration * 60); % Projector runs at 60Hz

% Generate the list of outMap indices that will be projected during the experiment
mapIdx = zeros(numFlips,1);
[ mapIdx ] = generateOutMapList( numFlips, mapIdx, framesOn, frameOff, framesPerUp, shiftsPerDim);
mapIdx = int32(mapIdx);

% %% TESTING WITH CONSTANT STIMULUS
% 
% % This is just a temporary hack. Don't use outside of this single test
% mapIdx(mapIdx == 1) = 2; % kills the partial stimulus after the flip
% 
% %% 

if alignmentOn
    
    % Add alignment pulses and the corresponding buffers to the mapIdx list
    buffer = zeros(startEndAlignmentBuffer*60,1); % Projector runs at 60Hz
%     buffer2 = zeros(ceil(frameOff/framesPerUp),1); % Old: Buffer between alignment stimulus and first activation. Duration of a normal off cycle
    mapIdx = [buffer;-1;buffer;mapIdx;buffer;-1;buffer];

    % Generate the alignPulseOutMap
    [ alignPulseOutMap ] = getAlignPulseOutMap( varargin{:} );
    
end

end

