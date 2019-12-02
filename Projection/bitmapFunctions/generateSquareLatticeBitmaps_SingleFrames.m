function [ bitMap, framesOn, frameOff, timesOn, timeOff, framesPerUp, shiftsPerDim ] = generateSquareLatticeBitmaps_SingleFrames( varargin )
% This function generates square lattice bitmaps which contain single
% frames for each of the desired conditions. It replaces
% generateSquareLatticeBitmaps.m and will be used with
% generateOutmapsFromSingleFrameBitmaps.m for projection.

% 604 --> Lightcrafter Projector Height in Pixels
% 684 --> Lightcrafter Projector Width in Pixels
% 912  --> LCR 4500 Projector Height in Pixels
% 1140 --> LCR 4500 Projector Width in Pixels
% 512 --> Camera Height in Pixels
% 640 --> Camera Width in Pixels

tic;
framesPerCycle = 60;
framesPerUp = 24;
frameRate = framesPerCycle * framesPerUp; % We project at 1440 Hz

% Set defaults
h = 512; % px 
w = 640; % px
cDist = 2; % px, distance between circle centers %64
cRad = 5; % px, radius of the circle
timesOn = 30; % ms (s * (10^(-3))) [3,6,12]
timeOff = 4000; % ms
% timesOn = 10000; % ms (s * (10^(-3)))
% timeOff = 0; % ms
shiftsPerDim = 1;
shiftSize = floor(cDist/shiftsPerDim);% px

% % % Parameters for positive control % % %
% h = 512; % px 
% w = 640; % px
% cDist = 2; % px, distance between circle centers %64
% cRad = 5; % px, radius of the circle
% timesOn = [3,6,12]; % ms (s * (10^(-3))) [3,6,12]
% timeOff = 1500; % ms
% % timesOn = 10000; % ms (s * (10^(-3)))
% % timeOff = 0; % ms
% shiftsPerDim = 1;
% shiftSize = floor(cDist/shiftsPerDim);% px
% % % End parameters for positive control % % %

% % % Constant Stimulus Test
% h = 512; % px 
% w = 640; % px
% cDist = 100; % px, distance between circle centers %64
% cRad = 10; % px, radius of the circle
% timesOn = 1000; % ms (s * (10^(-3))) *makes for no remainder
% timeOff = 0; % ms
% shiftsPerDim = 1;
% shiftSize = floor(cDist/shiftsPerDim);% px

% % % % Testing 
% cDist = 64; % px, distance between circle centers %64
% cRad = 4; % px, radius of the circle
% timesOn = 150000; % ms (s * (10^(-3)))
% timeOff = 0; % ms
% shiftsPerDim = 1;
% shiftSize = floor(cDist/shiftsPerDim);% px



% Evaluate any name-value argument pairs
for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

% Convert seconds to frames
framesOn = floor(timesOn * (10^(-3)) * frameRate);
frameOff = floor(timeOff * (10^(-3)) * frameRate);

% Update to get the true timesOn and timeOff
timesOn = (framesOn/frameRate)*(10^3);
timeOff = (frameOff/frameRate)*(10^3);

% Calculate the remainder frames needed for each duration
remainFramesOn = mod(framesOn,framesPerUp);

% Add an full frame entry for each of the 
remainFramesOn = [remainFramesOn,framesPerUp];

% Set the height and width of the map in pixels
% NOTE: Along each dimension the you must have a minimum of h and w pixels 
% respectively. The dimensions will be integer numbers of unitmap tiles so
% that the map can be easily circularly permuted.
if (mod(h,cDist) == 0) && (mod(w,cDist) == 0)
    mapH = h;
    mapW = w;
elseif (mod(h,cDist) == 0) && (mod(w,cDist) ~= 0)
    mapH = h;
    mapW = (floor(w/cDist) + 1)*cDist;
elseif (mod(w,cDist) == 0)
    mapW = w;
    mapH = (floor(h/cDist + 1))*cDist;
else
    mapH = (floor(h/cDist) + 1)*cDist;
    mapW = (floor(w/cDist) + 1)*cDist;
end

% Create a circular activation with radius cRad in a box with length cDist
cy = floor(cDist/2);
cx = floor(cDist/2);
[yy, xx] = meshgrid(1 : cDist);
unitmap = (sqrt((yy-cy).^2+(xx-cx).^2)<=cRad);

% Tile this activation to create a full resolution single frame
% circMap = uint8(repmat(unitmap, (mapH/cDist), (mapW/cDist)));
circMap = repmat(unitmap, (mapH/cDist), (mapW/cDist));

% Shuffle the map to produce the desired patterns 
% Default is to tile the map such that all activations would touch each other along the outer radius
counter = 1;
for idxR = 1:shiftsPerDim
    for idxC = 1:shiftsPerDim
        % Generate the desired shift
        circMap(:,:,counter) = circshift(circMap(:,:,1),shiftSize*(idxR-1),1);
        circMap(:,:,counter) = circshift(circMap(:,:,counter),shiftSize*(idxC-1),2);
        % Increment the counter
        counter = counter + 1;
    end    
end

% Crop the map to the correct size 
bitMap = circMap(1:h,1:w,:);

fprintf('%d bitmap frames generated in %f seconds.\n',size(bitMap,3),toc);


end

