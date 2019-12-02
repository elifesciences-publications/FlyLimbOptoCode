function [ bitMap, framesOn, frameOff, timesOn, timeOff, framesPerUp, shiftsPerDim ] = generatePoissonBitmap_SingleFrames( varargin )
% This function generates Poisson bitmpas which contain single
% frames for each of the desired conditions. It is an alternative to 
% generateSquareLatticeBitmaps_SingleFrames.m and will be used with
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
cRad = 5; % px, radius of the circle % (5 Originally)
cDist = 64; % (64 Originally)
lambda = 75; % Lambda of Poisson Distribution (75 Originally)
timesOn = 30; % ms (s * (10^(-3))) % (10 Originally)
timeOff = 500; % ms % (500 Originally)
shiftsPerDim = 3;
shiftSize = floor(h/(shiftsPerDim*2)); % px

% % Alignment Test Parameters
% h = 512; % px
% w = 640; % px
% cRad = 5; % px, radius of the circle
% cDist = 64;
% lambda = 75; % Lambda of Poisson Distribution
% timesOn = 100000; % ms (s * (10^(-3)))
% timeOff = 0; % ms
% shiftsPerDim = 1;
% shiftSize = floor(h/(shiftsPerDim*2)); % px

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

% Create a matrix for holding all the activations
circMap = zeros(mapH, mapW);

% Create a circular activation with radius cRad in a box with length cDist
cy = floor(cDist/2);
cx = floor(cDist/2);
[yy, xx] = meshgrid(1 : cDist);
unitmap = (sqrt((yy-cy).^2+(xx-cx).^2)<=cRad);

% Get the number of activations in this region
numActive = poissrnd(lambda,1,1);

% Generate a vector of coordinates for each of the activations
x_locs = randi([1 mapH], 1, numActive);
y_locs = randi([1 mapW], 1, numActive);

% Get the top left corner of each activation patch
x_start = x_locs - floor(cDist/2);
y_start = y_locs - floor(cDist/2);

x_end = x_start + (cDist-1);
y_end = y_start + (cDist-1);

for n = 1: numActive 
    
    % Check the boundaries and adjust the size of the activation patch
    cur_unitmap = unitmap;
    if x_end(n) > mapH
        cur_unitmap = cur_unitmap(1:end-(x_end(n)-mapH),:);
        x_end(n) = mapH; 
    elseif x_start(n) < 1
        cur_unitmap = cur_unitmap(1+(1-x_start(n)):end,:);
        x_start(n) = 1;
    end
        
    if y_end(n) > mapW
        cur_unitmap = cur_unitmap(:,1:end-(y_end(n)-mapW));
        y_end(n) = mapW;
    elseif y_start(n) < 1
        cur_unitmap = cur_unitmap(:,1+(1-y_start(n)):end);
        y_start(n) = 1;
    end
    
    % Place the activation in the map
    circMap(x_start(n):x_end(n),y_start(n):y_end(n)) = (circMap(x_start(n):x_end(n),y_start(n):y_end(n))|cur_unitmap);
    
end

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


