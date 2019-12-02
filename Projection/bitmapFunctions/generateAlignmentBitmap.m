function [ alignBitMap, alignOutMap, xList, yList, framesPerUp ] = generateAlignmentBitmap( varargin )
% This function generates the bitmap frames necessary for aligning the
% projector and the camera.

% 604 --> Lightcrafter Projector Height in Pixels
% 684 --> Lightcrafter Projector Width in Pixels
% 912  --> LCR 4500 Projector Height in Pixels
% 1140 --> LCR 4500 Projector Width in Pixels
% 512 --> Camera Height in Pixels
% 640 --> Camera Width in Pixels

% Set defaults
framesPerUp = 24;
h = 912; % px 
w = 1140; % px

% % R with +i spacing
% % (NOTE: In projector bitmap, y corresponds to rows, x corresponds to columns)
% y0 = 380; %380 %270
% x0 = 220; %320 (for i=8)
% i = 40;
% yList = [ y0, y0,   y0,       y0+i, y0+i,     y0+(2*i), y0+(2*i), y0+(2*i), y0+(3*i), y0+(3*i), y0+(4*i), y0+(4*i) ];
% xList = [ x0, x0+i, x0+(2*i), x0,   x0+(2*i), x0,       x0+i,     x0+(2*i), x0,       x0+i,     x0,       x0+(2*i) ];

% R with +i spacing and corners
% (NOTE: In projector bitmap, y corresponds to rows, x corresponds to columns)
y0 = 460; %500
x0 = 385; %370 (for i=8)
i = 45; %40
yList = [ y0, y0,   y0,       y0+i, y0+i,     y0+(2*i), y0+(2*i), y0+(2*i), y0+(3*i), y0+(3*i), y0+(4*i), y0+(4*i) ];
xList = [ x0, x0+i, x0+(2*i), x0,   x0+(2*i), x0,       x0+i,     x0+(2*i), x0,       x0+i,     x0,       x0+(2*i) ];
% Add three points in each of the other corners
y1 = y0-(4*i); %260  
y2 = y0-(4*i);
y3 = y0+(4*i); %540
x1 = x0; %220
x2 = x0+210; %420
x3 = x0+210; %420
yPnts = [ y1, y1,   y1+i, y2, y2,   y2+i, y3, y3-i, y3   ]; 
xPnts = [ x1, x1+i, x1  , x2, x2-i, x2,   x3, x3,   x3-i ];
yList = [yList,yPnts];
xList = [xList,xPnts];

% % Get the range of values that you're projecting into
% minX = min(xList)
% maxX = max(xList)
% minY = min(yList)
% maxY = max(yList)
        
% Evaluate any name-value argument pairs
for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

% Initialize an empty bitMap 
alignBitMap = zeros(h, w, 'logical');

% Convert the locations into indices
pixels = sub2ind(size(alignBitMap), yList, xList);

% Add the activations to the bitmap
alignBitMap(pixels) = true;

% Repeat to get the correct size
alignBitMap = repmat(alignBitMap,1,1,framesPerUp);

% Take the projector bitMap and generate the necessary outMap(s) for projection
[ alignOutMap ] = generateOutmaps( alignBitMap, framesPerUp );

% Save the necessary variables in the correct location
saveAlignmentVariables( alignBitMap, xList, yList );

end

