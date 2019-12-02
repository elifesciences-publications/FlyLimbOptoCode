function [ alignPulseOutMap ] = getAlignPulseOutMap( varargin )
% This function generates an outMap for the temporal alignment pulse for
% the projector and camera

% 604 --> Projector Height in Pixels
% 684 --> Projector Width in Pixels
% 512 --> Camera Height in Pixels
% 640 --> Camera Width in Pixels

% Define the alignPulseOutMap size
h = 604; % px 
w = 684; % px
framesPerUp = 24;

% Generate an alignment bitMap
alignBitMap = ones(h,w,framesPerUp,'uint8');

% Use the current bitmap to generate the appropriate outMap
alignPulseOutMap = createOutmap(alignBitMap,framesPerUp);


end

