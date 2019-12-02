function [ stim ] = rotateBitMap( stim )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Get the bitmap
rawBitMap = stim.cameraBitMap;
n = size(rawBitMap,4);
rawBitMap = repmat(rawBitMap, 1, 1, 1, stim.numRotations);
cameraBitMap = zeros(size(rawBitMap), 'uint8');

% Generate the rotation angles
rotationAngle = 360.*rand(stim.numRotations,1);
rotationIndex = kron(uint32((1:stim.numRotations)'), ones(n,1, 'uint32'));
rotationAngle = rotationAngle(rotationIndex);

% Get the indices to crop
cx = floor(size(cameraBitMap,1)/2)+1;
cy = floor(size(cameraBitMap,2)/2)+1;
idxH = (cx - (stim.height/2)):(cx + (stim.height/2-1));
idxW = (cy - (stim.width/2)):(cy + (stim.width/2-1));

% Rotate each frame
for flipInd=1:size(cameraBitMap,4)
    for frameInd=1:size(cameraBitMap,3)
        rotatedFrame = imrotate(squeeze(rawBitMap(:,:,frameInd,flipInd)), rotationAngle(flipInd), stim.interpMethod, 'crop');
        cameraBitMap(:,:, frameInd, flipInd) = rotatedFrame;
    end
end

% Crop the bitmap to the correct size
cameraBitMap = cameraBitMap(idxH, idxW,:,:);
% Store the output
stim.cameraBitMap = cameraBitMap;
stim.rotationAngle = rotationAngle;

end

