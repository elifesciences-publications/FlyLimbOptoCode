function [ projectorBitMap, cameraBitMap ] = generateProjectorBitmap_4D_dense( cameraBitMap, tform, h, w )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[cam1,cam2,cam3] = ind2sub([size(cameraBitMap,1), size(cameraBitMap,2), size(cameraBitMap,3)], (1:numel(squeeze(cameraBitMap(:,:,:,1))))' );

% Calculate the projector coordinates for each of the non-zero elements in the camera bitmap

% % Built-in Affine Method
% projXY1 = [camX,camY,ones(length(camY),1)] * tform.T; % For built in affine transformation

% Linear regression method with square terms
projXY1 = [cam2.^2, cam1.^2, cam2.*cam1, cam2, cam1, ones(length(cam1),1)] * tform;

proj1 = uint16(round(projXY1(:,1)));
proj2 = uint16(round(projXY1(:,2)));
proj3 = cam3;

% Eliminate elements that are outside the bounds of the projector bitmap
% Remove Invalid Row Values
outOfBounds = (proj2<1 | proj2>h | proj1<1 | proj1>w);
proj2(outOfBounds)   = [];
proj1(outOfBounds)   = [];
proj3(outOfBounds)   = [];

% Remove the entries from the list of non-zero positions in curBitMap to align indexing
cam1(outOfBounds) = [];
cam2(outOfBounds) = [];
cam3(outOfBounds) = [];

% Initialize a projector bitMap of the appropriate size
projectorBitMap = zeros(h,w, size(cameraBitMap,3), size(cameraBitMap,4), 'uint8');

% Get the appropriate one-dimensional index
idx = sub2ind(size(projectorBitMap), proj2, proj1, proj3);

for ind = 1:size(projectorBitMap,4)
    curBitMap = squeeze(cameraBitMap(:,:,:,ind));
    
    curProjBitMap = zeros(h,w, size(cameraBitMap,3), 'uint8');
    curProjBitMap(idx) = curBitMap(~outOfBounds);
    projectorBitMap(:,:,:,ind) = curProjBitMap;
    cameraBitMap(:,:,:,ind) = curBitMap;
end

end

