function [ projectorBitMap, cameraBitMap ] = generateProjectorBitmap_4D( cameraBitMap, tform, h, w )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Initialize a projector bitMap of the appropriate size
projectorBitMap = zeros(h,w, size(cameraBitMap,3), size(cameraBitMap,4), 'uint8');

for ind = 1:size(projectorBitMap,4)
    curBitMap = squeeze(cameraBitMap(:,:,:,ind));
    
    % Find all the indices of the non-zero elements in curBitMap
    logicals = find(curBitMap ~= 0);
    values = curBitMap(logicals);
    [cam1,cam2,cam3] = ind2sub(size(curBitMap),logicals);
    
    % Calculate the projector coordinates for each of the non-zero elements in the camera bitmap
    
    % % Built-in Affine Method
    % projXY1 = [camX,camY,ones(length(camY),1)] * tform.T; % For built in affine transformation
    
    % Linear regression method with square terms
    projXY1 = [cam2.^2, cam1.^2, cam2.*cam1, cam2, cam1, ones(length(cam1),1)] * tform;
    
    proj1 = projXY1(:,1);
    proj2 = projXY1(:,2);
    proj3 = cam3;
    % projYX1 = [camY,camX,ones(length(camY),1)] * tform.T;
    % projY = projYX1(:,1);
    % projX = projYX1(:,2);
    % projZ = camZ;
    
    
    % Eliminate elements that are outside the bounds of the projector bitmap
    % Remove Invalid Row Values
    rRemove = find(proj2<1 | proj2>h);
    proj2(rRemove)   = [];
    proj1(rRemove)   = [];
    proj3(rRemove)   = [];
    values(rRemove)  = [];
    
    % Check if activations in the camera frame were removed and modify the camera bitmap
    if ~isempty(rRemove)
        
        % Display a message saying that the curBitMap changed
        fprintf('Certain activations could not be generated with the projector alignment.\nRegenerating curBitMap.\n');
        
        % Remove unprojectable regions from the curBitMap
        idx = sub2ind(size(curBitMap), cam1(rRemove), cam2(rRemove), cam3(rRemove));
        curBitMap(idx) = 0;
        
        % Remove the entries from the list of non-zero positions in curBitMap to align indexing
        cam1(rRemove) = [];
        cam2(rRemove) = [];
        cam3(rRemove) = [];
    end
    
    % Remove Invalid Column Values
    cRemove = find(proj1<1 | proj1>w);
    proj2(cRemove)   = [];
    proj1(cRemove)   = [];
    proj3(cRemove)   = [];
    values(cRemove)  = [];
    
    % Check if activations in the camera frame were removed and modify the camera bitmap
    if ~isempty(cRemove)
        
        % Display a message saying that the curBitMap changed
        fprintf('Certain activations could not be generated with the projector alignment.\nRegenerating curBitMap.\n');
        
        % Remove unprojectable regions from the curBitMap
        idx = sub2ind(size(curBitMap), cam1(cRemove), cam2(cRemove), cam3(cRemove));
        curBitMap(idx) = 0;
        
        % Remove the entries from the list of non-zero positions in curBitMap to align indexing
        cam1(cRemove) = [];
        cam2(cRemove) = [];
        cam3(cRemove) = [];
    end
    
    % Put the non-zero elements from curBitMap into the new appropriate locations in projBitMap
    % NOTE: The projected coordinate value are not guaranteed to be integers so
    % we will need to round to integer pixels
    idx = sub2ind(size(projectorBitMap), round(proj2), round(proj1), round(proj3));
    
    curProjBitMap = zeros(h,w, size(cameraBitMap,3), 'uint8');
    curProjBitMap(idx) = values;
    
    projectorBitMap(:,:,:,ind) = curProjBitMap;
    cameraBitMap(:,:,:,ind) = curBitMap;
end

end

