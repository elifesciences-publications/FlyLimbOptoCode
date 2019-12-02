function [ outMap ] = generateOutmapsFromSingleFrameBitmaps( bitMap, framesPerUp, framesOn, frameOff )
% This function generates an outmap structure identical to
% generateOutmaps.m but uses a bitmap structure that includes only single
% frames for each stimulus condition.

% NOTE: To fill the indexing used to present outmaps to the projector, we 
% expect the following ordering in the outMaps structure (Along Dim4).
% Indexing Phase by p and Stimulus Duration by t (time) we get a list of
% elements:
% [p1t1, p1t2, ..., p1tN, p1tFull, p2t1, p2t2, ..., p2tN, p2tFull, ..., pNtN, pNtFull]
% where tFull refers to a flip that is completely filled with frames of
% phase p rather than partially filled.

tic;

% Define the number of distinct outMaps created
% (#Durations+1) * #Phases; +1 is for generating a fully filled framesPerUp of the stimulus
numMaps = (length(framesOn)+1)*(size(bitMap,3));

% Initialize the outMap variable
outMap = zeros(size(bitMap,1),size(bitMap,2),4,numMaps,'uint8');

% Calculate the number of remainder framesOn for each temporal duration
remainFramesOn = mod(framesOn,framesPerUp);

% Add an full frame entry for each duration 
remainFramesOn = [remainFramesOn,framesPerUp];

counter = 1;
% Loop over the number of phases
for phaseIdx = 1:size(bitMap,3)
    
    % Loop over the number of durations
    for timeIdx = 1:length(remainFramesOn)
    
        % Initialize an empty bitmap corresponding to the current outmap
        curBitMap = zeros(size(bitMap,1),size(bitMap,2),framesPerUp,'uint8');
        
        % Populate the current bitMap with the appropriate frames
        curBitMap(:,:,1:remainFramesOn(timeIdx)) = repmat(bitMap(:,:,phaseIdx),1,1,remainFramesOn(timeIdx));
        
        % Use the current bitmap to generate the appropriate outMap
        outMap(:,:,:,counter) = createOutmap(curBitMap,framesPerUp);
        
        % Iterate the counter
        counter = counter+1;
        
    end
end

fprintf('%d outMaps generated in %f seconds.\n',numMaps,toc);

end

