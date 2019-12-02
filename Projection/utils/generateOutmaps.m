function [ outMap ] = generateOutmaps( bitMap, framesPerUp )
% This function takes bitMaps produced by generateSquareLatticeBitmaps and
% converts them to outMaps used by projectorWrapper for presenting stimuli
% to the flies

tic;

% Define the number of distinct outMaps created
numMaps = size(bitMap,3)/framesPerUp;

% Initialize the outMap variable
outMap = zeros(size(bitMap,1),size(bitMap,2),4,numMaps,'uint8');

% Loop over numMaps breaking the bitMap into framesPerUp sized chunks
for n = 1:numMaps
    
    % Define the correct indices for this iteration
    startIdx = (n-1)*framesPerUp + 1;
    stopIdx  = n*framesPerUp;
    
    % Grab the correct region of the bitMap
    curBitMap = uint8(bitMap(:,:,startIdx:stopIdx));
    
    % Populate the outmap
    outMap(:,:,:,n) = createOutmap(curBitMap,framesPerUp);
    
end

fprintf('%d outMaps generated in %f seconds.\n',numMaps,toc);

end

