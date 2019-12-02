function [ mapIdx ] = generateOutMapList( numFlips, mapIdx, framesOn, frameOff, framesPerUp, shiftsPerDim, varargin )
% This function generates a list of outMap indices for a given experiment.
% It does this by calculating the number of flips associated with each
% stimulus event and building a list of indices into dim 4 of outMap for
% presentation

% NOTE: To fill mapIdx we expect the following ordering in the outMaps structure (Along Dim4).
% Indexing Phase by p and Stimulus Duration by t (time) we get a list of
% elements:
% [p1t1, p1t2, ..., p1tN, p1tFull, p2t1, p2t2, ..., p2tN, p2tFull, ..., pNtN, pNtFull]
% where tFull refers to a flip that is completely filled with frames of
% phase p rather than partially filled.

% Pull these from a structure that stores all these variables
numPhases = shiftsPerDim^2;
numTimes = length(framesOn);

% Convert all stimulus presentations into number of flips
flipsOn = ceil(framesOn/framesPerUp);

% Determine the total time needed for each stimulus presentation and following delay
flipsOnOff = flipsOn + ceil(frameOff/framesPerUp);

% Determine how many stimulations are needed to guarantee we fill the experimental duration
numStims = ceil(numFlips/min(flipsOnOff));

% NOTE: Assumes we want a random phase of the lattice presented at each stimulation event
phaseIDs = datasample([1:numPhases],numStims);
timeIDs  = datasample([1:numTimes], numStims);

count = 1;
curIdx = 1;

% Loop over the stimulation durations and offsets to the desired list of map indices
while curIdx  <= (numFlips-max(flipsOnOff)+1) % change what you are checking here to look at the last index of a stimulus rather than the first

    % Determine the length of the current stimulus in number of flips
    flips = flipsOnOff(timeIDs(count)); 
    
    % Determine what indicies must fill flips
    % NOTE: index of zero will reference the inactive outMap
    outMapList = zeros(flips,1);
    lastStimFlip = flipsOn(timeIDs(count));
    outMapList(1:lastStimFlip-1) = phaseIDs(count)*(numTimes+1);
    outMapList(lastStimFlip) = (phaseIDs(count)-1)*(numTimes+1) + timeIDs(count);
    
    % For each index in phaseIDs and timeIDs, populate the mapIdx array
    outMapIdx = (phaseIDs(count)-1)*numTimes + timeIDs(count); 
    
    % Populate the mapIdx array with the correct indicies for this stimulation
    mapIdx(curIdx:curIdx+flips-1) = outMapList;
    
    % Increment the counter
    count = count+1;
    
    % Update curIdx
    curIdx = curIdx + flips;
    
end

end

