function [ timestamp ] = projectorWrapper( outMap, framesPerUp, alignPulseOutMap, mapIdx, framesOn, frameOff, projBitMap, cameraBitMap, varargin )
% This function is a wrapper function for performing all the necessary
% steps to control the projector for single leg activation experiments.
% It uses psychtoolbox's OpenGL bindings to project a given bitmap 
% (processed into the required outMap format) using the DLP.

% NOTE: timestamp will be the timestamps of each projector flip

% Set the patternBitDepth based on framesPerUp
switch framesPerUp
    case 3
        patternBitDepth = 7;
    case 6 
        patternBitDepth = 4;
    case 12
        patternBitDepth = 2;
    case 24
        patternBitDepth = 1;
    otherwise
        error('Invalid framesPerUp value.')
end

% Define some defaults
saveResults = true;
startDir = 'D:\';
% startDir = 'C:\Users\labuser\Desktop\Brian_SL_Tests\Tests';

% Evaluate any name-value argument pairs
for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

% If no save path is defined, prompt the user to select one
if saveResults && exist('saveDir') == 0
    saveDir = uigetdir(startDir, 'Select the folder which will store the corresponding video.');
end
    
% Initialize a timestamp array
numFlips = length(mapIdx); % Two additional pulses for alignment
timestamp = nan(numFlips,1);

% Initialize the variable used to return the projector from alignment to projection
alignToProjFlag = false;

% Define an inactive outMap
inactiveOutMap = zeros(size(outMap,1),size(outMap,2),size(outMap,3),'uint8');
% NOTE: Having non-zero alpha values (which shouldn't change the image when
% applying only a single texture fixes timing issues that are otherwise
% present. It's non-sensical but it works.)
inactiveOutMap(:,1:2:end,4) = 100;
inactiveOutMap(:,2:2:end,4) = 200;

% Set a cleanup task for graceful exit
deadManSwitch = onCleanup(@() sca);

% Tell PTB not to display the welcome screen while calibrating
Screen('Preference', 'VisualDebuglevel', 3);

% Initialize the DLP
% DLP = initializeLCR4500('patternBitDepth', patternBitDepth, 'patternColor','red'); %Red
DLP = initializeLCR4500('patternBitDepth', patternBitDepth, 'patternColor','blue'); %Blue

% Open a new PTB window
screens=Screen('Screens');
screenNumber=max(screens);
[windowPtr, ~] = Screen('OpenWindow', screenNumber, 0);

% Set the priority level to max
topPriorityLevel = MaxPriority(windowPtr);
Priority(topPriorityLevel);

% Clear the buffer
[ ~ ] = Screen('Flip', windowPtr);

for cycleIdx = 1:numFlips
    
    
    % Select the correct outMap based on the stimulation list
    if mapIdx(cycleIdx) == 0
        curOutMap = inactiveOutMap;
    
    % Flash the alignment pulse
    elseif mapIdx(cycleIdx) == -1 % Alignment Pulses
        DLP.setPatternAttributes(patternBitDepth, 'blue');
        curOutMap = alignPulseOutMap;
    % Project an opto outmap
    else
        curOutMap = outMap(:,:,:,mapIdx(cycleIdx));
    end
    
    % Convert the image matrix into an OpenGL texture
    tex = Screen('MakeTexture',windowPtr, squeeze(curOutMap), [], 1);
    
    % Draw the texture into the buffer
    Screen('DrawTexture', windowPtr, tex);
    
    % Flip the buffer
    if cycleIdx == 1
        timestamp(cycleIdx) =  Screen('Flip', windowPtr );
    else
        timestamp(cycleIdx) =  Screen('Flip', windowPtr, timestamp(cycleIdx-1)+(1/120));
    end
    
    % Delete the current texture
    Screen('Close',tex);
    
    % Re-initialize the projector to red after a blue alignment pulse
    % NOTE: setPatternAttributes can cause issues with the proper
    % projection of frames for a currently poorly understood reason. It's
    % possible that this is clearing frame buffers on the projector. A
    % work-around for this issue is to pad 'inactive' frames around each
    % alignment pulse. In this case, any dropped presentations will be
    % blank frames and all stimulus and alignment frames will be properly
    % displayed with corresponding timestamps.  
    if mapIdx(cycleIdx) == -1
        % Make a black frame
        pause(0.05);
        tex = Screen('MakeTexture',windowPtr, squeeze(inactiveOutMap), [], 1);
        Screen('DrawTexture', windowPtr, tex);
        Screen('Flip', windowPtr, timestamp(cycleIdx)+(1/120));
        % Reset the DLP to red
        DLP.setPatternAttributes(patternBitDepth, 'red'); % Comment out to align photodiode
    end

%     % It's unclear exactly why this alternative solution does NOT work. 
%     if cycleIdx > 10 && mapIdx(cycleIdx-10) == -1
%         % Reset the DLP to red
%         DLP.setPatternAttributes(patternBitDepth, 'red');
%     end

end

% Ensure everything is closed
Screen('CloseAll');

% Save the projector files associated with this run
if saveResults
    
    save(fullfile(saveDir,'projectorData.mat'), '-v7.3');
    
end

end

