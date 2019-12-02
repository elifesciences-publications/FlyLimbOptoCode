function alignmentProjectorWrapper( outMap, framesPerUp, varargin  )
% This function is a wrapper function for controlling the projection of the 
% alignment frame for mapping the projector to the camera pixels. It uses 
% psychtoolbox's OpenGL bindings to project a given bitmap  (processed into 
% the required outMap format) using the DLP.

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

% Set defaults
expDuration = 2700; % seconds

% Evaluate any name-value argument pairs
for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

% Calculate the number of projector flips in the experiment
% NOTE: rounded down in case the expDuration given is not an integer
numFlips = floor(expDuration * 60); % Projector runs at 60Hz

% Set a cleanup task for graceful exit
deadManSwitch = onCleanup(@() sca);

% Tell PTB not to display the welcome screen while calibrating
Screen('Preference', 'VisualDebuglevel', 3);

% Initialize the DLP
initializeLCR4500('patternBitDepth', patternBitDepth,'patternColor','red');

% Open a new PTB window
screens=Screen('Screens');
screenNumber=max(screens);
[windowPtr, ~] = Screen('OpenWindow', screenNumber, 0);

% Clear the buffer
[ ~ ] = Screen('Flip', windowPtr);

% Initialize a timestamp array
timestamp = nan(numFlips,1);

missCount = 0;

for cycleIdx = 1:numFlips
    
    % Convert the image matrix into an OpenGL texture
    tex = Screen('MakeTexture',windowPtr, squeeze(outMap), [], 1);
    
    % Draw the texture into the buffer
    Screen('DrawTexture', windowPtr, tex);
    
    % Flip the buffer
    if cycleIdx == 1
        [timestamp(cycleIdx), ~,~, missed] =  Screen('Flip', windowPtr );
    else 
        [timestamp(cycleIdx), ~,~, missed] =  Screen('Flip', windowPtr, timestamp(cycleIdx-1)+(1/120));
    end
    
    if missed > 0
        missCount = missCount + 1;
    end
    
    % Delete the current texture
    Screen('Close',tex);

end

fprintf('Missed %g percent of frames',100*missCount/numFlips)

% Ensure everything is closed
Screen('CloseAll');

end

