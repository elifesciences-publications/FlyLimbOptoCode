function DLP = initializeLCR4500(varargin)
% This function initializes the DLP LCR4500 to run using psych-toolbox

% NOTE: framesPerUp should have the following relationship to patternBitDepth
% framesPerUp --> patternBitDepth
% 3 --> 7
% 6 --> 4
% 12 --> 2
% 24 -->1

patternBitDepth = 1; 
patternColor = 'red';

for ii = 1:2:length(varargin)-1
    eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

DLP = Lcr4500();
DLP.connect();
pause(0.01);
DLP.wakeup();
pause(0.01);
DLP.setMode(LcrMode.PATTERN);
pause(0.01);
DLP.setPatternAttributes(patternBitDepth, patternColor);

end