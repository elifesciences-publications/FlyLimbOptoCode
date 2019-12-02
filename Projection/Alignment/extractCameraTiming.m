function [ timeStamps, gpioPixels, meanPixels ] = extractCameraTiming()
% This function takes as input a path to the uncompressed video and
% extracts the timing information contained in the underlying pixel values

% NOTE: Use this when running outside of the wrapper function
% Pick the location of the videofile
[videoFile, pathName] = uigetfile('C:\Users\ClarkLab\Desktop\SingleLegActivation\*.*');
videoPath = fullfile(pathName,videoFile);

% Import the AVI video
im = VideoReader(videoPath); %RGB vs. Grayscale determination is read-only in VideoReader.  Need to figure out why it defaults to RGB

% Rename the video properties for convenience.
nFrames = int64(im.Duration*im.FrameRate);
vidHeight = im.Height;
vidWidth = im.Width;

% Initialize the variables for storing various stages of the time data
temp = uint8(0);
timePixels(nFrames, 4) = temp;
gpioPixels(nFrames, 4) = temp;
timeStamps = zeros(nFrames,1);
meanPixels(nFrames,1) = temp;

% Initialize variables for converting from binary to timeStamps
nBits = 8;
secsBits = uint32(hex2dec('7F'));
cycleCountBits = uint32(hex2dec('1FFF'));
cycleOffsetBits = uint32(hex2dec('FFF'));

% Initialize a temporary variable for the movie frames
mov = struct('cdata', zeros(vidHeight, vidWidth, 'uint8'),'colormap', []);

n = 1;
% Read one frame at a time.
while hasFrame(im)
    
    % Read in the current frame
    mov.cdata = readFrame(im);
    
    % % % For Testing of Temporal Alignment % % %
    meanPixels(n) = mean(mov.cdata(:)); 
    % % % For Testing of Temporal Alignment % % %
    
    % Extract the first n pixels that contain the timing information
    timePixels(n,:) = uint8(mov.cdata(1,1:4,1));
    
    % Extract the GPIO information from the video
    gpioPixels(n,:) = uint8(mov.cdata(1,5:8,1));
    
    % Convert the pixel intensities to binary
    str = dec2bin(timePixels(n,:),nBits);
    str = reshape(str',1,32);
    timeBinary = uint32(bin2dec(str));

    % Convert the binary into a timestamps    
    nSecond       = bitand(bitshift(timeBinary,-25, 'uint32'), secsBits);
    nCycleCount   = bitand(bitshift(timeBinary,-12, 'uint32'), cycleCountBits);
    nCycleOffset  = bitand(bitshift(timeBinary,-0, 'uint32'),  cycleOffsetBits);
    timeStamps(n) = double(nSecond) + ((double(nCycleCount)+(double(nCycleOffset)/3072.0))/8000.0);
    
    % Iterate the counter
    n = n + 1;  
        
end

% % NOTE: Corresponding Example C++ code
% imageTimeStampToSeconds(unsigned int uiRawTimestamp)
% {
% 
%    int nSecond      = (uiRawTimestamp >> 25) & 0x7F;   // get rid of cycle_* - keep 7 bits
%    int nCycleCount  = (uiRawTimestamp >> 12) & 0x1FFF; // get rid of offset
%    int nCycleOffset = (uiRawTimestamp >>  0) & 0xFFF;  // get rid of *_count
% 
%    return (double)nSecond + (((double)nCycleCount+((double)nCycleOffset/3072.0))/8000.0);
% }

end

