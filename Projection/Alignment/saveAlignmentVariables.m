function saveAlignmentVariables( alignBitMap, xList, yList )
% This function saves the necessary alignment information for spatial 
% alignment of the camera and projector. 

% Define the default parent path
parentFolder = 'F:\Data\_FeatureExtractionFiles\ProjectorAlignmentFiles\Alignment_Files';

% Request the folder name
folderName = input('Provide the Date in YYYYMMDD format:','s');

% Create and move to the folder storing the alignment data
folderFullPath = fullfile(parentFolder,folderName);
mkdir(folderFullPath);
cd(folderFullPath);

% Save each of the variables in the appropriate format
ProjectorFrame = alignBitMap(:,:,1);
save('ProjectorFrame.mat','ProjectorFrame','xList','yList');

end

