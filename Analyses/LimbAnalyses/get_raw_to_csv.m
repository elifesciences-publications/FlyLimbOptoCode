% get_raw_to_csv.m: This script extracts out all the variables of needed
% for limb analyses and writes them to a .csv which can be read by python

% Define the parent path
parent_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/';

% Define the save_path
save_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/csv/';

% Define the list of filenames
filenames = {'20190121_R38B08(Bristle)-Chrimson_20171103-20171112', ...
'20190617_Chrimson_Control_20190607_Merged', ...
'20190703_Gr5a-Chrimson_20180314-20190702_Merged(Large)'};

%% Define the variable lists

% Define the list of limb prefixes
limbList = {'L1','L2','L3','R1','R2','R3'};

% Define the the list of body velocity related variables
velVarList = {'angVel_radPerSec','forwardSpeed_mmPerSec',...
              'translationalSpeed_mmPerSec'};

% Define the list of hit variables
hitVarList = strcat(limbList, '_hit');

% Define lists for the limb position and limb phase related variables
limbVarListX = strcat(limbList, '_xPlot_mm');
limbVarListY = strcat(limbList, '_yPlot_mm');

% Define the body variables of interest
comVarList = {'Orient_Rad_FullRotation', 'xCOM','yCOM'};

% Define the list of limb hit location variables
hitLocVarListX = strcat(limbList, '_xHitLoc');
hitLocVarListY = strcat(limbList, '_yHitLoc');

% Define the list of variables that will be included in the csv
varList = ['uniqueFlyTrajID', comVarList, velVarList, hitVarList, ...
    limbVarListX, limbVarListY, hitLocVarListX, hitLocVarListY ];

%% Loop through each of the files in the list and generate a .csv of only 
% the selected variables

for idx = 1:length(filenames)
    
    % Load newData from the .mat file
    load(strcat(parent_path,filenames{idx},'.mat'), 'newData')
    
    % Restrict the dataset to only the variables we want in the .csv
    newData = newData(:, varList);
    
    % Write the new table to disk
    writetable(newData,strcat(save_path,filenames{idx},'.csv'))
    
    % Status update
    fprintf('Done with %d of %d. \n', idx, length(filenames))
    
end
    
    
    
