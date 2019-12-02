function [ xPos, yPos ] = findLocsTD(Usel, Xsel, Ysel, postStimWindow)
% Get x and y positions of at first touchdown post stimulus

% Replicate the indicator function for postStimWindow period
postStimWindow = repmat(postStimWindow,1,6);

% Stores the x-values and y-values
xPos = nan(size(Usel,2), size(Usel,3));
yPos = nan(size(Usel,2), size(Usel,3));

% Loop though each trial
for j = 1:size(Usel,3)
    
    % Get the current trial variables
    u = Usel(:,:,j);
    x = Xsel(:,:,j);
    y = Ysel(:,:,j);
    
    % Compute the difference in logicals and select only touchdown locs
    % NOTE: Need to shift all values by one timepoint
    loc = [zeros(1,6); diff(u,1,1)==1];
    
    % Restrict to only events post stimulus
    loc = loc & postStimWindow;
    
    % Generate a logical for the position of the first post-stim entry
    idx = cumsum(loc,1)==1 & loc;
    
    % Get the locations that have values
    % NOTE: This deals with the case where a touchdown isn't found for any of the limbs
    hasVal = any(idx,1);
    
    % Store the limb positions associated with the current frame
    xPos(hasVal,j) = x(idx);
    yPos(hasVal,j) = y(idx);
    
end

% Convert xPos and yPos into lists of positions mixing all limbs
xPos = xPos(:);
yPos = yPos(:);

% % Test: Sanity check for the locations of each column in final image
% xPos = xPos(1,:);
% yPos = yPos(1,:);

end