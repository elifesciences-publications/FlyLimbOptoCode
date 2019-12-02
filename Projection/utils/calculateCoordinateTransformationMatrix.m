function [ tform, status ] = calculateCoordinateTransformationMatrix( coordinates )
% This function takes paired lists of coordinates from the camera and
% projector and generates a matrix that maps the camera coordinates to the
% projector coordinates. This matrix will be used to generate the projector
% bitmaps used for singleLegActivations

% Extract the camera and projector coordinates from the coordinate matrix
camXYPnts  = coordinates(:,3:4);
projXYPnts = coordinates(:,1:2);

% NOTE: camYXPnts and projYXPnts are [n 2] matrices of coordinates

% %% Built in Affine Transformation Calculation
% % Calculate the transformation matrix to convert from camera pixels to projector pixels
% % NOTE: status = 0 means that the fitting ran without errors
% [tform,~,~,status] = estimateGeometricTransform(camXYPnts,projXYPnts,'affine','MaxNumTrials',10000,'MaxDistance',.001);
% % % NOTE: Use this version to check if the projected pixel locations are correct
% % [tform,newCam,newProj,status] = estimateGeometricTransform(camXYPnts,projXYPnts,'projective');

%% Linear Regression Method
% Calculate the transformation matrix to convert from camera pixels to projector pixels
% NOTE: Check that the appropriate call is made in generateProjectorBitmap.m when transforming points

% Write a default value for the status output
status = 0;

% Define the input features
camX = camXYPnts(:,1);
camY = camXYPnts(:,2);
X = [camX.^2, camY.^2, camX.*camY, camX, camY, ones(length(camY),1)];

% Define the output features
y = [projXYPnts, ones(size(projXYPnts,1),1)]; 

% Fit the transformation. tForm will be a 3x6 matrix with last column = [0;0;0;0;0;1]
tform = X\y;


end

