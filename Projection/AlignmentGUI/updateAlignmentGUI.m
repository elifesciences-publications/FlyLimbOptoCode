function [ handles, himage ] = updateAlignmentGUI( handles )
% This function updates the plots in the alignmentGUI

%% Update the ProjectorFrame
axes(handles.axes1);
hold off;
imshow(handles.projectorFrame);
hold on;

% Loop through all the pixels up to the current frame (Plot as Dots)
for i = 1:handles.curIdx-1 
    plot(handles.pointArray(i,1),handles.pointArray(i,2),'linestyle','none','marker','.','color','g');
end

% Plot the current frame as a circle
plot(handles.pointArray(handles.curIdx,1),handles.pointArray(handles.curIdx,2),'linestyle','none','marker','o','color','r');

%% Update the CameraFrame
axes(handles.axes2);
hold off;
himage = imshow(handles.cameraFrame);
hold on;

% Loop through all the pixels up to the current frame (Plot as Dots)
for i = 1:handles.curIdx-1 
    plot(handles.pointArray(i,3),handles.pointArray(i,4),'linestyle','none','marker','.','color','g');
end

end

