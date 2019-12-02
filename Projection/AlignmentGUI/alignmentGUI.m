function varargout = alignmentGUI(varargin)
% ALIGNMENTGUI MATLAB code for alignmentGUI.fig
%      ALIGNMENTGUI, by itself, creates a new ALIGNMENTGUI or raises the existing
%      singleton*.
%
%      H = ALIGNMENTGUI returns the handle to a new ALIGNMENTGUI or the handle to
%      the existing singleton*.
%
%      ALIGNMENTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGNMENTGUI.M with the given input arguments.
%
%      ALIGNMENTGUI('Property','Value',...) creates a new ALIGNMENTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alignmentGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alignmentGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alignmentGUI

% Last Modified by GUIDE v2.5 20-Feb-2017 13:20:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alignmentGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @alignmentGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before alignmentGUI is made visible.
function alignmentGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alignmentGUI (see VARARGIN)

% Choose default command line output for alignmentGUI
handles.output = hObject;

% Define a default directory for searching for alignment files
start_path = 'D:\Data\_FeatureExtractionFiles\ProjectorAlignmentFiles\Alignment_Files';

% Select the directory that contains the projector bitmap, camera image, and projector points
[PathName] = uigetdir(start_path,'Select the directory that contains the alignment files.');
handles.InputPath  = PathName;
cd(handles.InputPath);

% Load the Projector Frame data
load(fullfile(handles.InputPath,'ProjectorFrame.mat'));

% Display the reference frame (projector bitmap) in axes1
axes(handles.axes1);
handles.projectorFrame = ProjectorFrame;
himage = imshow(handles.projectorFrame);

% Display the camera frame in axes2
axes(handles.axes2);
handles.cameraFrame = imread(fullfile(handles.InputPath,'CameraFrame.tif'));
himage = imshow(handles.cameraFrame);

% Set the function that is called with the plot is clicked
set(himage,'ButtonDownFcn',@axes2_ButtonDownFcn);

% Create an array that stores the paired sets of points from the projector and the camera
% NOTE: The pointArray will have the following columns [xProjector, yProjector, xCamera, yCamera]
handles.pointArray = zeros(length(xList),4);
handles.pointArray(:,1) = xList';
handles.pointArray(:,2) = yList';

% Set the initial conditions for the GUI
handles.curIdx = 1;

% Display the current location
axes(handles.axes1);
hold on;
plot(handles.pointArray(handles.curIdx,1),handles.pointArray(handles.curIdx,2),'linestyle','none','marker','o','color','r');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes alignmentGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alignmentGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pass in the current handles structure from guidata
handles = guidata(hObject); 

% Capture mouse clicks on the plot
coordinates = get(handles.axes2,'CurrentPoint'); 
xCamLoc = coordinates(1,1);
yCamLoc = coordinates(1,2);

% Store the location in the labels matrix
handles.pointArray(handles.curIdx,3) = xCamLoc;
handles.pointArray(handles.curIdx,4) = yCamLoc;

% Print a statement describing the current data captured by the GUI
fprintf('(col,row): (%.2f, %.2f) in projector frame --> (%.2f, %.2f) in camera frame \n', ...
    handles.pointArray(handles.curIdx,1),handles.pointArray(handles.curIdx,2), ...
    handles.pointArray(handles.curIdx,3),handles.pointArray(handles.curIdx,4));

% Check if this is the last element in the list
% No --> Continue with labeling
if handles.curIdx ~= size(handles.pointArray,1)
    
% Iterate the counter
handles.curIdx = handles.curIdx + 1;

% Update the GUI
[ handles, himage ] = updateAlignmentGUI( handles );

% Set the image properties so ButtonDownFcn will work
set(himage,'ButtonDownFcn',@axes2_ButtonDownFcn);

% Pass the updated handles structure back to guidata
guidata(handles.axes2,handles);
% guidata(hObject,handles); 

% Yes --> Save and exit the GUI
else
    
    % Extract the pointArray structure from the dataset and save
    coordinates = handles.pointArray;
    save('coordinates.mat','coordinates');
    
    % Calculate the transformation matrix
    [ tform, status ] = calculateCoordinateTransformationMatrix( coordinates );
    
    % Check if this transformation matrix was fit properly
    if status == 0
        save('tform.mat','tform');
    else
        fprintf('An error occurred while calculating the transformation matrix. Please try matching points again.');
    end
    
    % Close the GUI
    close;
    
end
