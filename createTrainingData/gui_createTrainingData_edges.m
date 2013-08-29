function varargout = gui_createTrainingData_edges(varargin)
% GUI_CREATETRAININGDATA_EDGES MATLAB code for gui_createTrainingData_edges.fig
%      GUI_CREATETRAININGDATA_EDGES, by itself, creates a new GUI_CREATETRAININGDATA_EDGES or raises the existing
%      singleton*.
%
%      H = GUI_CREATETRAININGDATA_EDGES returns the handle to a new GUI_CREATETRAININGDATA_EDGES or the handle to
%      the existing singleton*.
%
%      GUI_CREATETRAININGDATA_EDGES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CREATETRAININGDATA_EDGES.M with the given input arguments.
%
%      GUI_CREATETRAININGDATA_EDGES('Property','Value',...) creates a new GUI_CREATETRAININGDATA_EDGES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_createTrainingData_edges_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_createTrainingData_edges_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_createTrainingData_edges

% Last Modified by GUIDE v2.5 29-Aug-2013 16:58:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_createTrainingData_edges_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_createTrainingData_edges_OutputFcn, ...
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


% --- Executes just before gui_createTrainingData_edges is made visible.
function gui_createTrainingData_edges_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_createTrainingData_edges (see VARARGIN)

% Choose default command line output for gui_createTrainingData_edges
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_createTrainingData_edges wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_createTrainingData_edges_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_edgemap_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_edgemap_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_edgemap_path as text
%        str2double(get(hObject,'String')) returns contents of edit_edgemap_path as a double


% --- Executes during object creation, after setting all properties.
function edit_edgemap_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_edgemap_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_toggleEdgeState.
function pushbutton_toggleEdgeState_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_toggleEdgeState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the current position (data point)

% which edge (edgeListInd) is this pixel part of?

% toggle the state of this edge

% update the displayed edge map and the table


% --- Executes on button press in push_loadEdgeMap.
function push_loadEdgeMap_Callback(hObject, eventdata, handles)
% hObject    handle to push_loadEdgeMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.edit_edgemap_path, 'string');
load(str); % loads edgeMap (edgeMap as matrix)
handles.edgeMap = edgeMap;
imshow(handles.edgeMap,[]);




function edit_edges2pixels_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_edges2pixels_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_edges2pixels_path as text
%        str2double(get(hObject,'String')) returns contents of edit_edges2pixels_path as a double


% --- Executes during object creation, after setting all properties.
function edit_edges2pixels_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_edges2pixels_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_loadEdgePixelsList.
function push_loadEdgePixelsList_Callback(hObject, eventdata, handles)
% hObject    handle to push_loadEdgePixelsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.edit_edges2pixels_path, 'string');
load(str); % loads edges2pixels
handles.edgeIDs = edges2pixels(:,1);
edges2pixels(:,1) = [];
handles.edgepixelsList = edges2pixels; 



function edit_onEdgeStates_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_onEdgeStates_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_onEdgeStates_path as text
%        str2double(get(hObject,'String')) returns contents of edit_onEdgeStates_path as a double


% --- Executes during object creation, after setting all properties.
function edit_onEdgeStates_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_onEdgeStates_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_loadOnEdgeStateList.
function push_loadOnEdgeStateList_Callback(hObject, eventdata, handles)
% hObject    handle to push_loadOnEdgeStateList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str=get(handles.edit_onEdgeStates_path, 'string');
load(str); % loads onEdgeStates
handles.onEdgeStates = onEdgeStates;
% display onEdgeStates in the uitable
set(handles.uitable_onEdgeStates, 'Data',handles.onEdgeStates);
