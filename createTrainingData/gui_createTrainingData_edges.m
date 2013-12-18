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

% Last Modified by GUIDE v2.5 18-Dec-2013 10:15:59

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
pos = get(0,'userdata')
% which edge (edgeListInd) is this pixel part of?

% toggle the state of this edge

% update the displayed edge map and the table


% --- Executes on button press in push_loadEdgeMap.
function push_loadEdgeMap_Callback(hObject, eventdata, handles)
% hObject    handle to push_loadEdgeMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.edit_edgemap_path, 'string');
handles.edgeMap = load(str);
load(str); % loads edgeMap (edgeMap as matrix)
handles.edgeMap = edgeMap;
fh = imshow(handles.edgeMap,[]);

dcm_obj = datacursormode(handles.edgeMap);
datacursormode on;
set(dcm_obj,'UpdateFcn', @myupdatefcn )


function txt = myupdatefcn(~, event_obj)
  pos = event_obj.Position;
  disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);
  txt = {'Point to Compute'};
  set(0,'userdata',pos);



function rawImgPath_Callback(hObject, eventdata, handles)
% hObject    handle to rawImgPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rawImgPath as text
%        str2double(get(hObject,'String')) returns contents of rawImgPath as a double
% rawImgPath = get(hObject,'String');




% --- Executes during object creation, after setting all properties.
function rawImgPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rawImgPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_rawImg.
function pushbutton_rawImg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rawImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rawImgPath = get(handles.rawImgPath,'String');
labelImgPath = get(handles.labelImgPath,'String');

% obtain initial set of connected components
[c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,...
    c_extNodeInds,edgeListInds,edgepixels,nodeInds,ws,...
            inactiveEdgeLIDs,offWsIDs,setOfRegions,edges2regions,...
            rawImage,labelImage,strDataVisualization] = ...
    getInitialStructuredLabels(rawImgPath,labelImgPath);

% assign to handles:
handles.edgeListInds = edgeListInds;
handles.edgepixels = edgepixels;
handles.nodeIndsAll = nodeInds;
handles.ws = ws;
handles.edges2regions = edges2regions;
handles.edgeIDSetRegions = setOfRegions;
handles.edges2nodes = edges2nodes;

handles.edgeLID_confVector = edgeLID_confVector;
handles.wsID_confVector = wsID_confVector;
handles.nodeLID_confVector = nodeLID_confVector;


strDataVisualization = ones(100);

handles.componentMat = strDataVisualization;
imageHandle = imshow(handles.componentMat,[]);

% dataPassing = ApplicationData;
% set(imageHandle, 'ButtonDownFcn', @(x,y) mouseClickCapture(x,y,dataPassing));
% set(handles, 'ResizeFcn',     @(x,y) updateXY(x,y, dataPassing));

handles.imageHandle = imageHandle;
handles.setXY_fhandle = @setXYcoordinates;
% handles.getCompProps_fhandle = @getComponentProperties;

set(handles.imageHandle,'ButtonDownFcn',{@ImageClickCallback,handles});


function mouseClickCapture(hAxis, ignored, dataPassingClass)
mousePositionData = get(hAxis,'CurrentPoint');
dataPassingClass.x = mousePositionData(1,1);
dataPassingClass.y = mousePositionData(1,2);


function ImageClickCallback (objectHandle,eventData,handles)
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = floor(coordinates(1,1:2));
% set( gcf, 'UserData', coordinates(1:2) );
x = coordinates(1);
y = coordinates(2);
%setXYcoordinates(x,y,handles);
handles.setXY_fhandle(x,y,handles);
% handles.getCompProps_fhandle(x,y,handles);
% message     = sprintf('x: %d , y: %d',coordinates(1),coordinates(2));
% helpdlg(message);

function setXYcoordinates(x,y,handles)
x_str = int2str(x);
set(handles.edit_xcoordinate,'String',x_str);
y_str = int2str(y);
set(handles.edit_ycoordinate,'String',y_str);



function labelImgPath_Callback(hObject, eventdata, handles)
% hObject    handle to labelImgPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labelImgPath as text
%        str2double(get(hObject,'String')) returns contents of labelImgPath as a double


% --- Executes during object creation, after setting all properties.
function labelImgPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelImgPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function mainFigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mainFigure


% --- Executes on selection change in popupmenu_type.
function popupmenu_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_type


% --- Executes during object creation, after setting all properties.
function popupmenu_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_compState_Callback(hObject, eventdata, handles)
% hObject    handle to edit_compState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_compState as text
%        str2double(get(hObject,'String')) returns contents of edit_compState as a double

handles.compState = str2double(get(hObject,'String'));



% --- Executes during object creation, after setting all properties.
function edit_compState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_compState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_compInd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_compInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_compInd as text
%        str2double(get(hObject,'String')) returns contents of edit_compInd as a double


% --- Executes during object creation, after setting all properties.
function edit_compInd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_compInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xcoordinate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xcoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xcoordinate as text
%        str2double(get(hObject,'String')) returns contents of edit_xcoordinate as a double


% --- Executes during object creation, after setting all properties.
function edit_xcoordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xcoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ycoordinate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ycoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ycoordinate as text
%        str2double(get(hObject,'String')) returns contents of edit_ycoordinate as a double


% --- Executes during object creation, after setting all properties.
function edit_ycoordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ycoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get value set in state (conf)

% update state vector of relevant component type
% type 1: edge
% type 2: node
% type 3: region

% update label vectors

% update image



