function varargout = SaveImgs(varargin)
% SaveImgs(GUI)
% A GUI tool for SaveFlxImgs getting the user defined information
%(format,resolution,...) for exporting images created by NetDraw|KeggDraw
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% June 2016

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SaveImgs_OpeningFcn, ...
                   'gui_OutputFcn',  @SaveImgs_OutputFcn, ...
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


% --- Executes just before SaveImgs is made visible.
function SaveImgs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SaveImgs (see VARARGIN)
set(handles.ImgResPopupmenu,'String',{'300','350','400','450',...
        '500','550','600','1200'});
set(handles.ImgFormatPopupmenu,'String',{'PNG','JPEG','BMP',...
        'TIFF','PDF'});
% Choose default command line output for SaveImgs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SaveImgs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SaveImgs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ImgResPopupmenu.
function ImgResPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ImgResPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImgResPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImgResPopupmenu


% --- Executes during object creation, after setting all properties.
function ImgResPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgResPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImgFormatPopupmenu.
function ImgFormatPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ImgFormatPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImgFormatPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImgFormatPopupmenu


% --- Executes during object creation, after setting all properties.
function ImgFormatPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgFormatPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Donepushbutton1.
function Donepushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to Donepushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImgFrmt=get(handles.ImgFormatPopupmenu,'Value');
setappdata(0,'ImgFrmt',ImgFrmt);
ImgRes=get(handles.ImgResPopupmenu,'Value');
setappdata(0,'ImgRes',ImgRes);
close(gcf)
