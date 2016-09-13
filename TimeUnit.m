function varargout = TimeUnit(varargin)
% TimeUnit(GUI)
% is prompted for each visualization task in KeggDraw
% to determine the time-serie unit of input data (in dynamic conditions).

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Mar. 2015

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TimeUnit_OpeningFcn, ...
                   'gui_OutputFcn',  @TimeUnit_OutputFcn, ...
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


% --- Executes just before TimeUnit is made visible.
function TimeUnit_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = TimeUnit_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function InsertManuallyEdit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function InsertManuallyEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in ChooseUipanel.
function ChooseUipanel_SelectionChangeFcn(hObject, eventdata, handles)

% --- Executes on button press in OKPushbutton.
function OKPushbutton_Callback(hObject, eventdata, handles)
IMedit=get(handles.InsertManuallyEdit, 'String');
if strcmp(IMedit,'Or Insert manually')
    hSelectedObj=get(handles.ChooseUipanel, 'SelectedObject');
    selectedObjTag=get(hSelectedObj, 'Tag');
    switch selectedObjTag
        case 'DonotshowRadio'
            SelectString='NShow';
        case 'SecondRadio'
            SelectString='s';
        case 'MinuteRadio'
            SelectString='m';
        case 'HourRadio'
            SelectString='h';
        case 'NoRadio'
            SelectString='NO';
    end
else
    SelectString=IMedit;
end
setappdata(0,'SelectString',SelectString)
close(TimeUnit)
