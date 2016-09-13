function varargout = Bigg2KeggTable(varargin)
% BIGG2KEGGTABLE MATLAB code for Bigg2KeggTable.fig
%      BIGG2KEGGTABLE, by itself, creates a new BIGG2KEGGTABLE or raises the existing
%      singleton*.
%
%      H = BIGG2KEGGTABLE returns the handle to a new BIGG2KEGGTABLE or the handle to
%      the existing singleton*.
%
%      BIGG2KEGGTABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIGG2KEGGTABLE.M with the given input arguments.
%
%      BIGG2KEGGTABLE('Property','Value',...) creates a new BIGG2KEGGTABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bigg2KeggTable_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bigg2KeggTable_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bigg2KeggTable

% Last Modified by GUIDE v2.5 17-Mar-2016 11:43:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bigg2KeggTable_OpeningFcn, ...
                   'gui_OutputFcn',  @Bigg2KeggTable_OutputFcn, ...
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


% --- Executes just before Bigg2KeggTable is made visible.
function Bigg2KeggTable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bigg2KeggTable (see VARARGIN)

% Choose default command line output for Bigg2KeggTable
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Bigg2KeggTable wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Bigg2KeggTable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
A=get(handles.uitable1,'Data');
setappdata(0,'Finalvalidity',A);
close(gcf)


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
Diffdata = getappdata(0,'Diffdata');
set(hObject, 'Data', Diffdata);
set(hObject, 'ColumnName', {'BiGG', 'KEGG', 'Validity'})
set(hObject, 'ColumnEditable', [false, false, true]);
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
RowChanged=eventdata.Indices(1);
ColChanged=eventdata.Indices(2);
newValue=eventdata.NewData;
D =get(handles.uitable1,'Data');
D{RowChanged,ColChanged}=newValue;
