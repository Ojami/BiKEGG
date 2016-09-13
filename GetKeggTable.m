function varargout = GetKeggTable(varargin)
% GetKeggTable (GUI)
% Identifies missed (unknown) reaction correspondences for the input BiGG model.
% Part of GetKEGG
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Nov. 2015

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GetKeggTable_OpeningFcn, ...
                   'gui_OutputFcn',  @GetKeggTable_OutputFcn, ...
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


% --- Executes just before GetKeggTable is made visible.
function GetKeggTable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GetKeggTable (see VARARGIN)

% Choose default command line output for GetKeggTable
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GetKeggTable wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GetKeggTable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A=get(handles.Table1,'Data');
setappdata(0,'Final',A);
close(gcf)

% --- Executes during object creation, after setting all properties.
function Table1_CreateFcn(hObject, eventdata, handles)
FirstClm = getappdata(0,'NotData');
SecClm = cell(1,numel(FirstClm));
SecClm(1:end) = {'Unknown'};
TabData = cell(numel(FirstClm),2);
TabData(1:end,1) = FirstClm; TabData(1:end,2) = SecClm;
set(hObject, 'Data', TabData);
set(hObject, 'ColumnName', {'BiGG', 'KEGG'})
set(hObject, 'ColumnEditable', [true, true]);
% hObject    handle to Table1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in Table1.
function Table1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Table1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
RowChanged=eventdata.Indices(1);
ColChanged=eventdata.Indices(2);
newValue=eventdata.NewData;
D =get(handles.Table1,'Data');
D{RowChanged,ColChanged}=newValue;
