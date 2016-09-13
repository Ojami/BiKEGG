function varargout = KeggDrawTable(varargin)
% KeggDrawTable(GUI)
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Nov. 2015

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KeggDrawTable_OpeningFcn, ...
                   'gui_OutputFcn',  @KeggDrawTable_OutputFcn, ...
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


% --- Executes just before KeggDrawTable is made visible.
function KeggDrawTable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KeggDrawTable (see VARARGIN)
CurrentMap = getappdata(0,'CurrentMap');
CurrentMap1 = ['Map: ',CurrentMap];
set(handles.text2,'String',CurrentMap1)
RptRxns = getappdata(0,'RptRxns');
ListData = {0};
for ct = 1:numel(RptRxns)
    ListData{ct} = RptRxns{ct}{1};
    data1{ct} = 'Not slected';
end
set(handles.KeggListbox,'String',ListData);
set(handles.ResultListbox,'String',data1);
% Choose default command line output for KeggDrawTable
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KeggDrawTable wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KeggDrawTable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data1 = getappdata(0,'data1');
setappdata(0,'data1',data1);
close(gcf)

% --- Executes on selection change in KeggListbox.
function KeggListbox_Callback(hObject, eventdata, handles)
% hObject    handle to KeggListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RptFlx = getappdata(0,'RptFlx');
RptBiggs = getappdata(0,'RptBiggs');
ListChoice = get(handles.KeggListbox, 'Value');
ColNm = size(RptFlx{ListChoice},2);
for ct=1:ColNm
    ColNm1{ct} = num2str(ct);
end
columnname = {'check' , ColNm1{:}};
CurrentSize = size(RptFlx{ListChoice},1);
checkBox = num2cell(false(CurrentSize,1));
columnformat = {'logical','numeric'};
columneditable = true;
Y = num2cell(RptFlx{ListChoice});
MixedData = [checkBox Y];
set(handles.Rxnuitable,'Data',MixedData)
set(handles.Rxnuitable,'ColumnFormat',columnformat)
set(handles.Rxnuitable,'ColumnName',columnname)
set(handles.Rxnuitable,'ColumnEditable',columneditable)
set(handles.Rxnuitable, 'RowName', RptBiggs{ListChoice})

% URptRxns = unique(RptRxns

% Hints: contents = cellstr(get(hObject,'String')) returns KeggListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from KeggListbox


% --- Executes during object creation, after setting all properties.
function KeggListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KeggListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in Rxnuitable.
function Rxnuitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Rxnuitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
data=get(hObject,'Data'); % get the data cell array of the table
cols=get(hObject,'ColumnFormat'); % get the column formats
if strcmp(cols(eventdata.Indices(2)),'logical') % if the column of the edited cell is logical
    if eventdata.EditData % if the checkbox was set to true
        data{eventdata.Indices(1),eventdata.Indices(2)}=true; % set the data value to true
    else % if the checkbox was set to false
        data{eventdata.Indices(1),eventdata.Indices(2)}=false; % set the data value to false
    end
end
set(hObject,'Data',data); % now set the table's data to the updated data cell array

RptBiggs = getappdata(0,'RptBiggs');
KeggChoice = get(handles.KeggListbox, 'Value');
BiggSelect = data(:,1); BiggSelect = [BiggSelect{:}];
BiggInd = find(BiggSelect);
data1 = get(handles.ResultListbox,'String');
if numel(BiggInd) ==1    
    data1{KeggChoice} = RptBiggs{KeggChoice}{BiggInd};
    set(handles.ResultListbox,'String',data1)
    NoteStr = ['BiGG ID: ',data1{KeggChoice},' is selected'];
    set(handles.text1,'String',NoteStr)
    setappdata(0,'data1',data1)
elseif numel(BiggInd)>1
    NoteStr = 'Warning: More than one BiGG ID is selected!';
    set(handles.text1,'String',NoteStr)
else
    data1{KeggChoice} = 'Not slected';
    set(handles.ResultListbox,'String',data1)
    NoteStr = 'Warning: No BiGG ID is selected!';
    set(handles.text1,'String',NoteStr)
    setappdata(0,'data1',data1)
end
    

% --- Executes on selection change in ResultListbox.
function ResultListbox_Callback(hObject, eventdata, handles)
% hObject    handle to ResultListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ResultListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ResultListbox


% --- Executes during object creation, after setting all properties.
function ResultListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ResultListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
