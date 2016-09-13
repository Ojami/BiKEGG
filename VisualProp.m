function varargout = VisualProp(varargin)
% VisualProp(GUI)
% Customizing color properties for KeggDraw maps.
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Mar. 2015

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VisualProp_OpeningFcn, ...
                   'gui_OutputFcn',  @VisualProp_OutputFcn, ...
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


% --- Executes just before VisualProp is made visible.
function VisualProp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualProp (see VARARGIN)

% Choose default command line output for VisualProp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VisualProp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualProp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in TextColPush.
function TextColPush_Callback(hObject, eventdata, handles)
% hObject    handle to TextColPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty (getappdata(0,'TextCol'))
    rmappdata (0,'TextCol')
else
    waitfor(msgbox({'                      Choosing text color is optional';...
        'You can use default properties simply by closing color panel'},...
        'warning','warn'));
end
TextCol=uisetcolor;
if size(TextCol,2) == 3
    setappdata(0,'TextCol',TextCol)
end

% --- Executes on button press in TextFontPush.
function TextFontPush_Callback(hObject, eventdata, handles)
% hObject    handle to TextFontPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty (getappdata(0,'TextFont'))
    rmappdata (0,'TextFont')
else
    waitfor(msgbox({'                      Choosing text font is optional';...
        'You can use default properties simply by closing text panel'},...
        'warning','warn'));
end
TextFont=uisetfont;
if isstruct(TextFont);
    setappdata(0,'TextFont',TextFont)
end

% --- Executes on button press in BoxColPush.
function BoxColPush_Callback(hObject, eventdata, handles)
% hObject    handle to BoxColPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty (getappdata(0,'BoxtCol'))
    rmappdata (0,'BoxtCol')
else
    waitfor(msgbox({'                      Choosing box color is optional';...
        'You can use default properties simply by closing color panel'},...
        'warning','warn'));
end
BoxCol=uisetcolor;
if size(BoxCol,2) == 3
    setappdata(0,'BoxtCol',BoxCol)
end


% --- Executes on button press in Donepushbutton5.
function Donepushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to Donepushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
