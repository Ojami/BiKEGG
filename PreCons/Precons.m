function varargout = Precons(varargin)
% Precons(GUI)
% A GUI for PreconsMethod
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Match 2016

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Precons_OpeningFcn, ...
                   'gui_OutputFcn',  @Precons_OutputFcn, ...
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


% --- Executes just before Precons is made visible.
function Precons_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Precons (see VARARGIN)
Fileid = fopen('orgs.txt','r');
orgcell = textscan(Fileid,'%s %s %[^\n]');
fclose(Fileid);
Tcode = orgcell{1}; Scode = orgcell{2};
for ct = 1:numel(Tcode)
    tscode{ct} = [Scode{ct},'  ',Tcode{ct}];
end
set(handles.OrganPopup,'String',tscode)

% Choose default command line output for Precons
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Precons wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Precons_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in OrganPopup.
function OrganPopup_Callback(hObject, eventdata, handles)
% hObject    handle to OrganPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fileid = fopen('orgs.txt','r');
orgcell = textscan(Fileid,'%s %s %[^\n]');
fclose(Fileid);
Ncode = orgcell{3}; Scode = orgcell{2};
Cursel = get(hObject,'Value');
if ~isempty(Cursel)
    set(handles.Orgtext,'String',Ncode{Cursel})
    set(handles.OrganEdit,'String',Scode{Cursel})
end

stat = get(handles.Biographcheck,'Value');
if stat
    Org = get(handles.OrganEdit,'String');
    if ~strcmp(Org,'Organism code')
        Txt1 = 'Downloding pathway data...';
        fprintf(Txt1)
        [getpathways,stat] = urlread(['http://rest.kegg.jp/list/pathway/',Org]);
        if ~stat
            msgbox('There is a problem with the internet connection!');
            for i =1:numel(Txt1)
                fprintf('\b')
            end
            return
        end
        for i =1:numel(Txt1)
            fprintf('\b')
        end
        pths = textscan(getpathways,'%s %[^\n]');
        pathnames = pths{1}; pathnames = strrep(pathnames,'path:','');
        pathdescription = pths{2};
        set(handles.Pathwaylistbox,'String',strcat(pathnames,'_',pathdescription))
    end
end

% --- Executes during object creation, after setting all properties.
function OrganPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OrganPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OrganEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OrganEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fileid = fopen('orgs.txt','r');
orgcell = textscan(Fileid,'%s %s %[^\n]');
fclose(Fileid);
Tcode = orgcell{1}; Scode = orgcell{2}; Ncode = orgcell{3};
OrgEditSelect = get(hObject,'String');
WhereOrgPop1 = find(ismember(Tcode,OrgEditSelect));
WhereOrgPop2 = find(ismember(Scode,OrgEditSelect));
if ~isempty(WhereOrgPop1)
    set(handles.OrganPopup,'Value',WhereOrgPop1)
    set(handles.Orgtext,'String',Ncode{WhereOrgPop1})
elseif ~isempty(WhereOrgPop2)
    set(handles.OrganPopup,'Value',WhereOrgPop2)
    set(handles.Orgtext,'String',Ncode{WhereOrgPop2})
else
    set(handles.Orgtext,'String','Wrong organism code')
end
stat = get(handles.Biographcheck,'Value');
if stat
    Org = get(handles.OrganEdit,'String');
    if ~strcmp(Org,'Organism code')
        Txt1 = 'Downloding pathway data...';
        fprintf(Txt1)
        [getpathways,stat] = urlread(['http://rest.kegg.jp/list/pathway/',Org]);
        if ~stat
            msgbox('There is a problem with the internet connection!');
            for i =1:numel(Txt1)
                fprintf('\b')
            end
            return
        end
        for i =1:numel(Txt1)
            fprintf('\b')
        end
        pths = textscan(getpathways,'%s %[^\n]');
        pathnames = pths{1}; pathnames = strrep(pathnames,'path:','');
        pathdescription = pths{2};
        set(handles.Pathwaylistbox,'String',strcat(pathnames,'_',pathdescription))
    end
end

% --- Executes during object creation, after setting all properties.
function OrganEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OrganEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DonePush.
function DonePush_Callback(hObject, eventdata, handles)
% hObject    handle to DonePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PreconsMethod(handles)


% --- Executes on selection change in Pathwaylistbox.
function Pathwaylistbox_Callback(hObject, eventdata, handles)
% hObject    handle to Pathwaylistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Pathwaylistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Pathwaylistbox


% --- Executes during object creation, after setting all properties.
function Pathwaylistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pathwaylistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Biographcheck.
function Biographcheck_Callback(hObject, eventdata, handles)
% hObject    handle to Biographcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Biographcheck
stat = get(hObject,'Value');
if stat
    set(handles.Pathwaylistbox,'Enable','on')
    set(handles.G2Bradiobutton,'Enable','on')
    set(handles.G2Kradiobutton,'Enable','on')
    set(handles.G2K2Bradiobutton,'Enable','on')
    Org = get(handles.OrganEdit,'String');
    if ~strcmp(Org,'Organism code')
        Txt1 = 'Downloding pathway data...';
        fprintf(Txt1)
        [getpathways,stat] = urlread(['http://rest.kegg.jp/list/pathway/',Org]);
        if ~stat
            msgbox('There is a problem with the internet connection!');
            for i =1:numel(Txt1)
                fprintf('\b')
            end
            return
        end
        for i =1:numel(Txt1)
            fprintf('\b')
        end
        pths = textscan(getpathways,'%s %[^\n]');
        pathnames = pths{1}; pathnames = strrep(pathnames,'path:','');
        pathdescription = pths{2};
        set(handles.Pathwaylistbox,'String',strcat(pathnames,'_',pathdescription))
    end
else
    set(handles.Pathwaylistbox,'Enable','off')
    set(handles.G2Bradiobutton,'Enable','off')
    set(handles.G2Kradiobutton,'Enable','off')
    set(handles.G2K2Bradiobutton,'Enable','off')
end
