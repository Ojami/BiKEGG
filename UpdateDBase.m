function varargout = UpdateDBase(varargin)
% UpdateDBase(GUI)
% Updating the local database of BiKEGG for offline use.
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Mar. 2015

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UpdateDBase_OpeningFcn, ...
                   'gui_OutputFcn',  @UpdateDBase_OutputFcn, ...
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


% --- Executes just before UpdateDBase is made visible.
function UpdateDBase_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UpdateDBase (see VARARGIN)

% Choose default command line output for UpdateDBase
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UpdateDBase wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UpdateDBase_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in rxn2mappushbutton1.
function rxn2mappushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to rxn2mappushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end

Text1 = 'Updating rxn2map.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/link/path/rn');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','rxn2map.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('rxn2map.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end
% --- Executes on button press in map2cpdpushbutton.
function map2cpdpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to map2cpdpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating map2cpd.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/link/cpd/path');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','map2cpd.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('map2cpd.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end

% --- Executes on button press in KEGGmapspushbutton.
function KEGGmapspushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to KEGGmapspushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating KEGGmaps.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/list/pathway');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','KEGGmaps.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('KEGGmaps.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end

% --- Executes on button press in KEGGcpdpushbutton.
function KEGGcpdpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to KEGGcpdpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating KEGGcpd.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/list/cpd');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','KEGGcpd.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('KEGGcpd.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end

% --- Executes on button press in ec2rxnpushbutton.
function ec2rxnpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ec2rxnpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating ec2rxn.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/link/rn/ec');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','ec2rxn.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('ec2rxn.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end

% --- Executes on button press in cpd2weightpushbutton.
function cpd2weightpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cpd2weightpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating cpd2weight.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/find/compound/1-9000/mol_weight');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','cpd2weight.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('cpd2weight.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end

% --- Executes on button press in cpd2rxnpushbutton.
function cpd2rxnpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cpd2rxnpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating cpd2rxn.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/link/rn/cpd');
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','cpd2rxn.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('cpd2rxn.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end


% --- Executes on button press in Donepushbutton8.
function Donepushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to Donepushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


% --- Executes on button press in Weight4TransRxnspushbutton.
function Weight4TransRxnspushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Weight4TransRxnspushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating W4TranRxns2.txt and W4TranRxns4.txt...';
fprintf(Text1)
fprintf('\n')
Weight4TransRxns(2)
fprintf('\n')
Weight4TransRxns(4)
fprintf('\n')
waitfor(msgbox('W4TranRxns2.txt and W4TranRxns4.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end
clc

% --- Executes on button press in getKGMLpushbutton.
function getKGMLpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to getKGMLpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating KEGGmaps folder...';
fprintf([Text1,'\n'])
getKGML
waitfor(msgbox('KEGGmaps folder was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end

% --- Executes on button press in MultiRxnpushbutton.
function MultiRxnpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to MultiRxnpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating OfflineRxnData folder...';
fprintf(Text1)
GetData4MultiRxns
waitfor(msgbox('OfflineRxnData folder was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end
clc

% --- Executes on button press in Orgspushbutton.
function Orgspushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Orgspushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TestLink='http://rest.kegg.jp';
[~,Stat1]=urlread(TestLink);
if ~Stat1
    msgbox({'This step requires a stable internet connection';...
        'Seemingly you are offline!'});
    return
end
Text1 = 'Updating orgs.txt...';
fprintf(Text1)
Data=urlread('http://rest.kegg.jp/list/organism');
Pth1 = which ('Precons.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth1 = [Pth,'\','orgs.txt'];
Idf = fopen(Pth1,'w');
fprintf(Idf,Data);
fclose(Idf);
waitfor(msgbox('orgs.txt was successfully updated!','Update'))
for i =1:numel(Text1)
    fprintf('\b')
end
