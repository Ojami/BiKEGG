function varargout = MapAdjuster(varargin)
% MapAdjuster(GUI)
% Post-processing the created customized metabolic maps by NetDraw
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% June 2016

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MapAdjuster_OpeningFcn, ...
                   'gui_OutputFcn',  @MapAdjuster_OutputFcn, ...
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


% --- Executes just before MapAdjuster is made visible.
function MapAdjuster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MapAdjuster (see VARARGIN)

% Choose default command line output for MapAdjuster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MapAdjuster wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MapAdjuster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in DonePush.
function DonePush_Callback(hObject, eventdata, handles)
PostData = getappdata(0,'PostData');
% Title processing ////////////////////////////////////////////////////////
Title = get(handles.TitleEdit,'String');
FHndl = getappdata(0,'ParentFig');
AxHndl = getappdata(0,'ParentAx');
% TitleHndl = getappdata(hObject,'TitleHndl');
% set(TitleHndl,'Visible','off');
text('position',[50,50],'String',Title,'Unit','data','FontSize',22,...
    'Parent',AxHndl);
% setappdata(hObject,'TitleHndl',TitleHndl);
set(FHndl,'name',Title,'numbertitle','off')
% Show Rates///////////////////////////////////////////////////////////////
ShowRate = get(handles.ShowRatesCheck,'Value');
if ShowRate
    FlxStrHndls = getappdata(hObject,'FlxStrHndls');
    set(FlxStrHndls,'Visible','off');
    AxHndl = getappdata(0,'ParentAx');
    flx = PostData.flx(:,PostData.k);
    ArX = PostData.ArX;
    ArY = PostData.ArY;
    for i1=1:numel(ArX)
        FlxStrHndls(i1) = text('position',[mean(ArX{i1}),mean(ArY{i1})],...
         'String',num2str(abs(flx(i1)),'%.2f'),'FontSize',7,'Parent',AxHndl);
    end
    setappdata(hObject,'FlxStrHndls',FlxStrHndls);
else
    FlxStrHndls = getappdata(hObject,'FlxStrHndls');
    set(FlxStrHndls,'Visible','off');
end
% Show compounds //////////////////////////////////////////////////////////
CpdPanelTags = {'KeggCpdIdRadio','KeggCpdNRadio','BiGGCpdRadio','RmvCpdRadio'};
hCpdpanel = get(handles.CpdPanel,'SelectedObject');
Cpdpanel = find(strcmp(CpdPanelTags,get(hCpdpanel,'Tag')));
AxHndl = getappdata(0,'ParentAx');
KEGGCpdH = getappdata(hObject,'KEGGCpdH');
KEGGCpdNH = getappdata(hObject,'KEGGCpdNH');
BiGGCpdH = getappdata(hObject,'BiGGCpdH');
if any(ishandle(KEGGCpdH))
    set(KEGGCpdH,'Visible','off');
end
if any(ishandle(KEGGCpdNH))
    set(KEGGCpdNH,'Visible','off');
end
if any(ishandle(BiGGCpdH))
    set(BiGGCpdH,'Visible','off');
end
switch Cpdpanel
    case 1
        cpdnames = PostData.cpdnames;
        x_cpd = PostData.x_cpd;
        y_cpd = PostData.y_cpd;
        cpdxy = PostData.cpdxy;
        for i1=1:numel(cpdnames)
             KEGGCpdH(i1) = text('position',[x_cpd(i1)+cpdxy(1,4)/2,y_cpd(i1)-cpdxy(1,4)/2],...
                 'String',cpdnames{i1},'FontSize',7,'Parent',AxHndl);
        end
        setappdata(hObject,'KEGGCpdH',KEGGCpdH);
    case 2
        cpdnames = PostData.cpdnames;
        x_cpd = PostData.x_cpd;
        y_cpd = PostData.y_cpd;
        cpdxy = PostData.cpdxy;
        Fileid = fopen('KEGGcpd.txt','r');
        TempCpds = textscan(Fileid,'%s %s %[^\n]');
        fclose(Fileid);
        CpdIds = TempCpds{1}; CpdIds = strrep(CpdIds,'cpd:','');
        CpdNms = TempCpds{2}; CpdNms = strrep(CpdNms,';','');
        clear TempCpds
        [~,N1] = ismember(cpdnames,CpdIds);
         KEGGCpdNms = CpdNms(N1);
        for i1=1:numel(KEGGCpdNms)
             KEGGCpdNH(i1) = text('position',[x_cpd(i1)+cpdxy(1,4)/2,y_cpd(i1)-cpdxy(1,4)/2],...
                 'String',KEGGCpdNms{i1},'FontSize',5.5,'Parent',AxHndl);
        end
        setappdata(hObject,'KEGGCpdNH',KEGGCpdNH);
    case 3
        BiGGModel = getappdata(0,'BiGGModel');
        if isempty(BiGGModel)
            [pname,fname1]=uigetfile({'*.*'},'Select the model file (SBML)');
            Testname=pname(strfind(pname,'.')+1:end);
            if ~strcmp(Testname,'xml')
                msgbox('Only xml file format!','Error','error');
                return
            else
                Idf = fopen([fname1,pname]);
                D = textscan(Idf,'%s');
                fclose(Idf);
                CbModel = readCbModel([fname1,pname]);
                [Metkegg,Metbigg] = KEGG2BiGGMets(D,CbModel);
                BiGGModel.Metkegg = Metkegg;
                BiGGModel.CbModel = CbModel;
                BiGGModel.Metbigg = Metbigg;
                setappdata(0,'BiGGModel',BiGGModel);
            end
        end
        Metbigg = BiGGModel.Metbigg;
        Metkegg = BiGGModel.Metkegg;
        cpdnames = PostData.cpdnames;
        x_cpd = PostData.x_cpd;
        y_cpd = PostData.y_cpd;
        cpdxy = PostData.cpdxy;
        for m1 = 1:numel(cpdnames)
            if any(ismember(cpdnames{m1},Metkegg))
                whrc = find(ismember(Metkegg,cpdnames{m1}));
                cpdnames(m1) = Metbigg(whrc(1));
            end
        end
        
        for i1=1:numel(cpdnames)
             BiGGCpdH(i1) = text('position',[x_cpd(i1)+cpdxy(1,4)/2,y_cpd(i1)-cpdxy(1,4)/2],...
                 'Interpreter','none','String',cpdnames{i1},'FontSize',5.5,'Parent',AxHndl);
        end
        setappdata(hObject,'BiGGCpdH',BiGGCpdH);
end
% Show rxns ///////////////////////////////////////////////////////////////
RxnPanelTags = {'KeggRxnIdRadio','BiGGRxnRadio','RmvRxnRadio'};
hrxnpanel = get(handles.RxnPanel,'SelectedObject');
rxnpanel = find(strcmp(RxnPanelTags,get(hrxnpanel,'Tag')));
AxHndl = getappdata(0,'ParentAx');
KEGGRxnHndls = getappdata(hObject,'KEGGRxnHndls');
BiGGRxnHndls = getappdata(hObject,'BiGGRxnHndls');
if any(ishandle(KEGGRxnHndls))
    set(KEGGRxnHndls,'Visible','off');
end
if any(ishandle(BiGGRxnHndls))
    set(BiGGRxnHndls,'Visible','off');
end

switch rxnpanel % Show rxn IDs ////////////////////////////////////////////
    case 1 % KEGG ID ------------------------------------------------------
        [OverlayMe,Allx4Arrow,Ally4Arrow] = RxnAdjuster(PostData,[],'kegg');
        for i1=1:numel(Allx4Arrow)
            if std(Allx4Arrow{i1})<1e-1
                KEGGRxnHndls(i1) = text('position',...
                    [mean(Allx4Arrow{i1})-10,mean(Ally4Arrow{i1})+15],...
                 'String',OverlayMe{i1},'FontSize',6.5,'rotation',90,'Parent',AxHndl);
            else
                KEGGRxnHndls(i1) = text('position',...
                    [mean(Allx4Arrow{i1})-15,mean(Ally4Arrow{i1})-10],...
                 'String',OverlayMe{i1},'FontSize',6.5,'Parent',AxHndl);
            end
        end 
        setappdata(hObject,'KEGGRxnHndls',KEGGRxnHndls)
    case 2 % BiGG ID ------------------------------------------------------
        Mdl = getappdata(0,'Mdl');
        if isempty(Mdl)
            [pname,~]=uigetfile({'*.*'},'Select the model file (SBML)');
            Testname=pname(strfind(pname,'.')+1:end);
            if ~strcmp(Testname,'xml')
                msgbox('Only xml file format!','Error','error');
                return
            end
            Chk = which([pname(1:strfind(pname,'.')-1),'KEGG.mat']);
            if isempty(Chk)
                msgbox({[pname(1:strfind(pname,'.')-1),'KEGG.mat',' cannot be found!'];...
                    'Make sure the file is present in BiGG2KEGG folder'},'Error','error');
                return
            end
            Mdl = load(which([pname(1:strfind(pname,'.')-1),'KEGG.mat']));
            setappdata(0,'Mdl',Mdl);
        end
        [OverlayMe,Allx4Arrow,Ally4Arrow] = RxnAdjuster(PostData,Mdl,'bigg');
        % Rmv empty ones
        Chkempty = find(cellfun('isempty',regexp(OverlayMe,'\w','match')));
        OverlayMe(Chkempty) = {''};
        
        for i1=1:numel(Allx4Arrow)
            if std(Allx4Arrow{i1})<1e-1
                BiGGRxnHndls(i1) = text('position',...
                    [mean(Allx4Arrow{i1})-10,mean(Ally4Arrow{i1})+15],...
                 'String',OverlayMe{i1},'FontSize',6,'rotation',90,...
                 'Interpreter','none','Parent',AxHndl);
            else
                BiGGRxnHndls(i1) = text('position',...
                    [mean(Allx4Arrow{i1})-15,mean(Ally4Arrow{i1})-10],...
                 'String',OverlayMe{i1},'FontSize',6,'Parent',AxHndl,...
                 'Interpreter','none');
            end
        end 
        setappdata(hObject,'BiGGRxnHndls',BiGGRxnHndls)
end

function TitleEdit_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function TitleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TitleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TitlePush.
function TitlePush_Callback(hObject, eventdata, handles)



function LineWEdit_Callback(hObject, eventdata, handles)
% hObject    handle to LineWEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineWEdit as text
%        str2double(get(hObject,'String')) returns contents of LineWEdit as a double


% --- Executes during object creation, after setting all properties.
function LineWEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineWEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowRatesCheck.
function ShowRatesCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ShowRatesCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowRatesCheck


% --- Executes on button press in OverlapCheck.
function OverlapCheck_Callback(hObject, eventdata, handles)
% hObject    handle to OverlapCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OverlapCheck


% --- Executes on button press in ColormapCheck.
function ColormapCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ColormapCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ColormapCheck


% --- Executes on button press in RedrawPush.
function RedrawPush_Callback(hObject, eventdata, handles)
OverlapC = get(handles.OverlapCheck,'Value');
ColormapC = get(handles.ColormapCheck,'Value');
LineStr = str2double(get(handles.LineWEdit,'String'));
BackCol = getappdata(handles.Backgroundpush,'BackCol');
CpdCol = getappdata(handles.Cpdpush,'CpdCol');
inactiveCol = getappdata(handles.inactivespush,'inactiveCol');
NetPoster(LineStr,OverlapC,ColormapC,BackCol,CpdCol,inactiveCol)


% --- Executes on button press in Backgroundpush.
function Backgroundpush_Callback(hObject, eventdata, handles)
BackCol=uisetcolor;
setappdata(hObject,'BackCol',BackCol);


% --- Executes on button press in Cpdpush.
function Cpdpush_Callback(hObject, eventdata, handles)
CpdCol=uisetcolor;
setappdata(hObject,'CpdCol',CpdCol);


% --- Executes on button press in inactivespush.
function inactivespush_Callback(hObject, eventdata, handles)
inactiveCol=uisetcolor;
setappdata(hObject,'inactiveCol',inactiveCol);


% --- Executes on button press in ClosePush.
function ClosePush_Callback(hObject, eventdata, handles)
% hObject    handle to ClosePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
