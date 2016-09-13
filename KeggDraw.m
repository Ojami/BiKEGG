function KeggDraw (Outflx,RxnCds,MapChoice,flxType,InOpts)
% KeggDraw
% uses KGML files to extract KEGG rxn IDs and their
% corresponding coordinates on each map image (extracted from
% rest.kegg.jp). Based on reactions in each map, a comparison between
% current map reactions and input reaction codes is performed to choose
% those that are present in the selected pathway. Then, using flx2col and
% performing the normalization, flux values will be placed on corresponding
% rxn boxes in form of colour intensity. Model names provided in input flux
% data are used to check if data belong to one|two models. In latter, the
% box is divided in two separate and equal sections, each for one model.
% In case of dynamic data, the data are displayed on each pathway,
% consecutively. This would allow for a better understanding of the dynamic
% behaviour of target systems.
% 
% Inputs:
% Outflx = A mXn matrix containing input flux data, which m corresponds to
%        KEGG reaction IDs and n corresponds to time-series values.
% RxnCds = A cell array of strings comprising of KEGG reaction IDs (m rows).
% MapChoice = KEGG map ID.
% flxType = An integer: 1- Entirely fills the boxes and 2- divides them into
%           two identical parts.
% 
% Optional inputs:
% InOpts = A struct containing the following fields:
% SVal = Check if user wants to save images : 1 (yes) or 0 (no), default:0.
% mTime = Time points corresponding to time-series values.
% Thresh = Threshold below which the flux data won't be depicted (default:
%          1e-3).
% Timevalue = The time step between two consecutive images (default: 0.1 s).
% TimeUnit = A string representing time unit for time-series values
%           (default: 'No unit').
% Netstat = Cehck if user wants to use KEGG maps and KGML files online(1)
%           /offline(0) (default = 1). Data in KEGGmpas folder will be
%           used in case of offline mode.
% Svid = Check if user wants to convert a sequence of images to an animated
%        file: 1 (video), 2 (GIF), 0 (no); default: 0.
% Viso = Check if user wants to adjust visual properties (text color, box
%        color and font characteristics) through a GUI: 1 (activate GUI), 0
%        (deactivate it); default: 0.
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Jan. 2016

% Optional inputs ---------------------------------------------------------
if exist('InOpts','var')
    OptFields = fieldnames(InOpts);
    for count = 1:numel(OptFields)
        if strcmpi('SVal',OptFields{count})
            SVal = InOpts.(OptFields{count});
        elseif strcmpi('mTime',OptFields{count})
            mTime = InOpts.(OptFields{count});
        elseif strcmpi('Thresh',OptFields{count})
            Thresh = InOpts.(OptFields{count});
        elseif strcmpi('Timevalue',OptFields{count})
            Timevalue = InOpts.(OptFields{count});
        elseif strcmpi('TimeUnit',OptFields{count})
            TimeUnit = InOpts.(OptFields{count});
        elseif strcmpi('Netstat',OptFields{count})
            Netstat = InOpts.(OptFields{count});
        elseif strcmpi('Svid',OptFields{count})
            Svid = InOpts.(OptFields{count});
        elseif strcmpi('Viso',OptFields{count})
            Viso = InOpts.(OptFields{count});
        end
    end
end
% -------------------------------------------------------------------------
% Check internet connection
if ~exist('Netstat','var') || isempty(Netstat)
    Netstat = 1;
end
TestLink='http://rest.kegg.jp';
if Netstat
    [~,Stat1]=urlread(TestLink);
    if ~Stat1
        msgbox({'This step requires a stable internet connection';...
            'Seemingly you are offline!'});
        return
    end
end
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
%==========================================================================
BiGG4KeggDraw = getappdata(0,'BiGG4KeggDraw'); % For further modification
for mpct = 1:numel(MapChoice)
    h=waitbar(0,'Initializing KeggDraw...');
    if Netstat
        Link=['http://rest.kegg.jp/get/rn',MapChoice{mpct},'/kgml'];
        ImgLink=['http://www.kegg.jp/kegg/pathway/map/map',MapChoice{mpct},'.png'];
        I=imread(ImgLink);
        Dat=urlread(Link);
    else
        load([Pth,'\KEGGmaps\',MapChoice{mpct},'.mat']);
    end
    waitbar(20/100); % 20% of the process
    % Looking for KEGG rxn IDs and corresponding coordinates===============
    RawRNs = regexp(Dat,'(?<=reaction="rn:)[^"]*','match');
    IdRec = regexp(Dat,'\<rectangle" x="\d*" y="\d*" \w*="\d*" \w*="\d*"');
    IdLine = regexp(Dat, '(?<=type="line" coords=")[^"]*');
    CoordID = [IdRec,IdLine]; CoordID = sort(CoordID);
    Recseg = regexp(Dat,'\<rectangle" x="\d*" y="\d*" \w*="\d*" \w*="\d*"',...
        'match');
    RecsegCel = regexp(Recseg', '\d+', 'match');
    Reccoord = cellfun(@str2double, RecsegCel, 'UniformOutput', false);
    Lineseg = regexp(Dat, '(?<=type="line" coords=")[^"]*', 'match');
    LinesegCel = regexp(Lineseg', '\d+', 'match');
    Linecoord1 = cellfun(@str2double, LinesegCel, 'UniformOutput', false);
    Linecoord = cell(numel(Linecoord1),1);
    for linect = 1:numel(Linecoord1)
        TempLinecoord = Linecoord1{linect};
        TempLineX = TempLinecoord(1:2:end);
        TempLineY = TempLinecoord(2:2:end);
        Linecoord{linect} = [sum(TempLineX)./numel(TempLineX)+22,...
            sum(TempLineY)./numel(TempLineY)+8];
    end
    AllCoord = cell(1,numel(RawRNs));
    AllCoord(ismember(CoordID,IdLine)) = Linecoord;
    AllCoord(ismember(CoordID,IdRec)) = Reccoord;
    if ~isempty(Reccoord)
        TempW4Line = Reccoord{1}(3);
        TempH4Line = Reccoord{1}(4);
    else
        TempW4Line = 46;
        TempH4Line = 17;
    end
    idx = zeros(numel(AllCoord),1);idy = zeros(numel(AllCoord),1);
    idw = zeros(numel(AllCoord),1);idh = zeros(numel(AllCoord),1);
    for coordct = 1:numel(AllCoord)
        if numel(AllCoord{coordct}) < 4
            AllCoord{coordct}(3) = TempW4Line;
            AllCoord{coordct}(4) = TempH4Line;
        end
        idx(coordct) = round(AllCoord{coordct}(1));
        idy(coordct) = round(AllCoord{coordct}(2));
        idw(coordct) = AllCoord{coordct}(3);
        idh(coordct) = AllCoord{coordct}(4);    
    end
    LenRawRN = numel(RawRNs)+1;
    RN = RawRNs;
    if numel(RN)~=numel(idx)
        waitfor(msgbox(['Number of rxns and rectangles are not same ',...
            'for map:',MapChoice{mpct}]));
        close(h)
        continue
    end

    for rnct = 1:numel(RawRNs)
        if numel(RawRNs{rnct}) > 6
            TempCurRN = regexp(RawRNs{rnct}, 'R\d{5}', 'match');
            for rnct1 = 2:numel(TempCurRN)
                RN{LenRawRN} = TempCurRN{rnct1};
                idx(LenRawRN) = idx(rnct);
                idy(LenRawRN) = idy(rnct);
                idw(LenRawRN) = idw(rnct);
                idh(LenRawRN) = idh(rnct);
                LenRawRN = LenRawRN + 1;
            end
            RN{rnct} = RawRNs{rnct}(1:6);
        end
    end
    waitbar (40/100); % 40% of the process
    % By experience, the difference of .png images and coordinates in KGML
    idy=idy-8;
    idx=idx-22;
    %==========================================================================
    %% Analyzing flux data ====================================================
    % RxnCds should be compared to RN, in order to identify which reactions
    % from the Outflx data are present in current KEGG map
    if isempty(Outflx)
        msgbox('Flux data must be read first!','Error','error');
        return
    end

    % Cehck rxnCodes with HMR IDs and replace them with KEGG ones==========
    Hind = strfind(RxnCds,'HMR');
    if sum([Hind{:}])
        fname = 'HMR2KEGG.mat';
        if ~exist(fname,'file')
             waitfor(msgbox(['The ',fname,' Does not',...
                 ' exist in local directory.',...
                 'Open it manually'],'Warning','warn'));
             [pname,fname1]=uigetfile({'*.mat';'*.xlsx'},...
                 ['Select the ',fname,' file']);
             D=importdata([fname1,pname]);
        else
            D= importdata(fname);
        end
        KEGG=D.KEGG;HMR=D.HMR;
        [Fr,Sn]=ismember(RxnCds,HMR);
        for ind=1:length(Fr)
            if Fr(ind)
                RxnCds{ind}=KEGG{Sn(ind)};
            end
        end
    end 
    %======================================================================
    [rxnCds,flx1] = ConsistRxns(BiGG4KeggDraw,RxnCds,Outflx,MapChoice{mpct});
    waitbar (60/100); % 50% of the process
    [f1,~]=ismember(rxnCds,RN);
    flx=flx1(f1,:);
    % Check repeated rxns on this map and let user choose one,among several
    % BiGG IDs, for which the flux value will be displayed on the map.
    [flx,rxnCds] = ChckRptRxns(flx,f1,rxnCds,MapChoice{mpct});
    waitbar(80/100); % 80% of the process
    copt2 = 1; idxt=(0);idyt=(0);idwt=(0);idht=(0);flx2=zeros(1,size(flx,2));
    rxnCds1={0};
    for copt = 1:numel(rxnCds)
        Temprxn = find(ismember(RN,rxnCds{copt}));
        for copt1 = 1:numel(Temprxn)
            idxt(copt2) = idx(Temprxn(copt1));
            idyt(copt2) = idy(Temprxn(copt1));
            idwt(copt2) = idw(Temprxn(copt1));
            idht(copt2) = idh(Temprxn(copt1));
            flx2(copt2,:) = flx(copt,:);
            rxnCds1{copt2} = rxnCds{copt};
            copt2 = copt2 + 1;
        end
    end
    flx = flx2;
    % Check multi-step reactions and identify if all the child rxns are
    % present in the current pathway, otherwise they won't be considered.
    % Also, there are reactions for which there is more than one BiGG ID.
    % Note that this case is different from ChckRptRxns function, because
    % these reactions are parts of multi-step reactions, and therefore, are
    % not equivalent.
    [OutFlx,IdxN,IdyN,IdwN,IdhN] = ChckMultiRxns(RN,idx,idy,idw,idh);
    waitbar (90/100); % 90% of the process
    
    flx2 = [flx;OutFlx];
    idxt1=[idxt';IdxN]; idyt1=[idyt';IdyN]; idwt1=[idwt';IdwN];idht1=[idht';IdhN];
    flx =flx2; idxt=idxt1;idyt=idyt1;idwt=idwt1;idht=idht1;
    %----------------------------------------------------------------------
    if isempty(flx)
        uiwait(msgbox('The flux data are not is current pathway!',...
            'Warning','warn'));
        close (h);
        return
    end
    %======================================================================
    %% Visualization=======================================================
    % Threshold below which the flux data won't be depicted:
    if ~exist('Thresh','var') || isempty(Thresh)
        Thresh=1e-6;
    end
    Cm=jet(size(idx,1)).*255; % [Warning] Changing Cm matrix is not recommended

    % The time step between two consecutive images:
    if ~exist('Timevalue','var') || isempty(Timevalue)
        Timevalue=0.1;
    end
    waitbar(100/100); % 100% of the process
    close (h);

    % Time points:
    if ~exist('mTime','var') || isempty(mTime)
        if flxType == 2
            mTime=1:size(Outflx,2)/2;
        else
            mTime=1:size(Outflx,2);
        end
    end

    if ~exist('TimeUnit','var') || isempty(TimeUnit)
        TimeUnit='No unit';
    end
    if ~exist('Svid','var') || isempty(Svid)
        Svid = 0; %Does not convert images to video
    end
    % Check wether user wants to save images (using 'Save images checkbox')
    if ~exist('SVal','var') || isempty(SVal)
        SVal = 0; % Does not save images.
    end
    if SVal
        FolderNm=[MapChoice{mpct},' images files'];
        Tn=uigetdir(pwd,['Folder to save images. To save in current folder:',...
        'Cancel']);
        if ~Tn
            if ~exist(FolderNm,'file') % Check if folder already exists
                mkdir(FolderNm);
                FolderPath=fullfile(pwd,FolderNm);
            else
                FolderPath=fullfile(pwd,FolderNm);
            end
        else
            if ~exist(fullfile(Tn,FolderNm),'file') % Check if folder already exists
                mkdir(fullfile(Tn,FolderNm));
                FolderPath=fullfile(Tn,FolderNm);
            else
                FolderPath=fullfile(Tn,FolderNm);
            end
        end
    end
    %======================================================================
    %Draw images on this panel
    figure(mpct)
    H1 = gcf;
    h1 = gca;
    set(H1,'position',get(0,'screensize'));
    set(h1,'position',[0,0,1,1],'xtick',[],'ytick',[]);
    set(H1,'toolbar','figure');
    set(H1,'menubar','figure');
    set(H1,'name',MapChoice{mpct},'numbertitle','off')
    % Arbitrary color and font properties
    if ~exist('Viso','var') || ~Viso
        TextFont.FontName = 'Arial';
        TextFont.FontWeight = 'normal';
        TextFont.FontAngle = 'normal';
        TextFont.FontUnits = 'points';
        TextFont.FontSize = 8;
        TextCol = [0,0,0];
        BoxCol = [240,248,255];
    end
    if exist('Viso','var') && Viso
        waitfor(VisualProp)
        BoxCol = getappdata(0,'BoxtCol');
        if ~isempty(BoxCol)
            BoxCol = BoxCol.*255;
        else
            BoxCol = [240,248,255];
        end
        TextFont = getappdata(0,'TextFont');
        if isempty(TextFont)
            TextFont.FontName = 'Arial';
            TextFont.FontWeight = 'normal';
            TextFont.FontAngle = 'normal';
            TextFont.FontUnits = 'points';
            TextFont.FontSize = 8;
        end
        TextCol = getappdata(0,'TextCol');
        if isempty(TextCol)
            TextCol = [0,0,0];
        end
    end
    switch flxType
        case 1
            for k=1:size(flx,2)
                for m=1:size(flx,1)
                    % Normalize colors and insertion into enzyme boxes
                    %flxMax=max(abs(flx(:,k)));
                    %flxMin=min(abs(flx(:,k)));
                    flxMax=max(abs(flx(:)));
                    flxMin=min(abs(flx(:)));
                    if flxMax==flxMin % One row (i.e. for 1 rxn)
                        %flxMax=max(abs(flx(m,:)));
                        %flxMin=min(abs(flx(m,:)));
                        flxMax=max(abs(flx(:)));
                        flxMin=min(abs(flx(:)));
                    end
                    Tempx=(idyt(m)+1:idyt(m)+idht(m)-1);
                    I= flx2col(BoxCol,flx,Cm,I,k,flxMax,flxMin,m,idxt,idwt,Tempx,Thresh);
                    % Show the image:
                    imshow(I,'InitialMagnification', 'fit','Parent',h1); 
                    % Insert flux value on each reaction box
                    TextPropUDF (TextFont,TextCol,'One',idxt,idyt,idwt,...
                        idht,k,flx,h1,Thresh)

                    % Show time-series values on the corresponding image
                    TimePropUDF ('One',TextFont,TimeUnit,mTime,k,h1) 
                end
                % Show colorbar
                flxstr=linspace(flxMin,flxMax,7);
                flxstr1=cell(numel(flxstr),1);
                for ir=1:numel(flxstr)
                    flxstr1{ir,1}=num2str(flxstr(ir),2);
                end
                ch = colorbar('peer',h1);
                set(ch,'yticklabel',flxstr1)

                % Save to image files within a new folder======================
                if SVal
                    if k == 1
                        waitfor(SaveImgs)
                        ImgRes = getappdata(0,'ImgRes');
                        ImgFrmt = getappdata(0,'ImgFrmt');
                        if isempty(ImgRes)
                            ImgRes=1; % Corresponds to 300 dpi in SaveFlxImgs
                        elseif isempty(ImgFrmt)
                            ImgFrmt=1; % Corresponds to PNG in SaveFlxImgs
                        end
                    end
                    SaveFlxImgs(H1,FolderPath,k,ImgRes,ImgFrmt);
                end
                %==============================================================
                pause(Timevalue); % Time step between consecutive images
            end
        case 2
            for k=1:2:size(flx,2)
                for m=1:size(flx,1)
                    % Normalize colors and insertion into enzyme boxes
                    %flxMax=max(max(abs(flx(:,k))),max(abs(flx(:,k+1))));
                    %flxMin=min(min(abs(flx(:,k))),min(abs(flx(:,k+1))));
                    flxMax=max(abs(flx(:)));
                    flxMin=min(abs(flx(:)));
                    if size(flx,1)==1
                        %flxMax=max(abs(flx(m,:)));
                        %flxMin=min(abs(flx(m,:)));
                        flxMax=max(abs(flx(:)));
                        flxMin=min(abs(flx(:)));
                    end
                    Boxheight_mid=(idht(m)-2)/2;
                    Tempx=(idyt(m)+1:idyt(m)+idht(m)-1-Boxheight_mid);
                    Tempx1=round((idyt(m)+idht(m)-1-Boxheight_mid:idyt(m)...
                        +idht(m)-1));
                    % For 1st flx set:
                    I= flx2col(BoxCol,flx,Cm,I,k,flxMax,flxMin,m,idxt,idwt,Tempx,Thresh); 
                    % For 2nd flx set:
                    I= flx2col(BoxCol,flx,Cm,I,k+1,flxMax,flxMin,m,idxt,idwt,Tempx1,Thresh);
                    % Show the image:
                    imshow(I,'InitialMagnification', 'fit','Parent',h1); 

                    % Insert flux value on each reaction box for 1st set(upper)
                    TextPropUDF (TextFont,TextCol,'Two',idxt,idyt,...
                        idwt,idht,k,flx,h1,Thresh)
                    % Show time-series values on the corresponding image
                    TimePropUDF ('Two',TextFont,TimeUnit,mTime,k,h1) 
                end

                % Show colorbar
                flxstr=linspace(flxMin,flxMax,7);
                flxstr1=cell(numel(flxstr),1);
                for ir=1:numel(flxstr)
                    flxstr1{ir,1}=num2str(flxstr(ir),2);
                end
                ch = colorbar('peer',h1);
                set(ch,'yticklabel',flxstr1)

                % Save to image files within a new folder======================
                if SVal
                    if k == 1
                        waitfor(SaveImgs)
                        ImgRes = getappdata(0,'ImgRes');
                        ImgFrmt = getappdata(0,'ImgFrmt');
                        if isempty(ImgRes)
                            ImgRes=1; % Corresponds to 300 dpi in SaveFlxImgs
                        elseif isempty(ImgFrmt)
                            ImgFrmt=1; % Corresponds to PNG in SaveFlxImgs
                        end
                    end
                    SaveFlxImgs(H1,FolderPath,ceil(k/2),ImgRes,ImgFrmt);
                end
                %==============================================================
                pause(Timevalue); % Time step between consecutive images
            end
    end
if Svid
    Img2Animation(Svid)
end
end

% Subfunctions ------------------------------------------------------------
function [flxOut,rxnCdsI] = ChckRptRxns (flx,f1,rxnCds,CurrentMap)
% First, rxns are pruned based on manual assessment of pathways
BiGG4KeggDraw = getappdata(0,'BiGG4KeggDraw1');
rxnCdsI = rxnCds(f1);
BiGGI = BiGG4KeggDraw(f1);
Fileid1 = fopen('rxnbiggmap.txt','r');
Rxnbmap = textscan(Fileid1,'%s %s %s %d');
Rtemp = Rxnbmap{1}; Btemp = Rxnbmap{2}; Mtemp = Rxnbmap{3};
Dtemp = Rxnbmap{4};
if any(ismember(Mtemp,CurrentMap))
    ChR = Rtemp(ismember(Mtemp,CurrentMap));
    ChB = Btemp(ismember(Mtemp,CurrentMap));
    ChD = Dtemp(ismember(Mtemp,CurrentMap));
    for ct = 1:numel(ChR)
        n1_chr = find(ismember(rxnCdsI,ChR{ct}));
        n1_chb = find(ismember(BiGGI,ChB{ct}));
        if ~isempty(intersect(n1_chr,n1_chb))
            %Rmvit = find(ismember(BiGGI,ChB{ct}));
			Rmvit = intersect(n1_chr,n1_chb);
            for ct1 = 1:numel(Rmvit)
                if strcmp(BiGGI{Rmvit(ct1)},ChB{ct})
                    Tmpflx = flx(Rmvit(ct1),:);
                end
            end
            rxnCdsI(Rmvit) = [];
            BiGGI(Rmvit) = [];
            flx(Rmvit,:) = [];
            if ChD(ct)
                CurL1 = numel(rxnCdsI) + 1;
                rxnCdsI{CurL1} = ChR{ct};
                BiGGI{CurL1} = ChB{ct};
                flx(CurL1,:) = Tmpflx;
            end
        end
    end
end
%------
[R1,~,R3] = unique(rxnCdsI,'stable');
countR = 1;RptCount=(0);
while numel(R3)
    RptCount(countR) = numel(find(R3==R3(1)));
    R3(R3==R3(1)) = [];
    countR = countR+1;
end
FndRxnLoci = find(RptCount>1);
RptRxns = {0}; RptBiggs = {0};RptFlx= {0};
for tctr = 1:numel(FndRxnLoci)
    RptRxnsI = R1(FndRxnLoci(tctr));
    RptRxns{tctr} = rxnCdsI(ismember(rxnCdsI,RptRxnsI));
    RptBiggs{tctr} = BiGGI(ismember(rxnCdsI,RptRxnsI));
    RptFlxLoci = ismember(rxnCdsI,RptRxnsI);
    RptFlx{tctr} = flx(RptFlxLoci,:);
end
if any(RptCount>1) 
    setappdata(0,'RptBiggs',RptBiggs);
    setappdata(0,'CurrentMap',CurrentMap);
    setappdata(0,'RptRxns',RptRxns);
    setappdata(0,'RptFlx',RptFlx);
    waitfor(KeggDrawTable)
    SlctdIds = getappdata(0,'data1');
    if ~isempty(SlctdIds)
        NtSel = RptBiggs(ismember(SlctdIds,'Not slected'));
        OwnSel = StubSel(NtSel);
        SlctdIds(ismember(SlctdIds,'Not slected')) = OwnSel;
    else
        SlctdIds = StubSel(RptBiggs);
    end
    % Refine flx
    for count = 1:numel(SlctdIds)
        Fi = find(ismember(RptBiggs{count},SlctdIds{count}));
        for count1 = 1:size(RptFlx{count},1)
            RptFlx{count}(count1,:) = RptFlx{count}(Fi,:);
        end
    end
    for count1 = 1:size(RptRxns,2)
        Fi1 = find(ismember(rxnCdsI,RptRxns{count1}));
        flx(Fi1,:) = RptFlx{count1};
    end
    flxOut = flx;
    
else
    flxOut = flx;
end
% Remove all appdata
AT=getappdata(0);
AT1 = fieldnames(AT);
for co = 1:length(AT1)
    if ~strcmp(AT1{co},'BiGG4KeggDraw') && ~strcmp(AT1{co},'MultiB')...
            && ~strcmp(AT1{co},'MultiK') && ~strcmp(AT1{co},'MultiFlx')
        rmappdata(0,AT1{co});
    end
end
%--------------------------------------------------------------------------
function OutSel = StubSel(InSel)
OutSel = {0};
for count1 = 1:numel(InSel)
    Tinder = regexp(InSel{count1},'\w*p$'); %Peroxisome
    Tinder = find(~cellfun('isempty', Tinder));
    CurSel = InSel{count1};
    if (numel(CurSel)-numel(Tinder))==1
        OutSel{count1} = CurSel{~ismember(CurSel,CurSel(Tinder))};
    else %Nothing found or code cannot decide
        OutSel{count1} = CurSel{1};
    end
end
%--------------------------------------------------------------------------
function [OutFlx,IdxN,IdyN,IdwN,IdhN] = ChckMultiRxns(RN,idx,idy,idw,idh)
MultiB = getappdata(0,'MultiB'); 
MultiK = getappdata(0,'MultiK');
MultiFlx = getappdata(0,'MultiFlx');

OutFlx =[];
IdxN=[];IdyN=[];IdwN=[];IdhN=[];
while numel(MultiB)
    Bigg1Loci = ismember(MultiB,MultiB{1});
    KeggTemp = MultiK(Bigg1Loci);
    if all(ismember(KeggTemp,RN)) % All child rxns are present in map
        FlxTemp = MultiFlx(Bigg1Loci,:);
        [Tchck1,Tchck2] = ismember(RN,KeggTemp);
        if sum(Tchck1) > numel(KeggTemp) % Some KggTemp elements have more than
            % one equivalent in RN, meaning more than one set of coordinates
            % exists
            Tchck2(Tchck2==0)=[];
            U = unique(Tchck2);
            Nch = histc(Tchck2,U);
            RptFound = U(Nch>1);
            RptRxns = KeggTemp(RptFound); % RptRxns contains repeated rxns
            % in RN, which are shared between RN and KeggTemp.
            OtherRxns1 = KeggTemp(~ismember(KeggTemp,RptRxns));
            OtherIdx = idx(ismember(RN,OtherRxns1));
            OtherIdy = idy(ismember(RN,OtherRxns1));
            OtherIdw = idw(ismember(RN,OtherRxns1));
            OtherIdh = idh(ismember(RN,OtherRxns1));
            % Assumption: all child rxns are placed close together, so,
            % their coordinates are similar. Therefore, RptRxns having
            % similar coordinates to mean(OtherIdx & OtherIdy) wille be
            % selected.
            for count1 = 1:numel(RptRxns)
                MeanIdx = mean(OtherIdx); MeanIdy = mean(OtherIdy);
                RptIdx = idx(ismember(RN,RptRxns{count1}));
                RptIdy = idy(ismember(RN,RptRxns{count1}));
                RptIdw = idw(ismember(RN,RptRxns{count1}));
                RptIdh = idh(ismember(RN,RptRxns{count1}));
                RptDist = sqrt((RptIdx-MeanIdx).^2+(RptIdy-MeanIdy).^2);
                [~,Id1] = min(RptDist);
                CurLen = numel(OtherIdx)+1;
                OtherIdx(CurLen) = RptIdx(Id1);
                OtherIdy(CurLen) = RptIdy(Id1);
                OtherIdw(CurLen) = RptIdw(Id1);
                OtherIdh(CurLen) = RptIdh(Id1);
            end
        else
            OtherIdx = idx(ismember(RN,KeggTemp));
            OtherIdy = idy(ismember(RN,KeggTemp));
            OtherIdw = idw(ismember(RN,KeggTemp));
            OtherIdh = idh(ismember(RN,KeggTemp));
        end
        if size(OutFlx,1) == 0
            OutFlx = FlxTemp;
            IdxN = OtherIdx; IdyN = OtherIdy; IdwN = OtherIdw; 
            IdhN = OtherIdh;
        else
            OutFlx(size(OutFlx,1)+1:size(OutFlx,1)+size(FlxTemp,1),:)...
                = FlxTemp;
            IdxN(numel(IdxN)+1:numel(IdxN)+numel(OtherIdx)) = OtherIdx;
            IdyN(numel(IdyN)+1:numel(IdyN)+numel(OtherIdy)) = OtherIdy;
            IdwN(numel(IdwN)+1:numel(IdwN)+numel(OtherIdw)) = OtherIdw;
            IdhN(numel(IdhN)+1:numel(IdhN)+numel(OtherIdh)) = OtherIdh;
        end
    end
    MultiB(Bigg1Loci) = [];
    MultiK(Bigg1Loci) = [];
    MultiFlx(Bigg1Loci,:) = [];
end
% -------------------------------------------------------------------------
function [rxnOut,flxOut] = ConsistRxns(BiGG4KeggDraw,rxnIn,flxIn,MapChoice)

% Read rxn2map ========================
Fileid1 = fopen('rxn2map.txt','r');
rxn2map = textscan(Fileid1,'%s %s');
RawRxns1 = rxn2map{1};
RawMaps = rxn2map{2};
[Loci1,~] = regexp(RawMaps, 'path:rn');
Loci2 = ~cellfun('isempty', Loci1);
RawMaps(Loci2) = [];
RawRxns1(Loci2) = [];
Maps = strrep(RawMaps,'path:map','');
Rxns1 = strrep(RawRxns1,'rn:',''); 
% Read ec2rxn =========================
Fileid1 = fopen('ec2rxn.txt','r');
ec2rxn = textscan(Fileid1,'%s %s');
RawEc= ec2rxn{1};
RawRxns = ec2rxn{2};
EC = strrep(RawEc,'ec:','');
Rxns = strrep(RawRxns,'rn:','');
FoundRxnsID = (0); RxnChk = rxnIn;
% Read cpd2rxn =======================
Fileid = fopen('cpd2rxn.txt','r');
cpd2rxn = textscan(Fileid,'%s %s');
rxnData = strrep(cpd2rxn{2},'rn:','');
cpds = strrep(cpd2rxn{1},'cpd:','');
% Read cpd2weight ===================
Fileid = fopen('cpd2weight.txt','r');
cpd2weight = textscan(Fileid,'%s %f');
cpd4weight = strrep(cpd2weight{1},'cpd:','');
Weight = cpd2weight{2};
% Read map2cpd ===================
% Fileid = fopen('map2cpd.txt','r');
% map2cpd = textscan(Fileid,'%s %s');
% map4cpdT = map2cpd{1};
% cpd4mapT = map2cpd{2};
% Tid = regexp(map4cpdT,'path:map\d*');
% Tid = ~cellfun('isempty', Tid);
% map4cpd = map4cpdT(Tid);
% map4cpd = strrep(map4cpd,'path:map','');
% cpd4map = cpd4mapT(Tid);
% cpd4map = strrep(cpd4map,'cpd:','');
fprintf('\nRxn consistency for map%s(of %d): ',MapChoice,numel(rxnIn))
ctr = 1;
for count = 1:numel(rxnIn)
    Rns4Map = Maps(ismember(Rxns1,rxnIn{count}));
      if count>1
          for j=0:log10(count-1)
              fprintf('\b'); 
          end
      end
      fprintf('%d', count);
    if ~any(ismember(Rns4Map,MapChoice)) %Otherwise there is no need to identify
        % reactions for this map
        FoundRxnsID(ctr) = count;
        ctr = ctr+1;
        TempEC=EC(ismember(Rxns,rxnIn{count}));
        TempRns = Rxns(ismember(EC,TempEC));
        TempRns = unique(TempRns);
        TempRns(ismember(TempRns,rxnIn{count})) = [];
        TempFlx = flxIn(count,:);
        CurLen = numel(rxnIn) + 1;
        cpds4curRxn = cpds(ismember(rxnData,rxnIn{count}));
        W4curCpds = Weight(ismember(cpd4weight,cpds4curRxn));
        for count1 = 1:numel(TempRns)
            TempRns2Map1 = Maps(ismember(Rxns1,TempRns{count1}));
            ChkCur = any(ismember(RxnChk,TempRns{count1}));
            if any(ismember(TempRns2Map1,MapChoice)) && ~ChkCur
                cpds4TempRxn = cpds(ismember(rxnData,TempRns{count1}));
%                 DiffCpds4Temp = cpds4TempRxn(~ismember(cpds4TempRxn,cpds4curRxn));
%                 Maps4TempRxn = map4cpd(ismember(cpd4map,DiffCpds4Temp));
%                 ChckTempMap1 = sum(ismember(Maps4TempRxn,MapChoice));
%                 ChckTempMap = (numel(DiffCpds4Temp) == ChckTempMap1);
                W4TempCpds = Weight(ismember(cpd4weight,cpds4TempRxn));
                if numel(W4curCpds) == numel(W4TempCpds) &&...
                        isempty(setdiff(W4curCpds,W4TempCpds)) %&& ChckTempMap
                    rxnIn{CurLen} = TempRns{count1};
                    flxIn(CurLen,:) = TempFlx;
                    BiGG4KeggDraw{CurLen} = BiGG4KeggDraw{count};
                    CurLen = CurLen+1;
                end
            end
        end
    end
end
fprintf('\n')
if FoundRxnsID
    flxIn(FoundRxnsID,:) = [];
    rxnIn(FoundRxnsID) = [];
    BiGG4KeggDraw(FoundRxnsID) = [];
    setappdata(0,'BiGG4KeggDraw1',BiGG4KeggDraw);
end
rxnOut = rxnIn;
flxOut = flxIn;