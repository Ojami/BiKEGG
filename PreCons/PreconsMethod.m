function PreconsMethod(handles)
% PreconsMethod
% Employs Precons GUI tool for building reaction correspondences based on
% KEGG genome annotations. Output is a MS Excel spreadsheet  containing
% gene entries (KEGG), KEGG reaction identifier, BiGG reaction abbr., and 
% associated pathways (KEGG). 
% A MATLAB biograph object can also be generated for the selected organism.
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% June 2016

Org = get(handles.OrganEdit,'String');

[getpathways,stat] = urlread(['http://rest.kegg.jp/list/pathway/',Org]);
if ~stat
    msgbox('There is a problem with the internet connection!');
    return
end
pths = textscan(getpathways,'%s %[^\n]');
pathnames = pths{1}; pathnames = strrep(pathnames,'path:','');
% pathdescription = pths{2};

% Load KEGG2BiGG reactions
BiGGdata = load('UniModelKEGG.mat');
Bdata = BiGGdata.B2Kegg.B; Kdata = BiGGdata.B2Kegg.K;
Multir = find(ismember(Bdata,'MULTIR'));
Bdata(Multir:end) = []; Kdata(Multir:end) = [];
Allgenes = {[]}; Allrxns = {[]}; Allbigg = {[]}; Allpath = {[]};
yct = 1;
for cp = 1:numel(pathnames)
    fprintf('Pathway: %d of %d\n',cp,numel(pathnames))
    [getKGML,stat] = urlread(['http://rest.kegg.jp/get/',pathnames{cp},'/kgml']);
    if ~stat
        msgbox('There is a problem with the internet connection!');
        clc
        return
    end
    Ids = regexp(getKGML,'(?<=entry id=)[^>]*','match');
    if isempty(Ids)
        continue
    end
    for ct = 1:numel(Ids)
        rnsTemp = Ids{ct};
        checkgene = regexp(rnsTemp,'type="gene"', 'once');
        if ~isempty(checkgene)
            rnnms1 = regexp(rnsTemp,'(?<=reaction=")[^"]*','match');
            if ~isempty(rnnms1)
                genenms1 = regexp(rnsTemp,'(?<=name=")[^"]*','match');
                CurrxnNum1 = regexp(rnnms1, 'R\d{5}', 'match');
                CurrxnNum = CurrxnNum1{1};
                for mc = 1:numel(CurrxnNum)
                    biggTemp = Bdata(ismember(Kdata,CurrxnNum{mc}));
                    if isempty(biggTemp)
                        biggtemp = '';
                    else
                        if numel(biggTemp)>1
                             biggtemp = strjoin(biggTemp,';');
                        else
                            biggtemp = biggTemp{1};
                        end
                    end
                    % Collect all data together
                    rxntemp = CurrxnNum{mc};
                    genetemp = genenms1{1};
                    pathtemp = pathnames{cp};
                    if ~isempty(Allrxns{1})
                        if ~any(ismember(Allrxns,rxntemp))
                            Allgenes{yct} = genetemp;
                            Allrxns{yct} = rxntemp;
                            Allbigg{yct} =  biggtemp;
                            Allpath{yct} = pathnames{cp};
                        else
                            Aloci = find(ismember(Allrxns,rxntemp));
                            Allgenes{yct} = strjoin(unique({Allgenes{Aloci},...
                                genetemp}),'|');
                            Allpath{yct} = strjoin(unique({Allpath{Aloci},...
                                pathtemp}),';');
                            Allrxns{yct} = rxntemp;
                            Allbigg{yct} =  biggtemp;
                            Allgenes(Aloci) = []; Allpath(Aloci) = [];
                            Allrxns(Aloci) = []; Allbigg(Aloci) = [];
                            yct = yct - numel(Aloci);
                        end
                    else
                        Allgenes{yct} = genetemp;
                        Allrxns{yct} = rxntemp;
                        Allbigg{yct} =  biggtemp;
                        Allpath{yct} = pathnames{cp};
                    end
                    yct = yct + 1;
                end                
            end
        end
    end
end
for ct2 = 1:numel(Allgenes)
    temp1 = unique(strsplit(Allgenes{ct2},'|'));
    Allgenes{ct2} = strjoin(temp1,'|');
    temp2 = unique(strsplit(Allpath{ct2},';'));
    Allpath{ct2} = strjoin(temp2,';');
end
% Decomposing data on the basis of pathway code
AllgenesP = {[]}; AllrxnsP = {[]}; AllbiggP = {[]};
for ct = 1:numel(pathnames);
    pathloci = regexp(Allpath,pathnames{ct});
    pathloci = find(~cellfun('isempty', pathloci));
    AllgenesP{ct} = Allgenes(pathloci);
    AllrxnsP{ct} = Allrxns(pathloci);
    AllbiggP{ct} = Allbigg(pathloci);
end
Biostat = get(handles.Biographcheck,'Value');
if Biostat
    hSelectedObj=get(handles.Biographpanel, 'SelectedObject');
    BioType = get(hSelectedObj, 'Tag');
    Pathchoice = get(handles.Pathwaylistbox,'Value');
    gene2rxn(AllgenesP{Pathchoice},AllrxnsP{Pathchoice},AllbiggP{Pathchoice},BioType)
end


% Write data to an Excel file
Pth1 = which ('Precons.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Filename = [Pth,'\',Org,'.xlsx'];
if exist(Filename,'file')
    Txt3 = [Org,'.xlsx already exists in the directory.\n'];
    fprintf(Txt3)
    Txt4 = 'The file will be replaced with the new one\n';
    fprintf(Txt4)
    delete(Filename)  
end
Txt2 = 'Writing data to Xlsx file...\n';
fprintf(Txt2)
Alldat = cell(numel(Allgenes)+1,4);
Alldat{1,1} = 'Gene Entry'; Alldat{1,2} = 'KEGG'; Alldat{1,3} = 'BiGG';
Alldat{1,4} = 'Pathway';
for nt = 1:numel(Allgenes)
    Alldat{nt+1,1} = Allgenes{nt};
    Alldat{nt+1,2} = Allrxns{nt};
    Alldat{nt+1,3} = Allbigg{nt};
    Alldat{nt+1,4} = Allpath{nt};
end
xlswrite(Filename,Alldat)
clc
fprintf('\nDone!')


function gene2rxn(Allgenes,Allrxns,Allbigg,BioType)
AllgeneNO1 = strjoin(Allgenes,'|');
AllgeneNO2 = strsplit(AllgeneNO1,'|'); AllgeneNO2 = unique(AllgeneNO2);
AllgeneNO = numel(AllgeneNO2); % Number of all unique genes
% AllbiggU1 = unique(Allbigg);
AllbiggU = find(~cellfun('isempty', Allbigg));
AllbiggStr = Allbigg(AllbiggU);
AllgeneStr = num2str(1:numel(AllgeneNO2));
AllgeneStr = strsplit(AllgeneStr,' '); AllgeneStr = strcat('G',AllgeneStr);
AllbiggNostr = num2str(1:numel(AllbiggStr));
if ~isempty(AllbiggNostr)
    AllbiggNostr = strsplit(AllbiggNostr,' ');
    AllbiggNostr = strcat('B',AllbiggNostr);
end

if strcmp(BioType,'G2K2Bradiobutton') %================================
    GprS = zeros(AllgeneNO + numel(Allrxns) + numel(AllbiggU));
    for ct1 = 1:numel(Allrxns)
        currgenes = strsplit(Allgenes{ct1},'|');
        for ct2 = 1:numel(currgenes)
            currgeneloci = find(ismember(AllgeneNO2,currgenes{ct2}));
            GprS(currgeneloci,AllgeneNO+ct1) = 1;
        end
        if ~isempty(Allbigg{ct1})
            y_coor = find(ismember(AllbiggStr,Allbigg{ct1}));
            GprS(AllgeneNO+ct1,AllgeneNO+numel(Allrxns)+y_coor) = 1;    
        end
    end
    Allcum = [AllgeneStr,Allrxns,AllbiggNostr];
    GprS = sparse(GprS);
    bg1 = biograph(GprS,Allcum);
    dolayout(bg1);
    
    Tetta1 = 0:pi/(numel(AllgeneStr)-1):2*pi;
    Tetta2 = 0:pi/(numel(Allrxns)-1):2*pi;
    Tetta3 = 0:pi/(numel(AllbiggStr)-1):2*pi;
    R1 = 350; R2 = R1*1.5; R3 = R1*2;
    X1 = R1.*cos(Tetta1); Y1 = R1.*sin(Tetta1);
    X2 = R2.*cos(Tetta2); Y2 = R2.*sin(Tetta2);
    X3 = R3.*cos(Tetta3); Y3 = R3.*sin(Tetta3);
    for ct = 1:numel(Allcum)
        if ct<=numel(AllgeneStr)
            bg1.nodes(ct).Position = [X1(ct),Y1(ct)];
            set(bg1.nodes(ct),'Shape','circle')
            set(bg1.nodes(ct),'Label',AllgeneNO2{ct})
        elseif ct>numel(AllgeneStr) && ct<=(numel(Allrxns)+numel(AllgeneStr))
            bg1.nodes(ct).Position = [X2(ct-numel(AllgeneStr)),...
                Y2(ct-numel(AllgeneStr))];
            set(bg1.nodes(ct),'Shape','diamond')
        else
            bg1.nodes(ct).Position = [X3(ct-numel(AllgeneStr)-numel(Allrxns)),...
                Y3(ct-numel(AllgeneStr)-numel(Allrxns))];
            set(bg1.nodes(ct),'Label',AllbiggStr{ct-numel(AllgeneStr)-numel(Allrxns)})
        end
    end

elseif strcmp(BioType,'G2Bradiobutton') %====================================
    gene4biggRaw = Allgenes(AllbiggU);
    gene4bigg = strjoin(gene4biggRaw,'|');
    gene4bigg = strsplit(gene4bigg,'|'); gene4bigg = unique(gene4bigg);
    gene4biggNo = numel(gene4bigg); % Number of all unique genes
    AllgeneStr1 = num2str(1:numel(gene4bigg));
    AllgeneStr1 = strsplit(AllgeneStr1,' ');
    AllgeneStr1 = strcat('G',AllgeneStr1);
    GprS = zeros(gene4biggNo + numel(AllbiggU));
    for ct1 = 1:numel(gene4biggRaw)
        currgenes = strsplit(gene4biggRaw{ct1},'|');
        for ct2 = 1:numel(currgenes)
            currgeneloci = find(ismember(gene4bigg,currgenes{ct2}));
            GprS(currgeneloci,gene4biggNo+ct1) = 1;
        end
    end
    Allcum = [AllgeneStr1,AllbiggNostr];
    GprS = sparse(GprS);
    bg1 = biograph(GprS,Allcum);
    dolayout(bg1);
    Tetta1 = 0:pi/(numel(AllgeneStr1)-1):2*pi;
    Tetta2 = 0:pi/(numel(AllbiggNostr)-1):2*pi;
        R1 = 350; R2 = R1*1.5;
    X1 = R1.*cos(Tetta1); Y1 = R1.*sin(Tetta1);
    X2 = R2.*cos(Tetta2); Y2 = R2.*sin(Tetta2);
    for ct = 1:numel(Allcum)
        if ct<=numel(AllgeneStr1)
            bg1.nodes(ct).Position = [X1(ct),Y1(ct)];
            set(bg1.nodes(ct),'Shape','circle')
            set(bg1.nodes(ct),'Label',gene4bigg{ct})
        else
            bg1.nodes(ct).Position = [X2(ct-numel(AllgeneStr1)),...
                Y2(ct-numel(AllgeneStr1))];
            set(bg1.nodes(ct),'Label',AllbiggStr{ct-numel(AllgeneStr1)})
        end
    end
else %=====================================================================
    GprS = zeros(AllgeneNO + numel(Allrxns));
    for ct1 = 1:numel(Allrxns)
        currgenes = strsplit(Allgenes{ct1},'|');
        for ct2 = 1:numel(currgenes)
            currgeneloci = find(ismember(AllgeneNO2,currgenes{ct2}));
            GprS(currgeneloci,AllgeneNO+ct1) = 1;
        end
    end
    Allcum = [AllgeneStr,Allrxns];
    GprS = sparse(GprS);
    bg1 = biograph(GprS,Allcum);
    dolayout(bg1);
    Tetta1 = 0:pi/(numel(AllgeneStr)-1):2*pi;
    Tetta2 = 0:pi/(numel(Allrxns)-1):2*pi;
    R1 = 350; R2 = R1*1.5;
    X1 = R1.*cos(Tetta1); Y1 = R1.*sin(Tetta1);
    X2 = R2.*cos(Tetta2); Y2 = R2.*sin(Tetta2);
    for ct = 1:numel(Allcum)
        if ct<=numel(AllgeneStr)
            bg1.nodes(ct).Position = [X1(ct),Y1(ct)];
            set(bg1.nodes(ct),'Shape','circle')
            set(bg1.nodes(ct),'Label',AllgeneNO2{ct})
        else
            bg1.nodes(ct).Position = [X2(ct-numel(AllgeneStr)),...
                Y2(ct-numel(AllgeneStr))];
            set(bg1.nodes(ct),'Shape','diamond')
        end
    end
    
end

dolayout(bg1, 'Pathsonly', true);
view(bg1)