function [MapChoice, Outflx, RxnCds] = GetKegg(Bigg, Inflx, ModelName)
% GetKegg
% reads flux values obtained from COBRA optimization, and uses the 
% input BiGG reactions to retreive the equivalent KEGG IDs obtained from 
% Bigg2Kegg function. 
% Inputs:
% Bigg = A cell array of BiGG reaction IDs (can be obtained from readCbModel
%        function of COBRA toolbox) for which flux values to be displayed
%        on KEGG maps.
% Inflx = A mXn matrix containing input flux data, which m corresponds to
%         BiGG IDs and n corresponds to time-series values.
% ModelName = Model description according to BiGG database (e.g. RECON1)
% Output:
% MapChoice = A cell array of KEGG pathway identifiers chosen by the user
%             via a simple GUI. These set of KEGG pathway IDs will be used
%             as the input of KeggDraw function.
% RxnCds = KEGG reaction identifiers equivalent to input Bigg IDs.
% Outflx = Similar to Inflx, with the difference that m corresponds to
%          RxnCds (Note: if KEGG rxn IDs for some Bigg elements cannot be 
%          found, these IDs should be trimmed).
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Nov. 2015

% Retrieve BiGG2KEGG rxn identifiers from Bigg2Kegg function---------------
Pth1 = which('GetKegg.m');
tind = find(Pth1=='\',1,'last');
Pth2 = Pth1(1:tind-1);
Nme = [ModelName,'KEGG.mat'];
TargetPth = fullfile(Pth2,'BIGG2KEGG',Nme);
load(TargetPth);
Brxns1 = B2Kegg.B; Krxns1 = B2Kegg.K;
% -------------------------------------------------------------------------
% Find common rxns between Brxns and Bigg
% First, seperate multi-step rxns------------------------------------------
MultiB1 = Brxns1(find(strcmp(Brxns1,'MULTIR'))+1:end);
MultiK1 = Krxns1(find(strcmp(Brxns1,'MULTIR'))+1:end);
MultiK = MultiK1(ismember(MultiB1, Bigg));
MultiB = MultiB1(ismember(MultiB1, Bigg));
[~,Fx22] = ismember(MultiB1, Bigg);
Fx22(Fx22==0) = [];
MultiFlx = Inflx(Fx22,:);
setappdata(0,'MultiB',MultiB); % For later use in KeggDraw
setappdata(0,'MultiK',MultiK);
setappdata(0,'MultiFlx',MultiFlx);
% -------------------------------------------------------------------------
Brxns = Brxns1(1:find(strcmp(Brxns1,'MULTIR'))-1);
Krxns = Krxns1(1:find(strcmp(Brxns1,'MULTIR'))-1);
RxnCdsI = Krxns(ismember(Brxns, Bigg));
BiGGI = Brxns(ismember(Brxns, Bigg));
[~,Fx2] = ismember(Brxns, Bigg);
Fx2(Fx2==0) = [];
FlxI = Inflx(Fx2,:);
% If any, ask user to fill missed KEGG rxn IDs 
if any(~ismember(Bigg,Brxns))
    disp('Some KEGG rxns cannot be found for input BiGG IDs')
    disp('Insert them manually')
    disp('NOTE: Exchange reactions are not displayed!')
    NotData = Bigg(~ismember(Bigg,Brxns));
    % Exclude Exchange reactions : EX_w\*
    NoDataIdf = regexp(NotData,'\<EX_\w*');
    NoDataIdf1 = ~cellfun('isempty', NoDataIdf);
    NotData(NoDataIdf1) = [];
    %------------------------------------
    setappdata(0,'NotData',NotData)
    waitfor(GetKeggTable)
    ModfdTab = getappdata(0,'Final');
    FindUndfnd = ModfdTab(:,2);
    BiGGUndfnd =  ModfdTab(:,1);
    DfinedLoci = ~strcmpi(FindUndfnd,'Unknown');
    RxnCdsII = FindUndfnd(DfinedLoci);
    BiGGII = BiGGUndfnd(DfinedLoci);
    Flxs = Inflx(~ismember(Bigg,Brxns),:);
    FlxII = Flxs(DfinedLoci,:);
    RxnCds = [RxnCdsI';RxnCdsII];
    BiGG4KeggDraw = [BiGGI';BiGGII];
    setappdata(0,'BiGG4KeggDraw',BiGG4KeggDraw);
    Outflx = [FlxI;FlxII];
else
    RxnCds = RxnCdsI;
    Outflx = FlxI;
end
%--------------------------------------------------------------------------
% Find all KEGG pathways for RxnCds
Fileid1 = fopen('rxn2map.txt','r');
rxn2map = textscan(Fileid1,'%s %s');
RawRxns = rxn2map{1};
RawMaps = rxn2map{2};
[Loci1,~] = regexp(RawMaps, 'path:rn');
Loci2 = ~cellfun('isempty', Loci1);
RawMaps(Loci2) = [];
RawRxns(Loci2) = [];
Maps = strrep(RawMaps,'path:map','');
Rxns = strrep(RawRxns,'rn:','');
RxnTemp = Rxns(ismember(Rxns,RxnCds));
MapTemp = Maps(ismember(Rxns,RxnCds));
count=1;
RxnList = {0}; MapList = {0}; SortTemp = (0);
RmvMaps = {'01100','01110','01120','01130','01200','01210','01212',...
    '01220','01230','00121'};
while numel(MapTemp)
    Loci3 = ismember(MapTemp,MapTemp(1));
    if ~ismember(MapTemp(1),RmvMaps)
        RxnList{count} = RxnTemp(Loci3);
        SortTemp(count) = numel(RxnList{count});
        MapList(count) = MapTemp(1);
        count=count+1;
    end
    RxnTemp(Loci3)=[];
    MapTemp(Loci3)=[];
end
% Sort MapList based on rxns in RxnList
[~,ST] = sort(SortTemp,'descend');
RxnList = RxnList(ST);
MapList = MapList(ST);
% Add pathway names to map codes of MapList
Fileid2 = fopen('KEGGmaps.txt','r');
TempMaps = textscan(Fileid2,'%s %[^\n]');
TempMaps1 = cell(size(TempMaps{1},1),1);
for count = 1:size(TempMaps{1},1)
    TempMaps1{count}=[TempMaps{1}{count},' ',TempMaps{2}{count}];
end
TempMaps1 = strrep(TempMaps1,'path:map','');
MapNames2 = cell(numel(MapList),1);
for count = 1:numel(MapList)
    Loci4 = regexp(TempMaps1,MapList{count});
    MapNames2{count} = TempMaps1{~cellfun('isempty', Loci4)};
end

setappdata(0,'MapNames2',MapNames2)
setappdata(0,'RxnList',RxnList)
waitfor(GetKeggMaps)
MapChoiceLoci = getappdata(0,'AllChoices');
MapChoice = MapList(MapChoiceLoci);
% Remove previous data in appdata------------------------------------------
AT=getappdata(0);
AT1 = fieldnames(AT);
for co = 1:length(AT1)
    if strcmp(AT1{co},'BiGGModel') || strcmp(AT1{co},'Mdl') || ...
            strcmp(AT1{co},'CRxns')
        rmappdata(0,AT1{co});
    end
end