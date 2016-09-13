function [ArX,ArY,Cntrs,Director] = ArrowheadProcess(Xin1,Yin1,min_x,min_y,...
    idx,basemap,JumpData,RefOverlapData,ID,flx_idx)
% ArrowheadProcess
% is a subfunction of NetDraw for identifiying the
% appropriate coordinates for drawing arrowheads on metabolic maps.
% Inputs:
% Xin1, Yin1: Reaction coordinates
% min_x, min_y: Minimum of coordinates on the current customized map
% basemap: KGML data of KEGG global metabolic pathway
% JumpData: Details (ID) for overlapping flux carrying reactions
% RefOverlapData: Details (ID) for all reactions in the current customized
% map.
% ID: Reaction directionality based on either KEGG or BiGG (string)
% flx_idx: Details (ID) for all flux carrying reactions
% 
% Outputs:
% ArX,ArY: Appropriate coordinates for overlaying arrowheads on map.
% Cntrs: Centers' coordinates corresponding to ArX and ArY for MapTrimmer
% or MapTrimmerF functions.
% Director: Reaction directionality for arrowhead function.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016
%--------------------------------------------------------------------------

ArX = ({}); ArY = ArX; Cntrs = ArX; Jumpidx = []; Director = ArX;

if ~isempty(RefOverlapData)
    [~,Idx1,~] = unique(RefOverlapData);
    Chck1 = setdiff(1:size(RefOverlapData,1),Idx1);
    Jumpidx = JumpData.Id;
    RefJumpidx = Jumpidx;
    for m2 = 1:numel(Chck1)
        WhrC = find(ismember(RefOverlapData,RefOverlapData{Chck1(m2)}));
        RedundantIdx = idx(WhrC);
        if any(ismember(RedundantIdx,RefJumpidx))
            Addme = setdiff(RedundantIdx,RefJumpidx);
            for m3 = 1:numel(Addme)
                CurrLen = numel(Jumpidx);
                Jumpidx{CurrLen+1} = Addme{m3};
            end
        end
    end
    Jumpidx = unique(Jumpidx);
end
if strcmp(ID,'bigg') % Rxns directionality based on BiGG-------------------
    Bgcoord = getBiGGrxns(idx,basemap);
end
% -------------------------------------------------------------------------
for m1 = 1:numel(idx)
    if ~isempty(Jumpidx)
        if any(ismember(idx{m1},Jumpidx))
            ArX{m1} = ''; ArY{m1} = ''; Cntrs{m1} = ''; Director{m1} = '';
            continue
        end
    end
    [~,in1] = intersect(basemap.rc.rxnid,str2double(idx{m1}(7:end)));
    Xin = Xin1{m1};
    Yin = Yin1{m1};
    if isempty(in1)
        ArX{m1} = ''; ArY{m1} = ''; Cntrs{m1} = ''; Director{m1} = '';
        continue
    end
    if strcmp(ID,'bigg')
        if strcmp(Bgcoord{m1},'None')
            continue
        end
        if ~isempty(Bgcoord{m1})
            ProdCoords = Bgcoord{m1}.P; SubCoords = Bgcoord{m1}.S;
        else
            [ProdCoords,SubCoords] = getCoords4KEGG(basemap,in1);
        end
    else
        [ProdCoords,SubCoords] = getCoords4KEGG(basemap,in1);
    end

    if isempty(SubCoords)
        MixCoords = ProdCoords;
        Director{m1} = ones(size(MixCoords,1),1);
    else
        MixCoords = [ProdCoords;SubCoords];
        Director{m1} = [ones(size(ProdCoords,1),1);-1.*ones(size(SubCoords,1),1)];
    end
    MixCoords(:,1) = MixCoords(:,1) - min_x - 18;
    MixCoords(:,2) = MixCoords(:,2) - min_y;
    X = (0); Y =(0); Cntr =(0);
    for cnt = 1:size(MixCoords,1)
         DistancesToCpd = hypot(Xin - MixCoords(cnt,1), Yin - MixCoords(cnt,2));
         ind3 = find(DistancesToCpd >= 0.7*MixCoords(cnt,4)-1 & DistancesToCpd <= MixCoords(cnt,4));
         if isempty(ind3)
             [~,ind3] = min(abs(DistancesToCpd));
         end
         InrangeX = Xin(ind3); InrangeY = Yin(ind3);
         DistancesToCpd1 = hypot(InrangeX - MixCoords(cnt,1), InrangeY - MixCoords(cnt,2));
         [~,ind4] = min(abs(DistancesToCpd1));
          X(cnt,2) = InrangeX(ind4); Y(cnt,2) = InrangeY(ind4);  
         ind5 = find(DistancesToCpd >= MixCoords(cnt,4)+7);
         OutrangeX = Xin(ind5); OutrangeY = Yin(ind5);
         DistancesToPint1_temp = hypot(OutrangeX - X(cnt,2), OutrangeY - Y(cnt,2));
         [~,ind6] = min(abs(DistancesToPint1_temp));
         X(cnt,1) = OutrangeX(ind6); Y(cnt,1) = OutrangeY(ind6);
         Cntr(cnt,1) = MixCoords(cnt,1); Cntr(cnt,2) = MixCoords(cnt,2);
    end
    ArX{m1} = X; ArY{m1} = Y; Cntrs{m1} = Cntr;
end
if ~isempty(flx_idx)
    [ArX,ArY,Cntrs,Director] = PostMod(ArX,ArY,Cntrs,Director,...
        idx,RefOverlapData,flx_idx);
end
% Subfunction =============================================================
function AllCoords = getBiGGrxns(idx,basemap)
% Read BiGG model ---------------------------------------------------------
BiGGModel = getappdata(0,'BiGGModel'); Mdl = getappdata(0,'Mdl');
if isempty(BiGGModel) || isempty(Mdl)
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
    Chk = which([pname(1:strfind(pname,'.')-1),'KEGG.mat']);
    if isempty(Chk)
        msgbox({[pname(1:strfind(pname,'.')-1),'KEGG.mat',' cannot be found!'];...
            'Make sure the file is present in BiGG2KEGG folder'},'Error','error');
        return
    end
    Mdl = load(which([pname(1:strfind(pname,'.')-1),'KEGG.mat']));
    setappdata(0,'Mdl',Mdl);
end

Metbigg = BiGGModel.Metbigg;
Metkegg = BiGGModel.Metkegg;
CbModel = BiGGModel.CbModel;
% -------------------------------------------------------------------------
modelmets = CbModel.mets;
modelmets = strrep(modelmets,'[','_'); modelmets = strrep(modelmets,']','');
modelmets = regexp(modelmets,'\w*(?=_\w{1}\>)','match'); % Remove compartments
modelmets = [modelmets{:}];
modelrxns = CbModel.rxns; S = CbModel.S; K = Mdl.B2Kegg.K; B = Mdl.B2Kegg.B;
Rev=CbModel.rev;

for m1 = 1:numel(idx)
    Isrxn = B(ismember(K,idx{m1}(1:6))); 
    if ~isempty(Isrxn)% If rxn has a corresponding BiGG ID
        IsMdlRxn = find(ismember(modelrxns,Isrxn(1)));
        RevState = Rev(IsMdlRxn);
        if ~isempty(IsMdlRxn) % If corresponding BiGG rxn exists in the model
            [~,in1] = intersect(basemap.rc.rxnid,str2double(idx{m1}(7:end)));
            if isempty(in1)
                AllCoords{m1} = 'None';
                continue
            end
            MetIdx = find(S(:,IsMdlRxn));
            Metsign = S(MetIdx,IsMdlRxn); Metsign(Metsign>0)=1;Metsign(Metsign<0)=-1;
            ctr = 1; foundMets = ({}); MetsignEx = (0);
            for m2 = 1:numel(MetIdx)
                Whrmt = find(ismember(Metbigg,modelmets(MetIdx(m2))));
                for m3 = 1:numel(Whrmt)
                    foundMets{ctr} = Metkegg{Whrmt(m3)};
                    MetsignEx(ctr) = Metsign(m2);
                    ctr = ctr + 1;
                end
            end
            rxnSubs = foundMets(MetsignEx<0);
            rxnProds = foundMets(MetsignEx>0);
            ProdNames = basemap.rc.prod{in1};
            SubNames = basemap.rc.sub{in1};
            if RevState %Reversible rxn
                if any(ismember(rxnProds,SubNames)) || any(ismember(rxnSubs,ProdNames))
                    metCoords = getCoords(basemap,in1,'Rev_back');
                else
                    metCoords = getCoords(basemap,in1,'Rev_forward');
                end
            else % Irreversible rxn
                if any(ismember(rxnProds,SubNames)) || any(ismember(rxnSubs,ProdNames))
                    metCoords = getCoords(basemap,in1,'Irr_back');
                elseif any(ismember(rxnSubs,SubNames)) || any(ismember(rxnProds,ProdNames))
                    metCoords = getCoords(basemap,in1,'Irr_forward');
                else % All KEGG rxns are reversible!
                    metCoords = getCoords(basemap,in1,'Rev_forward');
                end
            end   
            AllCoords{m1}.P = metCoords.P; AllCoords{m1}.S = metCoords.S;
        end
    else
        AllCoords{m1} = [];
    end
end
% -------------------------------------------------------------------------
function metCoords = getCoords(basemap,indx,Stat)
switch Stat
    case 'Rev_forward' %////////////////////////////////
        ProdId = basemap.rc.prodid{indx};
        SubsId = basemap.rc.subid{indx};
        [~,in1] = intersect(basemap.c.cpdid,ProdId);
        [~,in2] = intersect(basemap.c.cpdid,SubsId);
        ProdCoords = basemap.c.cpdxy(in1,:);
        SubCoords = basemap.c.cpdxy(in2,:);
        ProdCoords1 = []; SubCoords1 = [];
        if isempty(in1)
            [~,in1] = intersect(basemap.c.glid,ProdId);
            ProdCoords = basemap.c.glxy(in1,:);
        elseif numel(in1) < numel(ProdId)
            [~,in1] = intersect(basemap.c.glid,ProdId);
            ProdCoords1 = basemap.c.glxy(in1,:);
        end
        if isempty(in2)
            [~,in2] = intersect(basemap.c.glid,SubsId);
            SubCoords = basemap.c.glxy(in2,:);
        elseif numel(in2) < numel(SubsId)
            [~,in2] = intersect(basemap.c.glid,SubsId);
            SubCoords1 = basemap.c.glxy(in2,:);
        end
        ProdCoords = [ProdCoords;ProdCoords1];
        SubCoords = [SubCoords;SubCoords1];
    case 'Rev_back' %////////////////////////////////
        ProdId = basemap.rc.prodid{indx};
        SubsId = basemap.rc.subid{indx};
        [~,in1] = intersect(basemap.c.cpdid,ProdId);
        [~,in2] = intersect(basemap.c.cpdid,SubsId);
        ProdCoords = basemap.c.cpdxy(in1,:);
        SubCoords = basemap.c.cpdxy(in2,:);
        ProdCoords1 = []; SubCoords1 = [];
        if isempty(in1)
            [~,in1] = intersect(basemap.c.glid,ProdId);
            ProdCoords = basemap.c.glxy(in1,:);
        elseif numel(in1) < numel(ProdId)
            [~,in1] = intersect(basemap.c.glid,ProdId);
            ProdCoords1 = basemap.c.glxy(in1,:);
        end
        if isempty(in2)
            [~,in2] = intersect(basemap.c.glid,SubsId);
            SubCoords = basemap.c.glxy(in2,:);
        elseif numel(in2) < numel(SubsId)
            [~,in2] = intersect(basemap.c.glid,SubsId);
            SubCoords1 = basemap.c.glxy(in2,:);
        end
        ProdCoordsR = [ProdCoords;ProdCoords1];
        SubCoordsR = [SubCoords;SubCoords1];
        SubCoords = ProdCoordsR; ProdCoords = SubCoordsR; % Switch coords
    case 'Irr_forward' %////////////////////////////
        ProdId = basemap.rc.prodid{indx};
        [~,in2] = intersect(basemap.c.cpdid,ProdId);
        ProdCoords = basemap.c.cpdxy(in2,:);
        ProdCoords1 = [];
        if isempty(in2)
            [~,in2] = intersect(basemap.c.glid,ProdId);
            ProdCoords = basemap.c.glxy(in2,:);
        elseif numel(in2) < numel(ProdId)
            [~,in2] = intersect(basemap.c.glid,ProdId);
            ProdCoords1 = basemap.c.glxy(in2,:);
        end
        ProdCoords = [ProdCoords;ProdCoords1];
        SubCoords = [];
    case 'Irr_back' %//////////////////////////////
        ProdId = basemap.rc.subid{indx}; % Switch sub/prod
        [~,in2] = intersect(basemap.c.cpdid,ProdId);
        ProdCoords = basemap.c.cpdxy(in2,:);
        ProdCoords1 = [];
        if isempty(in2)
            [~,in2] = intersect(basemap.c.glid,ProdId);
            ProdCoords = basemap.c.glxy(in2,:);
        elseif numel(in2) < numel(ProdId)
            [~,in2] = intersect(basemap.c.glid,ProdId);
            ProdCoords1 = basemap.c.glxy(in2,:);
        end
        ProdCoords = [ProdCoords;ProdCoords1]; SubCoords = [];
end
metCoords.P = ProdCoords; metCoords.S = SubCoords;
%--------------------------------------------------------------------------
function [ProdCoords,SubCoords] = getCoords4KEGG(basemap,indx)
ProdId = basemap.rc.prodid{indx};
SubsId = basemap.rc.subid{indx};
[~,in1] = intersect(basemap.c.cpdid,ProdId);
[~,in2] = intersect(basemap.c.cpdid,SubsId);
ProdCoords = basemap.c.cpdxy(in1,:);
SubCoords = basemap.c.cpdxy(in2,:);
ProdCoords1 = []; SubCoords1 = [];
if isempty(in1)
    [~,in1] = intersect(basemap.c.glid,ProdId);
    ProdCoords = basemap.c.glxy(in1,:);
elseif numel(in1) < numel(ProdId)
    [~,in1] = intersect(basemap.c.glid,ProdId);
    ProdCoords1 = basemap.c.glxy(in1,:);
end
if isempty(in2)
    [~,in2] = intersect(basemap.c.glid,SubsId);
    SubCoords = basemap.c.glxy(in2,:);
elseif numel(in2) < numel(SubsId)
    [~,in2] = intersect(basemap.c.glid,SubsId);
    SubCoords1 = basemap.c.glxy(in2,:);
end
ProdCoords = [ProdCoords;ProdCoords1];
SubCoords = [SubCoords;SubCoords1];
%--------------------------------------------------------------------------
function [ArX,ArY,Cntrs,Director] = PostMod(ArX,ArY,Cntrs,Director,...
    idx,RefOverlapData,flx_idx)

[~,Idx1,~] = unique(RefOverlapData);
Chck1 = setdiff(1:size(RefOverlapData,1),Idx1);
for ctr = 1:numel(Chck1)
     WhrC = find(ismember(RefOverlapData,RefOverlapData{Chck1(ctr)}));
     if any(ismember(idx(WhrC),flx_idx)) % Overlap with flux carrying rxns
         WhchrxnID = find(ismember(idx(WhrC),flx_idx));
         WhchX = ArX(WhrC(WhchrxnID)); sizeX = (0);
         WhchY = ArY(WhrC(WhchrxnID));
         WhchC = Cntrs(WhrC(WhchrxnID));
         WhchD = Director(WhrC(WhchrxnID));
         for i1 = 1:numel(WhchX)
             sizeX(i1) = size(WhchX{i1},1);
         end
         [~,maxId] = max(sizeX);
         for i1 = 1:numel(WhrC)
             ArX{WhrC(i1)} = WhchX{maxId};
             ArY{WhrC(i1)} = WhchY{maxId};
             Cntrs{WhrC(i1)} = WhchC{maxId};
             Director{WhrC(i1)} = WhchD{maxId};
         end 
     end    
end