function [OverlayMe,Allx4Arrow,Ally4Arrow] = RxnAdjuster(PostData,Mdl,RefID)
% RxnAdjuster
% is a subfunction of MapAdjuster for displaying reaction
% details on created maps by NetDraw.
%
% Inputs:
% PostData: All necessary data from NetDraw
% Mdl: COBRA model for which in silico simulations are performed.
% RefID: A string; 'bigg' or 'kegg'.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016

Allx4Arrow = PostData.Allx4Arrow;
Ally4Arrow = PostData.Ally4Arrow;
field_intsect = PostData.field_intsect;
Idxs = PostData.field_intsect_overlay;
rxnCds = PostData.rxnCds;
[~,n2] = ismember(Idxs,field_intsect);
OrxnCds = rxnCds(n2);
ArXo = PostData.ArXo;
ArYo = PostData.ArYo;
RefOverlapData = PostData.RefOverlapData;
[~,Idx1,~] = unique(RefOverlapData);
Chck1 = setdiff(1:size(RefOverlapData,1),Idx1);

for m2 = 1:numel(Chck1)
    WhrC = ismember(RefOverlapData,RefOverlapData{Chck1(m2)});
    RedundantIdx = field_intsect(WhrC);
    if any(ismember(Idxs,RedundantIdx))
        if sum(ismember(Idxs,RedundantIdx)) < numel(RedundantIdx)
            Allx4Arrow(WhrC) = {''}; Ally4Arrow(WhrC) = {''}; 
            field_intsect(WhrC) = {''};
            WhrAdd = find(ismember(Idxs,RedundantIdx));
            for m3 = 1:numel(WhrAdd)
                Allx4Arrow{end+1} = ArXo{WhrAdd(m3)};
                Ally4Arrow{end+1} = ArYo{WhrAdd(m3)};
                field_intsect(end+1) = Idxs(WhrAdd(m3));
            end
        end      
    end   
end
emrmv = ~cellfun('isempty',field_intsect);
Allx4Arrow = Allx4Arrow(emrmv);
Ally4Arrow = Ally4Arrow(emrmv);
field_intsect = field_intsect(emrmv); NewRefData = ({});
for i1 =1:numel(Allx4Arrow)
    NewRefData{i1} = [num2str(numel(Allx4Arrow{i1})),',',...
        num2str(mean(Allx4Arrow{i1})),',',num2str(numel(Ally4Arrow{i1})),',',...
        num2str(mean(Ally4Arrow{i1}))];
end
[~,Idx1,~] = unique(NewRefData); clear WhrC
Chck1 = setdiff(1:size(NewRefData,2),Idx1);
[~,N2] = ismember(field_intsect,PostData.field_intsect);
switch RefID
    case 'kegg'
        OverlayMe = rxnCds(N2);
    case 'bigg'
        K1 = Mdl.B2Kegg.K;
        B1 = Mdl.B2Kegg.B;
        Chk = which('UniModelKEGG.mat');
        if isempty(Chk)
            msgbox({'UniModelKEGG.mat cannot be found!';...
                'Make sure the file is present in BiGG2KEGG folder'},'Error','error');
            return
        end
        CRxns = getappdata(0,'CRxns');
        OverlayMe = rxnCds(N2);
        if ~isempty(CRxns) % Correct for consistent reaction correspondences
            for cnter = 1:numel(CRxns.C)
                OverlayMe(ismember(OverlayMe,CRxns.C{cnter})) = ...
                    CRxns.O(cnter);
                OrxnCds(ismember(OrxnCds,CRxns.C{cnter})) = ...
                    CRxns.O(cnter);
            end
        end
        
        Uni = load(which('UniModelKEGG.mat'));
        K = Uni.B2Kegg.K;
        B = Uni.B2Kegg.B;
        for i1 = 1:numel(OverlayMe)
            if any(strcmp(OverlayMe{i1},OrxnCds))
                WhrK = find(ismember(K1,OverlayMe{i1}));
                if isempty(WhrK)
                    WhrK = find(ismember(K,OverlayMe{i1}));
                    Tmpchk = B(WhrK);
                else
                    Tmpchk = B1(WhrK);
                end
            else
                WhrK = find(ismember(K,OverlayMe{i1}));
                Tmpchk = B(WhrK);
            end
            if isempty(WhrK)
                OverlayMe(i1) = {''};
                continue
            end
            if numel(WhrK)>1                
                C1 = regexp(Tmpchk,'\w*p\>'); C1 = ~cellfun('isempty',C1);
                if sum(C1) < numel(Tmpchk)
                    Tmpchk(C1) = [];
                end
                C1 = regexp(Tmpchk,'\w*copy\d\>'); C1 = ~cellfun('isempty',C1);
                if sum(C1) < numel(Tmpchk)
                    Tmpchk(C1) = [];
                end
                C1 = regexp(Tmpchk,'\w*f\>'); C1 = ~cellfun('isempty',C1);
                if sum(C1) < numel(Tmpchk)
                    Tmpchk(C1) = [];
                end
                C1 = regexp(Tmpchk,'\w*m\>'); C1 = ~cellfun('isempty',C1);
                if sum(C1) < numel(Tmpchk)
                    Tmpchk(C1) = [];
                end
                C1 = regexp(Tmpchk,'ICDHhr'); C1 = ~cellfun('isempty',C1);
                if sum(C1) < numel(Tmpchk)
                    Tmpchk(C1) = [];
                end
            end
            OverlayMe(i1) = Tmpchk(1);
        end
end
Allx4Arrow = PostData.Allx4Arrow(N2);
Ally4Arrow = PostData.Ally4Arrow(N2);
overlayMetemp = OverlayMe;
for m2 = 1:numel(Chck1)
    WhrC = find(ismember(NewRefData,NewRefData{Chck1(m2)}));
    for m3 = 1:numel(WhrC)
        OverlayMe{WhrC(m3)} = strjoin(overlayMetemp(WhrC),',');
        if ~isempty(OverlayMe{WhrC(m3)})       
            if strcmp(OverlayMe{WhrC(m3)}(1),',')
                OverlayMe{WhrC(m3)}(1) = [];
            end
        end
        if ~isempty(OverlayMe{WhrC(m3)})
            if strcmp(OverlayMe{WhrC(m3)}(end),',')
                OverlayMe{WhrC(m3)}(end) = [];
            end
        end
    end
end
[~,N1] = unique(OverlayMe);
OverlayMe = OverlayMe(N1); Allx4Arrow = Allx4Arrow(N1); Ally4Arrow = Ally4Arrow(N1);