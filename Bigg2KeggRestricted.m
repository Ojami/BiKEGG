function [BiGGID,AllRxn] = Bigg2KeggRestricted(D,Metkegg,rxnData,cpds,cpd2rxn)
% Bigg2KeggRestricted 
% is a subfunction of Bigg2Kegg and employs restricted
% conditions for identifying reaction correspondences.
% Full documentation is available in the user manual.
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Oct. 2015

load('UniModelKEGG.mat')
UniB = B2Kegg.B; UniK = B2Kegg.K;
Umulti = find(ismember(UniB,'MULTIR'));
UniB = UniB(1:Umulti-1); UniK = UniK(1:Umulti-1);
clear B2Kegg
Fileid12 = fopen('PrunedRxns.txt','r');
Prunedata12 = textscan(Fileid12,'%s %s %d %[^\n]');
Testpruned = Prunedata12{1};
fclose(Fileid12);
RxnNames = D.rxns;
Rxns = D.S;
Keggrm = (0);
BiGGID = {0};
AllRxn = {0};
Uctr = 1;
for rn = 1:size(Rxns,2) % for all rxns
    fprintf('Rxn %d of %d\n',rn,size(Rxns,2))
    
    % ---------------------------------------------------------------------
    %if current reaction exists in pruned reactions, or UniModelKEGG there 
    %is no need to further search for potential equivalent reactions.
    if any(ismember(Testpruned,RxnNames{rn}))
        continue
    end
    if any(ismember(UniB,RxnNames{rn}))
        UniBloci = find(ismember(UniB,RxnNames{rn}));
        for u1 = 1:numel(UniBloci)
            BiGGID{Uctr} = UniB{UniBloci(u1)};
            AllRxn{Uctr} = UniK{UniBloci(u1)};
            Uctr = Uctr + 1;
        end    
        continue
    end
    %----------------------------------------------------------------------
    Prd = find(Rxns(:,rn)); % Find substrates and products for each rxn
    CheckTransport = [Metkegg{Prd}]; % Exclude transport rxns
    if (numel(CheckTransport)/numel(unique(CheckTransport)) == 2)
        continue
    end
    % 1st condition: All reactants and products have metKeggIDs.
    % 2nd condition: Exclude all one element rxns
    if (~sum(strcmp (Metkegg(Prd),'')) && (length(Prd) > 1))
        Cbs = MetCombs(Metkegg,Prd); % All combinations (permutations) of
        % KEGG compounds, regarding the fact that there are more than one
        % KEGG cpd for some metabolites in the model.
        for Outctr = 1:size(Cbs,1)
            % Contains all KEGG cpds in a certain rxn
            TempMet = unique(Cbs(Outctr,:));
            % Ignore transport reactions with water molecule: c00001
            if (numel(TempMet) == 2) && sum(strcmp(TempMet,'C00001'))
                continue
            end
            if length(TempMet) > 1 % Exclude all exchange rxns
                KeggrmMethod = 'off'; %Offline
                switch KeggrmMethod
                    case 'off'
                        Id1 = ismember(cpds,TempMet);
                        Keggrm = str2num(rxnData(Id1,5:end));
                    case 'on'
                        Keggrm = KeggApiMatch (TempMet);
                end
                Uq = unique(Keggrm);
                N = histc(Keggrm, Uq);
                if isempty(max(N))
                    continue
                end
                if (max(N) == numel(TempMet))
                    UqTemp = Uq(N==max(N));
                    for count = 1:numel(UqTemp)
                        UqTemp1{count} = strcat('rn:R',repmat('0',...
                            [1,4- floor(log10(UqTemp(count)+eps))])...
                            ,num2str(UqTemp(count)));
                    end
                    if numel(UqTemp) > 1
                        for count = 1:length(UqTemp1)
                            UqTempSum(count) = sum(ismember(cpd2rxn{2},...
                                UqTemp1{count}));
                        end
                        if min(UqTempSum)==numel(TempMet)
                            [~,In1] = min(UqTempSum);
                            AllRxn{Uctr} = UqTemp1{In1}(4:end);
                            BiGGID {Uctr} = D.rxns{rn};
                            Uctr = Uctr + 1;
                        end
                    else
                        UqTempSum = sum(ismember(cpd2rxn{2},UqTemp1{1}));
                        if min(UqTempSum)==numel(TempMet)
                            AllRxn{Uctr} = UqTemp1{1}(4:end);
                            BiGGID {Uctr} = D.rxns{rn};
                            Uctr = Uctr + 1;
                        end
                        
                    end
                    UqTempSum = 0;
                    UqTemp1 = {0};
                    
                end
            end
            Keggrm = 0;
        end
    end
end
% Identify BiGG rxns for which there is more than one KEGG rxns, and prune
% BiGGIDs and AllRxns
[BiGGID,AllRxn,MultiBID,MultiRxn] = RxnPrune(BiGGID, AllRxn, RxnNames);
AllRxnTemp = {[]}; BiGGIDTemp = {[]}; Tctr = 1;
if any ([MultiBID{:}])
    for ct = 1:numel(MultiBID)
        for ci = 1:numel(MultiRxn{ct})
            BiGGIDTemp{Tctr} = MultiBID{ct};
            AllRxnTemp{Tctr} = MultiRxn{ct}{ci};
            Tctr = Tctr + 1;
        end
    end
    BiGGID = [BiGGID,BiGGIDTemp];
    AllRxn = [AllRxn,AllRxnTemp];
end



%% Subfunctions
function Cbs = MetCombs(Metkegg,Prd)
SizeCell = cell(1,numel(Prd));
for count = 1:numel(Prd)
    TempCounter = numel(Metkegg{Prd(count)});
    SizeCell{count} = 1:TempCounter;
end
MetkeggT = Metkegg(Prd);
MetkeggTemp = [MetkeggT{:}];
OutMet = combvec(SizeCell{:}).';
for count = 1:size(OutMet,2)-1 % Corresponds to linearized Metkegg(Prd)
    OutMet(:,count+1)=OutMet(:,count+1)+max(OutMet(:,count));
end
CbsT = MetkeggTemp(OutMet);
if any(strcmp(CbsT(:),'C00080'))
    Cbs1 = CbsT;
    [N1,N2] = find(strcmp(Cbs1,'C00080'));
    LenCbs = 1:numel(Cbs1);
    for count = 1:numel(N1)
        NonH2Cpd = find(~ismember(LenCbs,N2(count)));
        Cbs1{N1(count),N2(count)} = CbsT{N1(count),NonH2Cpd(1)};
    end
    Cbs = [CbsT;Cbs1];
else 
    Cbs1 = CbsT;
    Cbs2 = Cbs1;
    for count = 1:size(CbsT,1)
        Cbs1{count,size(CbsT,2)+1} = 'C00080'; % Add H+
        % Keep dimensional consistency:
        Cbs2{count,size(CbsT,2)+1} = Cbs2{count,size(CbsT,2)};
    end
    Cbs = [Cbs1;Cbs2];
end
%--------------------------------------------------------------------------
function [Bg,Ar,Bg2,Ar2] = RxnPrune(BiGGID, AllRxn, RxnNames)
Fileid = fopen('PrunedRxns.txt','r');
PrAll = textscan(Fileid,'%s %s %d %[^\n]');
PrB = PrAll{3};
PrBY = find(PrB); % To be added
PrBN = ~PrB; % To be removed
TBrxns = PrAll{1}; TKrxns = PrAll{2};
BrxnY = TBrxns(PrBY); KrxnY = TKrxns(PrBY);
BrxnN = TBrxns(PrBN);
% Which BiGGID elements are in BrxnY?
BiGGInBrxnY = BiGGID(ismember(BiGGID,BrxnY)); 
KeggOfBrxnY = KrxnY(ismember(BrxnY,BiGGID));
BrxnYTemp = BrxnY(ismember(BrxnY,BiGGID));
SharedRxns = cell(1,numel(BiGGInBrxnY));
for count = 1:numel(BrxnYTemp)
    T= ismember(BiGGInBrxnY,BrxnYTemp{count});
    SharedRxns(T) = KeggOfBrxnY(count);
end
AllRxn(ismember(BiGGID,BrxnY)) = SharedRxns;
% Which BrxnY elements are not in BiGGID and should be added separately?
NotInBiGG = BrxnY(~ismember(BrxnY,BiGGID));
RxnNot = KrxnY(~ismember(BrxnY,BiGGID));
% Are NotInBiGG elements exist in the current model (Note: BrxnY contains
% universal BiGG reactions) ?
YaLoci = find(ismember(NotInBiGG,RxnNames));
Yf1 = NotInBiGG(YaLoci);
AddedRxns = RxnNot(YaLoci);
ModLength = (1:numel(Yf1))+numel(BiGGID);
BiGGID(ModLength) = Yf1;
AllRxn(ModLength) = AddedRxns;
% Remove other found rxns in BiGGID which their absence is manually
% confirmed.
TobeRemovedRxns = ismember(BiGGID,BrxnN);
BiGGID(TobeRemovedRxns)=[]; AllRxn(TobeRemovedRxns)=[]; % Remove other rxns
% Find BiGG IDs with more than one rxn id
ct=1;cr=1;
Bg2={0};Ar2={0};Bg={0};Ar={0};
BiGGT = BiGGID; RxnT = AllRxn;
while numel(BiGGT)
    T = find(ismember(BiGGT,BiGGT(1))); R = RxnT(T);
    if numel(unique(R))~=1 % Reactions with more than one KEGG ID
        Bg2{ct} = BiGGT{1};
        Ar2{ct} = unique(R);      
        ct = ct+1;
    else
        Bg{cr}=BiGGT{1};
        Ar{cr}=cell2mat(unique(R));
        cr=cr+1;
    end
    BiGGT(T)=[]; %Remove
    RxnT(T)=[];
end
%--------------------------------------------------------------------------
function Keggrm = KeggApiMatch(TempMet)
ctr=1;
Keggrm=(0);
for mt = 1:length(TempMet)
    Metstr = urlread(['http://rest.kegg.jp/link/rn/',TempMet{mt}]);
    Id1 =  regexp(Metstr,'rn:R'); % All rxns containing a cpd
    if Id1 % KEGG rxn IDs are provided for this cpd
        for di = 1:length(Id1) % Extracts rxns
            Keggrm(ctr) = str2double(Metstr(Id1(di)+4:Id1(di)+9));
            ctr = ctr+1;
        end
    else
        break % Don't generate Metstr for other compounds
    end
end
% -------------------------------------------------------------------------
function Metkegg = ModelMoidfy (Metkegg,MetAbr)
% Read Metkegg data for RECON1.xml
Mdat = load('RECON1Metkegg.mat');
RECON1k = Mdat.RECON1meta.K;
RECON1m = Mdat.RECON1meta.B;
[Fx1,Fx2] = ismember(RECON1m,MetAbr);
Fx2(Fx2==0) = [];
% Metkegg(Fx2) = RECON1k(Fx1);
Fx3 = find(Fx1);
for i1 = 1:numel(Fx2)
    if numel(Metkegg{Fx2(i1)}) > numel(RECON1k{Fx3(i1)})
       Tempmet = Metkegg{Fx2(i1)};
       whrc = ismember(Tempmet,RECON1k{Fx3(i1)});
       if isempty(whrc)
        Metkegg{Fx2(i1)} = RECON1k(Fx3(i1));
        end
    else
        Metkegg{Fx2(i1)} = RECON1k{Fx3(i1)};
    end
end
% -------------------------------------------------------------------------
function [BiGGID,AllRxn] = AddLumpRxns(Steprxns,Insteprxns,BiGGID1,AllRxn1)
% ExcRxns contains rxns which need to be removed manually from Steprxns
ExcRxns = {'R03082','R01651','R01210'};
FindLocs = find(ismember(Steprxns,ExcRxns));
Steprxns(ismember(Steprxns,ExcRxns)) = [];
Insteprxns(FindLocs)=[];
PrevLen = numel(BiGGID1)+1;
BiGGID = BiGGID1; AllRxn = AllRxn1;
BiGGID{PrevLen} = 'MULTIR'; AllRxn{PrevLen} = 'MULTIR';
CurrLen = numel(BiGGID1)+2;
for count = 1:numel(Steprxns)
    if ~isempty(Insteprxns{count})
        if ~any(ismember(Insteprxns{count},AllRxn1)) %No common rxns
            BiggTemp = BiGGID1(ismember(AllRxn1,Steprxns{count}));
            for count1 = 1:numel(BiggTemp)
                for count2 = 1:numel(Insteprxns{count})
                    BiGGID{CurrLen} = BiggTemp{count1};
                    AllRxn{CurrLen} = Insteprxns{count}{count2};
                    CurrLen = CurrLen+1;
                end
            end
        else
            TempLoc = find(~ismember(Insteprxns{count},AllRxn1));
            BiggTemp = BiGGID1(ismember(AllRxn1,Steprxns{count}));
            for count1 = 1:numel(BiggTemp)
                for count2 = 1:numel(TempLoc)
                    BiGGID{CurrLen} = BiggTemp{count1};
                    AllRxn{CurrLen} = Insteprxns{count}{TempLoc(count2)};
                    CurrLen = CurrLen+1;
                end
            end
        end
    end
end