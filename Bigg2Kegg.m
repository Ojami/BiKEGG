function Bigg2Kegg(Idfier)
% BiKEGG
% generates corresponding reactions in KEGG and BiGG/HMR
% databases and saves the output as a *.mat file in BiGG2KEGG folder.
% Input:
%     - Idfier : a string of either 'bigg' or 'hmr'
% Full documentation can be found in the user manual
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Oct. 2015

Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
if strcmp(Idfier,'hmr')
    [pname,fname]=uigetfile('*.xlsx','Select the HMRdatabase xlsx file');
    [~,dat2]=xlsread([fname,pname],'RXNS');
    ind=dat2(:,1);
    ind1=ismember(ind,'#');
    ind2=find(ind1,3,'first');
    hdrs=dat2(1,:);
    % Reading the BIGG IDs:
    BIGG=dat2(2:ind2(3)-1,ismember(hdrs,'BIGG DATABASE ID')); 
    % Reading the KEGG rxn IDs:
    KEGG=dat2(2:ind2(3)-1,ismember(hdrs,'KEGG ID'));
    % Reading the rxn equations:
    Eqn=dat2(2:ind2(3)-1,ismember(hdrs,'EQUATION'));
    % Locating empty cells
    r1= ~ismember(KEGG,'');
    r2= ~ismember(BIGG,'');
    r3 = find(r1 & r2);
    KEGG=KEGG(r3);BIGG=BIGG(r3);Eqn=Eqn(r3);
    B2Kegg.Kegg=KEGG;B2Kegg.Bigg=BIGG;B2Kegg.Eqn=Eqn;
    fname = fullfile(Pth,'BiGG2KEGG','BiGG2KEGG_HMRbased-RECON1.mat');
    save (fname,'B2Kegg')
    return
end

if ~strcmp(Idfier,'bigg')
    msgbox('Unknown identifier: bigg/hmr','Error','error');
    return
end
% Read text file extracted from http://rest.kegg.jp/link/rn/cpd, which
% contains corresponding rxns for each metabolite.
Fileid = fopen(fullfile(Pth,'cpd2rxn.txt'),'r');
cpd2rxn = textscan(Fileid,'%s %[^\n]');
rxnData = cell2mat(cpd2rxn{2}); % KEGG rxn IDs
cpds = strrep(cpd2rxn{1},'cpd:',''); % KEGG Cpds
fclose(Fileid);
% Read the model file======================================================
[pname,fname1]=uigetfile({'*.*'},'Select the input model file (mat|xml)');
Testname=pname(strfind(pname,'.')+1:end);
if strcmp(Testname,'mat')
    D = importdata([fname1,pname]);
    keggloc1 = fieldnames(D);
    keggid = strcmpi(keggloc1,'metkeggid');
    if ~any(keggid)
        msgbox('There is no metKeggID field in the model!','Error','error');
        return
    end
    Metkegg = D.(keggloc1{keggid});
%     NOTE: following commented lines can be activated to retrieve the
%     Metkegg online form BiGG database
%     Mets = D.mets;
%     for counter = 1:length(Mets)
%         BiggMet = urlread(['http://bigg.ucsd.edu/api/v2/models/',...
%             D.description,'/metabolites/',Mets{counter}]);
%         KeggLoci = strfind(BiggMet,'KEGG');
%         if ~isempty(KeggLoci)
%             Metkegg{counter} = BiggMet(KeggLoci+62:KeggLoci+67);
%         else
%             Metkegg{counter} = '';
%         end
%     end
elseif strcmp(Testname,'xml')
    Idf = fopen([fname1,pname]);
    D = textscan(Idf,'%s');
    CbModel = readCbModel([fname1,pname]);
else
    msgbox('Unknown filetype!','Error','error');
    return
end
% =========================================================================
if strcmp(Testname,'xml')
    % Reading xml
    % Find listOfSpecies which contains metabolite annotations 
    SpeciesLoci1 = strfind(D{1},'listOfSpecies');
    SpeciesLoci = find(~cellfun('isempty', SpeciesLoci1));
    % Discard other data in the model:
    SpeciesData = D{1}(SpeciesLoci(1):SpeciesLoci(2));
    clear D
    % Find each metabolite-------------------------------------------------
    SpeciesStart1 = strfind(SpeciesData,'<species');
    SpeciesStart = find(~cellfun('isempty', SpeciesStart1));
    SpeciesEnd1 = strfind(SpeciesData,'</species>');
    SpeciesEnd = find(~cellfun('isempty', SpeciesEnd1));
    if numel(SpeciesStart) ~= numel(SpeciesEnd)
        msgbox('SBML format error:Species have wrong format!','Error',...
            'error');
        return
    end
    % ---------------------------------------------------------------------
    % Search for KEGG compounds for each BiGG metabolite. Note that each BiGG
    % metabolite may have more than one equivalent KEGG compund, and although
    % HMR and RECON2 (Thiele lab.) do not provide more than one KEGG compound
    % for each metabolite, here, we followed the format of BiGG models (i.e.
    % several equivalent KEGG compounds for a specific metabolite).
    Metkegg = {''}; % KEGG compounds will fill this cell!
    for count = 1:numel(SpeciesStart)
        SpeciesTempStr = SpeciesData(SpeciesStart(count):SpeciesEnd(count));
        SpeciesTempStr = [SpeciesTempStr{:}];
        KeggCpdLoci = regexp(SpeciesTempStr,...
            '(?<=http://identifiers.org/kegg.compound/)[CG]\d{5}','match');
        if isempty(KeggCpdLoci)
            Metkegg{count} = '';
            continue
        end
        for count1 = 1:numel(KeggCpdLoci)
            if strcmp(KeggCpdLoci{count1},'C00001') % H2O
                Metkegg{count}= {'C00001'};
                break
            elseif strcmp(KeggCpdLoci{count1},'C00027') % H2O2
                Metkegg{count}= {'C00027'};
                break
            else
                Metkegg{count}{count1}= KeggCpdLoci{count1};
            end
        end
    end
end
% Model modification: some metabolites have wrong/absent KEGG annotations.
% Therefore, on the basis of modified KEGG compound IDs for RECON1.xml,
% these metabolites will be added to Metkegg separately.
D = CbModel;
ModelName = D.description;
tind1 = find(ModelName=='\',1,'last');
ModelName1 = ModelName(tind1+1:end-4);
if ~strcmp(ModelName1,'RECON1')
    MetAbr = CbModel.mets;
    Metkegg = ModelMoidfy (Metkegg,MetAbr);
end

% Search for Metkegg components with numel > 1, compare their names to 
% BiGG metabolites, and remove other components with different names if a
% unique result is obtained. KEGGcpd.txt is used which contains all KEGG
% compounds. 
Fileid1 = fopen(fullfile(Pth,'KEGGcpd.txt'),'r');
Keggcpd = textscan(Fileid1,'%s %[^\n]');
CpdData = strrep(Keggcpd{1},'cpd:','');
CpdNames = regexp(Keggcpd{2},';','split');
BiGGmets = D.metNames;
for count = 1:numel(Metkegg)
    if all(ismember({'C01326','C00177'},Metkegg{count})) % For cyan
        continue
    end
% Second condition: L-carnitine definition discrepancy between KEGG and BiGG
    if numel(Metkegg{count})>1 && ~all(ismember({'C00318','C00487'},...
            Metkegg{count}))
        MatchLoci = find(ismember(CpdData,Metkegg{count}));
        for incount = 1:numel(MatchLoci)
            CpdNames{MatchLoci(incount)}=strtrim(CpdNames{MatchLoci(incount)});
            FoundL = find(ismember(BiGGmets,CpdNames{MatchLoci(incount)}));
            if FoundL
                % Check if found cpd participates in any rxn
                if any(strcmp(CpdData(MatchLoci(incount)),cpds))
                    Metkegg{count} = CpdData(MatchLoci(incount));
                    break
                end
            end
        end
    end
end

% Extract KEGG metabolites from BIGG API ==================================
Fileid3 = fopen(fullfile(Pth,'bigg_models_metabolites.txt'),'r');
biggmets = textscan(Fileid3,'%s %s %[^\n]');
fclose(Fileid3);
biggmetid = biggmets{1}; biggmetdes = biggmets{3}; clear biggmets
Metkegg1 = {''};
modelm = D.mets;
modelm=strrep(modelm,'[','_');modelm=strrep(modelm,']','');
for ct1 = 1:numel(modelm)
    tempdes = biggmetdes(strcmp(biggmetid,modelm{ct1}));
    temploci = regexp(tempdes,...
        'http://identifiers.org/kegg.compound/C\S{5}','match');
    temploci = temploci{1};
    if isempty(temploci)
        Metkegg1{ct1} = '';
    else
        for ct2 = 1:numel(temploci)
            regtemp = regexp(temploci{ct2},'C\S{5}','match');
            Metkegg1{ct1}{ct2} = regtemp{1};
        end
    end
end
% Consensus Metkegg: Metkegg1 from BiGG API and Metkegg from COBRA models
% are mixed to generate the most complete set of KEGG metabolites for a
% SBML model.
for cp1 = 1:numel(Metkegg1)
    if isempty(Metkegg1{cp1}) || (numel(Metkegg1{cp1}) ~= numel(Metkegg{cp1}))
        for cp2 = 1:numel(Metkegg{cp1})
            Metkegg1{cp1}{cp2} = Metkegg{cp1}{cp2};
        end
    end
end
clear Metkegg
Metkegg = Metkegg1;
% ==========================================================================
% Check the existence of metabolites participating in Bigg rxns and map the
% Bigg rxns to KEGG rxns if possible. To do so, metabolites in Bigg rxns are
% compared to KEGG metabolites, and KEGG metabolite IDs are located in KEGG
% databse to identify the KEGG rxn IDs they involve in KEGG databse. If all
% KEGG metabolites share a unique KEGG rxn ID, then that rxn ID corresponds
% to the Bigg rxn from which the metabolites were extracted (restricted 
% condition).
% To reduce the CPU time, rxns are first compared to rxns in PrunedRxns.txt
% and those in UniModelKEGG to exclude them from search process.
load('UniModelKEGG.mat')
UniB = B2Kegg.B; UniK = B2Kegg.K;
Umulti = find(ismember(UniB,'MULTIR'));
UniB = UniB(1:Umulti-1); UniK = UniK(1:Umulti-1);
clear B2Kegg
Fileid12 = fopen('PrunedRxns.txt','r');
Prunedata12 = textscan(Fileid12,'%s %s %d %[^\n]');
Testpruned = Prunedata12{1};
fclose(Fileid12);
RxnNames = D.rxns; % BiGG rxn names
Rxns = D.S; % It is assumed that this field is universal.
Keggrm = (0);
BiGGID = {0};
AllRxn = {0};
Uctr = 1; ChkTrans = 1; TransKegg = ({}); TransBiGG = ({});
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
    CheckTransport = [Metkegg{Prd}]; % For transport rxns
    if (numel(CheckTransport)/numel(unique(CheckTransport)) == 2) && ...
            numel(Prd) > 2
        continue
    end
    % 1st condition: All reactants and products have metKeggIDs.
    % 2nd condition: Exclude all one element rxns
    if (~sum(strcmp (Metkegg(Prd),'')) && (length(Prd) > 1))
        Cbs = MetCombs(Metkegg,Prd); % All combinations (permutations) of
        % KEGG compounds, regarding the fact that there is more than one
        % KEGG cpd for some metabolites in the model.
        for Outctr = 1:size(Cbs,1)
            % Contains all KEGG cpds in a certain rxn
            TempMet = unique(Cbs(Outctr,:));
            ChkTempMet = Cbs(Outctr,:);
            % Transport rxns with 2 cpds
            % 1st cond: finding a transport rxn; 3 is as a result of Cbs
            % structure.
            % 2nd cond: Transport reactions with 2 metabolites; the reason
            % for 1 is similar to the 1st cond.
            % 3rd cond: No redundancy
            if (numel(ChkTempMet)/numel(TempMet) == 3) ...
                    && (numel(TempMet) == 1) ...
                   && ~any(ismember(TransBiGG,RxnNames{rn}))
                fprintf('Transport rxn: %s\n',RxnNames{rn})
                TransKegg{ChkTrans}{1} = TempMet{1};
                TransKegg{ChkTrans}{2} = TempMet{1};
                TransBiGG{ChkTrans} = RxnNames{rn};
                ChkTrans = ChkTrans + 1;
                continue
            end
            if (numel(TempMet) == 2) && any(strcmp(TempMet,'C00080'))
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
                % Identify unique KEGG rxn ID
                Uq = unique(Keggrm);
                N = histc(Keggrm, Uq);
                if isempty(max(N))
                    continue
                end
                if (max(N) == numel(TempMet)) || ...
                        (max(N) == numel(TempMet)-1) % Found repeated rxns
                    UqTemp = Uq(N==max(N));
                    % Convert rxn codes to string
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
                        % Restrict condition: If only the number of metabolites
                        % in min(UqTempSum)reaction are identical to those
                        % of TempMet
                        if min(UqTempSum)==numel(TempMet) || ...
                               min(UqTempSum)+1 == numel(TempMet) || ...
                               min(UqTempSum)-1 == numel(TempMet)  
                            In1 = find(UqTempSum == min(UqTempSum(:)));
                            for Inctr = 1:numel(In1)
                                AllRxn{Uctr} = UqTemp1{In1(Inctr)}(4:end);
                                BiGGID {Uctr} = D.rxns{rn};
                                Uctr = Uctr + 1;
                            end
                        end
                    else
                        UqTempSum = sum(ismember(cpd2rxn{2},UqTemp1{1}));
                        if min(UqTempSum)==numel(TempMet) || ...
                                min(UqTempSum)+1 == numel(TempMet) || ...
                                min(UqTempSum)-1 == numel(TempMet)
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
% Identify transport reactions (diffusion)
[TransKeggT,TransBiggT] = CheckTransRxns(TransKegg,TransBiGG);

if ~isempty(TransBiggT{1})
    BiGGID = [BiGGID,TransBiggT];
    AllRxn = [AllRxn,TransKeggT];
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
% Check with restricted condition
[B1,K1] = Bigg2KeggRestricted(D,Metkegg,rxnData,cpds,cpd2rxn);
Bdiff = BiGGID(~ismember(BiGGID,B1));
Kdiff = AllRxn(~ismember(BiGGID,B1));
if any ([MultiBID{:}])
    Fileid = fopen(fullfile(Pth,'PrunedRxns.txt'),'r');
    Prunedata = textscan(Fileid,'%s %s %d %[^\n]');
    Prunedrxns = Prunedata{1};
    fclose(Fileid);
    MultiBIDtemp = MultiBID(~ismember(MultiBID,Prunedrxns));
    MultiRxntemp = MultiRxn(~ismember(MultiBID,Prunedrxns));
    Bmultidiff = MultiBIDtemp(~ismember(MultiBIDtemp,Bdiff));
    Kmultidiff = MultiRxntemp(~ismember(MultiBIDtemp,Bdiff));
    Tctr = 1; BmultidiffT = {[]}; KmultidiffT = {[]};
    for ct = 1:numel(Bmultidiff)
        for ci = 1:numel(Kmultidiff{ct})
            BmultidiffT{Tctr} = Bmultidiff{ct};
            KmultidiffT{Tctr} = Kmultidiff{ct}{ci};
            Tctr = Tctr + 1;
        end
    end
    if ~isempty(BmultidiffT{1})
        Bdiff = [Bdiff,BmultidiffT];
        Kdiff = [Kdiff,KmultidiffT];
    end
end
% Check with bigg rxn ref.
[bref,kref,Bdiff,Kdiff] = biggref(Bdiff,Kdiff);
% -----------------------------------
if ~isempty(Bdiff{1})
    Diffdata = cell(numel(Bdiff),3);
    Diffdata(:,1) = Bdiff; Diffdata(:,2) = Kdiff;
    for ct = 1:numel(Bdiff)
        Diffdata{ct,3} = 0;
    end
    setappdata(0,'Diffdata',Diffdata);
    waitfor(Bigg2KeggTable)
    Finalv=getappdata(0,'Finalvalidity');
    Fileid = fopen(fullfile(Pth,'PrunedRxns.txt'),'a');
    checkb = Finalv(:,1);
    checkv = Finalv(:,3);
    checkv = [checkv{:}]; checkb1 = ({});
    for ct = 1:size(Finalv,1)
        if any(strcmp(Finalv{ct,1},checkb1))
            continue
        end
        checkb1{ct} = Finalv{ct,1};
        checkloci = find(ismember(checkb,Finalv{ct,1}));
        checkfinal = checkv(checkloci);
        if ~any(checkfinal)
            for ct1 = 1:numel(checkfinal)
                fprintf(Fileid,'\n%s\t%s\t%d',Finalv{checkloci(ct1),1},...
                    Finalv{checkloci(ct1),2},...
                    0);
            end
        else
            for ct1 = 1:numel(checkfinal)
                if checkv(checkloci(ct1))
                    fprintf(Fileid,'\n%s\t%s\t%d',Finalv{checkloci(ct1),1},...
                        Finalv{checkloci(ct1),2},1);
                end
            end
        end
    end
    fclose(Fileid);
    AT=getappdata(0);
    AT1 = fieldnames(AT);
    for co = 1:length(AT1)
        if strcmp(AT1{co},'Diffdata') || strcmp(AT1{co},'Finalvalidity')
            rmappdata(0,AT1{co});
        end
    end
end
if ~isempty(bref{1})
    rmvtemp = ismember(BiGGID,bref);
    BiGGID(rmvtemp) = [];
    AllRxn(rmvtemp) = [];
    BiGGID = [BiGGID,bref];
    AllRxn = [AllRxn,kref];
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

% Extract multi-step reactions from KEGG. This step may need internet
% connection
[Steprxns,Insteprxns] = MultiRxns(AllRxn);
[BiGGID,AllRxn] = AddLumpRxns(Steprxns,Insteprxns,BiGGID,AllRxn);

% Save output data
B2Kegg.B=BiGGID; B2Kegg.K=AllRxn;
FileName = [ModelName1,'KEGG.mat'];
fname = fullfile(Pth,'BiGG2KEGG',FileName);
save (fname,'B2Kegg')

if any ([MultiBID{:}])
    fprintf('There exists %.0f BiGG rxns with multiple KEGG IDs:\n',...
        numel(MultiBID));
    for ct = 1:numel(MultiBID)
        fprintf('%s:\t',MultiBID{ct});
        fprintf('%s  ',MultiRxn{ct}{:});
        fprintf('\n')
    end
else
    disp('There exists no BiGG rxns with multiple KEGG IDs for this model')
end

%==========================================================================
%%===========================Subfunctions==================================
%==========================================================================
function Cbs = MetCombs(Metkegg,Prd)
% OutMet = A 2D cell with rows and columns corresponding to combinations
% and compounds, respectively.
% OutMet generates all combinations among Metkegg components. This is done 
% to search through KEGG REACTIONs for all possibilities.
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
% Check if H+ is peresent in Cbs. As KEGG reactions are not charge
% balanced, several BiGG reactions will not be found. To resolve this
% issue, new rows withouth H+ will be appended to Cbs.
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
    % On the contrary to above-mentioned case, the opposite case is also 
    % observed. In such a case, H+ will be added as an extra metabolite to
    % each row in CbsT, and final Cbs will contain previous set of
    % metabolites plus H+.
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
% This subfunction uses PrunedRxns.txt, a set of manually found equivalence
% KEGG rxn IDs, to prune current rxn mapping. This set of reactions are
% generated based on univeral reactions in BiGG databse, and can
% be utilized for every desired model.
Fileid = fopen('PrunedRxns.txt','r');
PrAll = textscan(Fileid,'%s %s %d %[^\n]');
PrB = PrAll{3};
PrBY = find(PrB); % To be added
PrBN = ~PrB; % To be removed
TBrxns = PrAll{1}; TKrxns = PrAll{2};
BrxnY = TBrxns(PrBY); KrxnY = TKrxns(PrBY);
BrxnN = TBrxns(PrBN);
% Which BiGGID elements are in BrxnY?
% [Fx1,Fx2] = ismember(BrxnY,BiGGID);
% Fx2(Fx2==0) = [];
% AllRxn(Fx2) = KrxnY(Fx1);
BiGGInBrxnY = BiGGID(ismember(BiGGID,BrxnY)); 
KeggOfBrxnY = KrxnY(ismember(BrxnY,BiGGID));
BrxnYTemp = BrxnY(ismember(BrxnY,BiGGID));
SharedRxns = cell(1,numel(BiGGInBrxnY));
Tcount = 1; TempSharedRxns = ({}); TempBiGGID = ({});
for count = 1:numel(BrxnYTemp)
    T= find(ismember(BiGGInBrxnY,BrxnYTemp{count}));
    T1 = find(ismember(BrxnYTemp,BiGGInBrxnY(T(1))));

    for count2 = 1:numel(T)
        % Extract only the first element:: T1(1); others are stored in
        % TempSharedRxns and will be added to AllRxn
        SharedRxns(T(count2)) = KeggOfBrxnY(T1(1));
        if numel(T1)>1 && ~any(ismember(TempBiGGID,BiGGInBrxnY(T(1))))
            for count3 = 2:numel(T1)
                TempSharedRxns{Tcount} = KeggOfBrxnY{T1(count3)};
                TempBiGGID{Tcount} = BrxnYTemp{T1(count3)};
                Tcount = Tcount + 1;
            end
        end
    end

end
AllRxn(ismember(BiGGID,BrxnY)) = SharedRxns;
if ~isempty(TempSharedRxns)
    Addlen = numel(AllRxn)+1:numel(AllRxn)+numel(TempSharedRxns);
    AllRxn(Addlen) = TempSharedRxns;
    BiGGID(Addlen) = TempBiGGID;
end
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
ExcRxns = {'R03082','R01651','R01210','R08767','R01651','R03171'};
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
%--------------------------------------------------------------------------
function [TransKeggT,TransBiggT] = CheckTransRxns(TransKegg,TransBiGG)
% Read cpd2weight ===================
Fileid = fopen('cpd2weight.txt','r');
cpd2weight = textscan(Fileid,'%s %f');
cpd4weight = strrep(cpd2weight{1},'cpd:','');
Weight = cpd2weight{2};
fclose(Fileid);
% Read rxns with 2, 4 and 6 components ========
Fileid = fopen('W4TranRxns2.txt','r');
rxn2t = textscan(Fileid,'%s %f %f');
rxn2 = rxn2t{1};
W21T = rxn2t{2};
W22T = rxn2t{3};
W2 = zeros(numel(W22T),2);
for ct = 1:numel(W21T)
    W2(ct,1) = W21T(ct);
    W2(ct,2) = W22T(ct);
end
fclose(Fileid);
Fileid = fopen('W4TranRxns4.txt','r');
rxn4t = textscan(Fileid,'%s %f %f %f %f');
rxn4 = rxn4t{1};
W41T = rxn4t{2};
W42T = rxn4t{3};
W43T = rxn4t{4};
W44T = rxn4t{5};
W4 = zeros(numel(W44T),2);
for ct = 1:numel(W44T)
    W4(ct,1) = W41T(ct);
    W4(ct,2) = W42T(ct);
    W4(ct,3) = W43T(ct);
    W4(ct,4) = W44T(ct);
end
fclose(Fileid);
CurLen = 1;
TransKeggT = {[]}; TransBiggT = {[]};
for ct = 1:numel(TransKegg)
    Text1 = ['Checking equivalent transport reaction for ',TransBiGG{ct},'\n'];
    fprintf(Text1)
    RxnComNo = numel(TransKegg{ct});
    RxnWeight = Weight(ismember(cpd4weight,TransKegg{ct}));
    if RxnComNo == 2
        rxn24 = rxn2;
        AllW = W2;
    elseif RxnComNo == 4
        rxn24 = rxn4;
        AllW = W4;
    end
    if isempty(RxnWeight)
        continue
    end
    Reslt = bsxfun(@eq,AllW,RxnWeight);
    [mID,~]=find(Reslt);
    if ~isempty(mID) && (RxnComNo.*numel(unique(mID)) == numel(mID))
        CandRxns1 = rxn24(unique(mID));
        for ct1 = 1:numel(CandRxns1)
            TransBiggT{CurLen} = TransBiGG{ct};
            TransKeggT{CurLen} = CandRxns1{ct1};
            CurLen = CurLen + 1;
        end
    end
end
%--------------------------------------------------------------------------
function[boutref,koutref,bout,kout] = biggref(bin,kin)
% Read BiGG data
bout = {[]}; kout = {[]};  boutref = {[]}; koutref = {[]};
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Fileid = fopen(fullfile(Pth,'bigg_models_reactions.txt'),'r');
biggrxns = textscan(Fileid,'%s %[^\n]');
fclose(Fileid);
biggid = biggrxns{1}; biggdes = biggrxns{2};
clear biggrxns
ct1 = 1; ct11 = 1;
while numel(bin)
    bintemp = bin{1};
    kintemp = kin(ismember(bin,bin{1}));
    bin1loc = find(ismember(bin,bin{1}));
    destemp = biggdes(ismember(biggid,bintemp));
    temploci = regexp(destemp,...
        'http://identifiers.org/kegg.reaction/R\S{5}','match');
    if isempty(temploci) || isempty(temploci{1})        
        for ct3 = 1:numel(bin1loc)
            bout{ct1} = bintemp;
            kout{ct1} = kintemp{ct3};
            ct1 = ct1 + 1;
        end
    else
        temploci = temploci{1};
        for ct2 = 1:numel(temploci)
            regtemp = regexp(temploci{ct2},'R\S{5}','match');
            boutref{ct11} = bintemp;
            koutref{ct11} = regtemp{1};
            ct11 = ct11 + 1;
        end
    end
    bin(bin1loc) = [];
    kin(bin1loc) = [];
end