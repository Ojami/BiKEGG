function NetDraw (Outflx,RxnCds,MapChoice,InOpts)
% NetDraw
% uses KGML file of global metabolic pathway of KEGG (map01100)
% to extract KEGG rxn IDs and their corresponding coordinates. 
% The overall structure of NetDraw is similar to that of KeggDraw
% with minor modifications on subfunctions. Further details can be found in
% the BiKEGG user manual
% 
% Inputs:
% Outflx = A mXn matrix containing input flux data, which m corresponds to
%        KEGG reaction IDs and n corresponds to time-series values.
% RxnCds = A cell array of strings comprising of KEGG reaction IDs (m rows).
% MapChoice = KEGG map ID.
% 
% Optional inputs:
% InOpts = A struct containing the following fields:
% SVal = Check if user wants to save images : 1 (yes) or 0 (no), default:0.
% Timevalue = The time step between two consecutive images (default: 0.1 s).
% Netstat = Cehck if user wants to use KEGG maps and KGML files online(1)
%           /offline(0) (default = 1). Data in KEGGmpas folder will be
%           used in case of offline mode.
% Svid = Check if user wants to convert a sequence of images to an animated
%		 GIF file: 2 (GIF), 0 (no); default: 0 :: Video option is not valid
%		 for NetDraw
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% June. 2016

% Optional inputs ---------------------------------------------------------
if exist('InOpts','var')
    OptFields = fieldnames(InOpts);
    for count = 1:numel(OptFields)
        if strcmpi('SVal',OptFields{count})
            SVal = InOpts.(OptFields{count});
        elseif strcmpi('Timevalue',OptFields{count})
            Timevalue = InOpts.(OptFields{count});
        elseif strcmpi('Netstat',OptFields{count})
            Netstat = InOpts.(OptFields{count});
        elseif strcmpi('Svid',OptFields{count})
            Svid = InOpts.(OptFields{count});
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
% Read basemap data ==================
basemap = BaseMapReader(Netstat);
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
fclose(Fileid1);
AllRN = ({});
AllRNids = [];
AllRNcol = ({});
Allflx = [];
AllrxnCds = ({}); 
AllBiGGI = ({});
for mpct = 1:numel(MapChoice)  
    rxns4thismap = Rxns1(ismember(Maps,MapChoice{mpct})); % rxns in the current map
    wherxn = find(ismember(basemap.r.rxn,rxns4thismap));
    RN = basemap.r.rxn(wherxn); % rxns of current map in basemap 
    RNids = basemap.r.rxnid(wherxn);
    RNcol = basemap.r.col(wherxn);
    %==========================================================================
    % Analyzing flux data ====================================================
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
    [MultiFlx,MultiRxn,MultiBigg] = ChckMultiRxns(BiGG4KeggDraw,RxnCds,rxns4thismap);
    MultiRxn1 = []; MultiFlx1 = []; MultiRN = []; MultiRNids = []; MultiRNcol = [];
    MultiBigg1 = [];
    for m1 = 1:numel(MultiRxn)
        where_multi = find(ismember(basemap.r.rxn,MultiRxn{m1}));
        if ~isempty(where_multi)
            MultiRxn1 = [MultiRxn1;MultiRxn(m1)];
            MultiFlx1 = [MultiFlx1;MultiFlx(m1,:)];
            MultiRN = [MultiRN,basemap.r.rxn(where_multi)];
            MultiRNids = [MultiRNids,basemap.r.rxnid(where_multi)];
            MultiRNcol = [MultiRNcol,basemap.r.col(where_multi)];
            MultiBigg1 = [MultiBigg1;MultiBigg(m1)];
        end
    end
    MultiRxn = MultiRxn1;
    MultiFlx = MultiFlx1;
    MultiBigg = MultiBigg1;
    if size(RN,2) > size(RN,1)
        RN = RN';
    end
    if size(MultiRN,2) > size(MultiRN,1)
        MultiRN = MultiRN';
    end
    if size(RNids,2) > size(RNids,1)
        RNids = RNids';
    end
    if size(MultiRNids,2) > size(MultiRNids,1)
        MultiRNids = MultiRNids';
    end
    if size(RNcol,2) > size(RNcol,1)
        RNcol = RNcol';
    end
    if size(MultiRNcol,2) > size(MultiRNcol,1)
        MultiRNcol = MultiRNcol';
    end
    if size(MultiRxn,2) > size(MultiRxn,1)
        MultiRxn = MultiRxn';
    end
    if size(MultiBigg,2) > size(MultiBigg,1)
        MultiBigg = MultiBigg';
    end
    AllrxnCds = [AllrxnCds;MultiRxn];
    AllRN = [AllRN;RN;MultiRN];
    AllRNids = [AllRNids;RNids;MultiRNids];
    Allflx = [Allflx;MultiFlx];
    AllBiGGI = [AllBiGGI;MultiBigg];
    AllRNcol = [AllRNcol;RNcol;MultiRNcol];
    
end
[rxnCds,flx1] = ConsistRxns(BiGG4KeggDraw,RxnCds,Outflx,MapChoice);
if size(rxnCds,2) > size(rxnCds,1)   
    rxnCds = rxnCds';
end
[f1,~]=ismember(rxnCds,AllRN);
flx=flx1(f1,:);
rxnCds = rxnCds(f1);
BiGG4KeggDraw1 = getappdata(0,'BiGG4KeggDraw1');
BiGGI = BiGG4KeggDraw1(f1);
if size(BiGGI,2) > size(BiGGI,1)
    BiGGI = BiGGI';
end
AllrxnCds = [AllrxnCds;rxnCds];
Allflx = [Allflx;flx];
AllBiGGI = [AllBiGGI;BiGGI];

MixBiggRxn = strcat(AllrxnCds,AllBiGGI);
[~, nflx] = unique(MixBiggRxn);
Allflx = Allflx(nflx,:);
AllrxnCds = AllrxnCds(nflx); AllBiGGI = AllBiGGI(nflx);
% Check repeated rxns on this map and let user choose one,among several
% BiGG IDs, for which the flux value will be displayed on the map.
[flx,rxnCds] = ChckRptRxns(Allflx,AllrxnCds,MapChoice,AllBiGGI);

[rxnCds,nflx1] = unique(rxnCds);
flx = flx(nflx1,:);
InrxnCds = rxnCds; Inflx = flx;
%----------------------------------------------------------------------
if isempty(flx)
    uiwait(msgbox('The flux data are not is current pathway!',...
        'Warning','warn'));
    close (h);
    return
end
%======================================================================
load(which('Pixdata.mat'))
fnames = fieldnames(pix_x);
copt2 = 1;
rxnCds1={0};rxnPlusid = ({}); rxnCol = ({}); these_ids1 = (0);
for copt = 1:numel(AllRN)
    Temprxn = find(ismember(AllRN,AllRN{copt}));
    for copt1 = 1:numel(Temprxn)
        these_ids1(copt2) = AllRNids(Temprxn(copt1));
        rxnPlusid{copt2} = [AllRN{Temprxn(copt1)},num2str(these_ids1(copt2))];
        rxnCol{copt2} = AllRNcol{Temprxn(copt1)}; % Identify rxns based on their color to remove rxns in other pathways
        rxnCds1{copt2} = AllRN{copt};
        copt2 = copt2 + 1;
    end
end
[rxnPlusid,plusn] = unique(rxnPlusid);
rxnCds = rxnCds1(plusn);
rxnCol = rxnCol(plusn);
these_ids1 = these_ids1(plusn);
rcid = basemap.rc.rxnid;
rc_subid = basemap.rc.subid;
rc_proid = basemap.rc.prodid;
[shared_ids,~,rcid_n] = intersect(these_ids1,rcid);
S_rxn = ({}); s3 = 1; S_mat = (0); S_rxnp = ({}); S_cpdcum = []; S_color = ({});
for s11 = 1:numel(shared_ids)
    S_sub_temp = rc_subid{rcid_n(s11)};
    S_pro_temp = rc_proid{rcid_n(s11)};
    S_cpdcum = [S_cpdcum,S_sub_temp,S_pro_temp];
end
S_cpdcum = unique(S_cpdcum); % Unique all compounds for these rxns

% Generate unsigned-stoichiometric matrix and subsequent rxn-rxn graph
for s1 = 1:numel(shared_ids)
    where_rxn = find(these_ids1==shared_ids(s1));
    S_sub_temp = rc_subid{rcid_n(s1)};
    S_pro_temp = rc_proid{rcid_n(s1)};
    for s2 = 1:numel(where_rxn)
        S_rxnp{s3} = rxnPlusid{where_rxn(s2)};
        S_rxn{s3} = rxnCds{where_rxn(s2)}; % Rxns for 
        S_color{s3} = rxnCol{where_rxn(s2)};
        [~,S_sub_temp_idx,~] = intersect(S_cpdcum,S_sub_temp);
        [~,S_pro_temp_idx,~] = intersect(S_cpdcum,S_pro_temp);
        if size(S_sub_temp_idx,1) > size(S_sub_temp_idx,2)
            S_sub_temp_idx = S_sub_temp_idx';
        end
        if size(S_pro_temp_idx,1) > size(S_pro_temp_idx,2)
            S_pro_temp_idx = S_pro_temp_idx';
        end
        S_mat([S_sub_temp_idx,S_pro_temp_idx],s3) = 1; 
        s3 = s3 + 1;
    end
end
S_mat = sparse(S_mat);
% Cpd-cpd graph------------------------------------------------------------
% S_graph = zeros(size(S_mat,1)); t_graph = ({});
% for i1=1:size(S_mat,1)
%     nr = find(S_mat(i1,:));
%     if ~isempty(nr)
%         for i2 = 1:numel(nr)
%             t_graph{i2} = find(S_mat(:,nr(i2)));
%             if size(t_graph{i2},1) > size(t_graph{i2},2)
%                 t_graph{i2} = t_graph{i2}';
%             end
%         end
%         t_con=[t_graph{:}]; t_con=unique(t_con);
%         clear t_graph
%         S_graph(i1,setdiff(t_con,i1)) = 1;
%     end
% end
% Rxn-rxn graph -----------------------------------------------------------
S_graph = zeros(size(S_mat,2)); t_graph = ({});
for i1=1:size(S_mat,2)
    nr = find(S_mat(:,i1));
    if ~isempty(nr)
        for i2 = 1:numel(nr)
            t_graph{i2} = find(S_mat(nr(i2),:));
            if size(t_graph{i2},1) > size(t_graph{i2},2)
                t_graph{i2} = t_graph{i2}';
            end
        end
        t_con=[t_graph{:}]; t_con=unique(t_con);
        clear t_graph
        S_graph(i1,setdiff(t_con,i1)) = 1;
    end
end
b_graph = biograph(S_graph);
[S_bins, where_bins] = conncomp(b_graph); % Connected parts of rxn-rxn graph
freq_bins = histc(where_bins,1:S_bins);
[~,max_freq_bin] = max(freq_bins);
Secplace_freq_bin = find(freq_bins>2);
 % All cpds present in the strongest (most connected) part of cpd-cpd
 % graph:
% core_cpds = S_cpdcum(where_bins == max_freq_bin);
core_rxns= S_rxnp(where_bins == max_freq_bin); % All cpds present in the strongest (most connected) part of rxn-rxn graph
Secids = [];
for sctr = 1:numel(Secplace_freq_bin)
    Secids = [Secids,find(where_bins == Secplace_freq_bin(sctr))];
end
Secondcore_rxns= S_rxnp(Secids);
[~,orph_id] = intersect(where_bins,find(freq_bins == 1));
% [~,double_id] = intersect(where_bins,find(freq_bins == 2));
% double_rxns = S_rxnp(double_id);
orphan_rxns = S_rxnp(orph_id); % Single rxns which are not connected to any part: These rxns should be removed.
% Identify repeated rxns on the final map and remove redundancies based on
% core_cpds (from graph connectivity assessment): Among
% different version of same rxn (same rxn name and different rxn id), only
% those with high connectedness are retained.
% Dump color check ////////////////////////////////////////////////////////
Prominent_col = S_color(where_bins == max_freq_bin);
Prominent_col = unique(Prominent_col);
%//////////////////////////////////////////////////////////////////////////
[~,S_rxn_idx,~] = unique(S_rxn);
S_rxn_rptd = S_rxn(setdiff(1:numel(S_rxn),S_rxn_idx)); 
rmv_rxnPlus_id = []; plus_id = 1;
for c1 = 1:numel(S_rxn_rptd)
    where_rptd = find(ismember(S_rxn,S_rxn_rptd{c1}));
    rmv_mat = []; rmv_id = 1;
    for c2 = 1:numel(where_rptd)
%         cpds4rptd = find(S_mat(:,where_rptd(c2))); % Only for cpd-cpd graph
        if ~isempty(setdiff(S_rxnp(where_rptd(c2)),core_rxns)) && ...
              ~isempty(setdiff(S_color(where_rptd(c2)),Prominent_col))     
            rmv_mat(rmv_id) = where_rptd(c2); % Those who failed the connectivity test
            rmv_id = rmv_id + 1;
        end
    end
    for c3 = 1:numel(rmv_mat)
        rmv_rxnPlus_id(plus_id) = rmv_mat(c3);
        plus_id = plus_id + 1;
    end
end
% Final check for orphan rxns after removal of rmv_rxnPlus_id ids
S_graph1 = S_graph;
S_graph1(rmv_rxnPlus_id,:) = [];
S_graph1(:,rmv_rxnPlus_id) = [];
b_graph1 = biograph(S_graph1);
[S_bins1, where_bins1] = conncomp(b_graph1); 
freq_bins1 = histc(where_bins1,1:S_bins1);
[~,orph_id1] = intersect(where_bins1,find(freq_bins1 == 1));
S_rxnp1 = S_rxnp;
S_rxnp1(rmv_rxnPlus_id) = [];
orphan_rxns1 = S_rxnp1(orph_id1);
%--------------------------------------------------------------------------
% Check for orphan multiple rxns or their children : These rxns do not have
% any cpd on basemap, so, do not appear in rxn-rxn graph and therefore, 
% cannot be identified. To identify these rxns, cpd2rxn is used and through
% comparison with cpds in core rxns, its orphanage state can be stablished.
multirxns_id = setdiff(these_ids1,shared_ids);
[~, multirxns_id] = intersect(these_ids1,multirxns_id);
multirxns = rxnCds(multirxns_id);
% Read cpd2rxn 
Fileid = fopen('cpd2rxn.txt','r');
cpd2rxn = textscan(Fileid,'%s %s');
rxnData = strrep(cpd2rxn{2},'rn:','');
cpds = strrep(cpd2rxn{1},'cpd:','');
fclose(Fileid);
% Identify cpds of rxns with good connectedness (>2)
SecCore_ids = these_ids1(ismember(S_rxnp,Secondcore_rxns));
[~,SecCore_idx] = intersect(basemap.rc.rxnid,SecCore_ids);
core2_subs = basemap.rc.sub(SecCore_idx); core2_subs = [core2_subs{:}];
core2_pros = basemap.rc.prod(SecCore_idx); core2_pros = [core2_pros{:}];
core2_cpds = [core2_subs,core2_pros]; core2_cpds = unique(core2_cpds);
m2 = 1; rmv_mult_id = [];
for m1 = 1:numel(multirxns)
    curr_cpds = cpds(ismember(rxnData,multirxns(m1)));
    if isempty(intersect(curr_cpds,core2_cpds)) || ...
            (~any(strcmp(MapChoice,'00640')) && strcmp(multirxns(m1),'R00928'))
        rmv_mult_id(m2) = multirxns_id(m1);
        m2 = m2 + 1;
    end
end
rxnCds(rmv_mult_id) = [];
rxnPlusid(rmv_mult_id) = [];
these_ids1(rmv_mult_id) = [];
%--------------------------------------------------------------------------
% Modify final rxns and remove orphans and redundants
rmv_rxnPlus = S_rxnp(rmv_rxnPlus_id); % Rxns to be removed from final map
rxnCds(ismember(rxnPlusid,rmv_rxnPlus)) = [];
these_ids1(ismember(rxnPlusid,rmv_rxnPlus)) = [];
rxnPlusid(ismember(rxnPlusid,rmv_rxnPlus)) = [];
rxnCds(ismember(rxnPlusid,orphan_rxns)) = [];
these_ids1(ismember(rxnPlusid,orphan_rxns)) = [];
rxnPlusid(ismember(rxnPlusid,orphan_rxns)) = []; % Remove orphan rxns
rxnCds(ismember(rxnPlusid,orphan_rxns1)) = [];
these_ids1(ismember(rxnPlusid,orphan_rxns1)) = [];
rxnPlusid(ismember(rxnPlusid,orphan_rxns1)) = []; % Remove orphan rxns after removal of redundant rxns
%Flux data arrangement-----------------------------------------------------
[On2,~] = ismember(rxnPlusid,fnames); % Idx of all rxns in these maps present in basemap
field_intsect = rxnPlusid(On2); % Reference to all rxns in these maps for reading coordinates
rxnCds = rxnCds(On2); % Rearrange all rxn IDs in these maps

overlay_rxnPlus = rxnPlusid(ismember(rxnCds,InrxnCds)); % Idx of rxns for which flux rates are provided 
[~,o2] = ismember(rxnCds,InrxnCds); % Rearrange rxn IDs and fluxes based on Idx
F_flx = [];
o2 = o2(find(o2)); OverlayRxnsIDs = InrxnCds(o2);
for m1 = 1:numel(o2)
    F_flx(m1,:) = Inflx(o2(m1),:);
end
[On1,~] = (ismember(overlay_rxnPlus,fnames)); % Find idxs of rxns for which flux rates are provided and present in basemap
field_intsect_overlay = overlay_rxnPlus(On1); % Reference to all rxns with flux rates
F_flx = F_flx(On1,:); % All flux rates
OverlayRxnsIDs = OverlayRxnsIDs(On1); % KEGG rxn IDs corresponding to field_intsect_overlay
flx = F_flx;
%--------------------------------------------------------------------------
LineStrength = 7;
pixx = []; pixy =[];
for i4 = 1:numel(field_intsect)
    pixx = [pixx;pix_x.(field_intsect{i4})];
    pixy = [pixy;pix_y.(field_intsect{i4})];
end
min_x = min(pixx); max_x = max(pixx); min_y = min(pixy); max_y = max(pixy);
pixx = []; pixy = []; Allx = ({}); Ally = ({}); Allx4Arrow = ({});Ally4Arrow = ({});
RefOverlapData = cell(numel(field_intsect),1);
I = ones(max_y-min_y+50,max_x-min_x+50,3);
I(:,:,1) = I(:,:,1).*0.8; I(:,:,2) = I(:,:,2).*0.992; I(:,:,3) = I(:,:,3).*0.8;
 
for i4 = 1:numel(field_intsect)
    pixx = pix_x.(field_intsect{i4})-min_x+10;
    pixy = pix_y.(field_intsect{i4})-min_y+10;
    RefOverlapData{i4} = [num2str(numel(pixx)),',',...
        num2str(mean(pixx)),',',num2str(numel(pixy)),',',num2str(mean(pixy))];
    Allx4Arrow{i4} = pixx; Ally4Arrow{i4} = pixy;
    [I,pixx_temp,pixy_temp,pixx_temp1,pixy_temp1] = pixFill(I,pixx,pixy,LineStrength,[0 0 0]);
    if numel(pixx_temp)>20; pixx_temp(1:10)=[]; pixx_temp(end-10:end)=[];end;
    if numel(pixy_temp)>20; pixy_temp(1:10)=[]; pixy_temp(end-10:end)=[];end;
    if numel(pixx_temp1)>20; pixx_temp1(1:10)=[]; pixx_temp1(end-10:end)=[];end;
    if numel(pixy_temp1)>20; pixy_temp1(1:10)=[]; pixy_temp1(end-10:end)=[];end;
    Allx{i4} = [pixx_temp;pixx_temp1];
    Ally{i4} = [pixy_temp;pixy_temp1];
    clear pixx pixy
end

% Compounds ---------------------------------------------------------------
[~,n1] = intersect(basemap.rc.rxnid,these_ids1);
a1 = basemap.rc.subid(n1); a1 = [a1{:}];
b1 = basemap.rc.prodid(n1); b1 = [b1{:}];
a=[a1,b1];
a=unique(a);
[~,n21]=intersect(basemap.c.cpdid,a);
[~,n22]=intersect(basemap.c.glid,a);
cpdnames1 = basemap.c.cpd(n21);
cpdnames2 = basemap.c.gl(n22);
cpdnames = [cpdnames1,cpdnames2];
cpdxy1 = basemap.c.cpdxy(n21,:);
cpdxy2 = basemap.c.glxy(n22,:);
cpdxy = [cpdxy1;cpdxy2];
x_cpd = cpdxy(:,1) - min_x - 18;
y_cpd = cpdxy(:,2) - min_y;

% Flux overlay ============================================================
% The time step between two consecutive images:
if ~exist('Timevalue','var') || isempty(Timevalue)
    Timevalue=0.1;
end

% Time points: % Disabled in current settings /////////////////////////////
% if ~exist('mTime','var') || isempty(mTime)
%     mTime=1:size(Outflx,2);
% end
% 
% if ~exist('TimeUnit','var') || isempty(TimeUnit)
%     TimeUnit=' ';
% end
%//////////////////////////////////////////////////////////////////////////
if ~exist('Svid','var') || isempty(Svid)
    Svid = 0; %Does not convert images to video
end
% Check wether user wants to save images (using 'Save images checkbox')
if ~exist('SVal','var') || isempty(SVal)
    SVal = 0; % Does not save images.
end
if SVal
    FolderNm='NetDraw output images';
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
Cm=jet(numel(field_intsect_overlay)); % [Warning] Changing Cm matrix is not recommended
if numel(MapChoice) > 4
    TempTitle = MapChoice(1:4);
    TempTitle = strjoin(TempTitle,',');
    TempTitle = [TempTitle,',...'];
else
    TempTitle = strjoin(MapChoice,',');
end
%Draw images on this panel
H1 = figure;
h1 = axes;
set(H1,'position',get(0,'screensize'));
set(h1,'position',[0,0,1,1],'xtick',[],'ytick',[]);
set(H1,'toolbar','figure');
set(H1,'menubar','figure');
set(H1,'name',TempTitle,'numbertitle','off')
% Arbitrary color and font properties % Disabled in current settings///////
% if ~exist('Viso','var') || ~Viso
%     TextFont.FontName = 'Arial';
%     TextFont.FontWeight = 'normal';
%     TextFont.FontAngle = 'normal';
%     TextFont.FontUnits = 'points';
%     TextFont.FontSize = 8;
% end
% if exist('Viso','var') && Viso
%     waitfor(VisualProp)
%     TextFont = getappdata(0,'TextFont');
%     if isempty(TextFont)
%         TextFont.FontName = 'Arial';
%         TextFont.FontWeight = 'normal';
%         TextFont.FontAngle = 'normal';
%         TextFont.FontUnits = 'points';
%         TextFont.FontSize = 8;
%     end
% end
%//////////////////////////////////////////////////////////////////////////
% Visualization -----------------------------------------------------------
flx_abs = abs(flx);
for k=1:size(flx,2)
    ArX = ({}); ArY = ArX; McmCollect = zeros(size(flx,1),3); 
    OverlapData = cell(size(flx,1),1);
    for m=1:size(flx,1)
        % Normalize colors and insertion into enzyme boxes
        flxMax=max(abs(flx(:)));
        flxMin=min(abs(flx(:)));
        if flxMax==flxMin % One row (i.e. for 1 rxn)
            flxMax=max(abs(flx(:)));
            flxMin=min(abs(flx(:)));
        end
        flxSl=(size(Cm,1)-1)./(flxMax-flxMin);
        if isinf(flxSl) % in case of flxMax == flxMin
            flxSl = 1e15.*(size(Cm,1)-1);
        end
        Mcm=ceil((flx_abs(m,k)-flxMax).*flxSl+size(Cm,1));

        % Show the image:
        pixx = pix_x.(field_intsect_overlay{m})-min_x+10;
        pixy = pix_y.(field_intsect_overlay{m})-min_y+10;
        OverlapData{m} = [num2str(numel(pixx)),',',...
            num2str(mean(pixx)),',',num2str(numel(pixy)),',',num2str(mean(pixy))];
        I = pixFill(I,pixx,pixy,LineStrength,[Cm(Mcm,1) Cm(Mcm,2) Cm(Mcm,3)]);
        ArX{m} = pixx; ArY{m} = pixy;
        McmCollect(m,:) = [Cm(Mcm,1) Cm(Mcm,2) Cm(Mcm,3)];               
        clear pixx pixy        
    end
     %---------------------------------------------------------------------
    [I,ArXo,ArYo,ArrowChek,RxnId] = OverlapChecker(I,ArX,ArY,OverlapData,...
        flx(:,k),McmCollect,field_intsect_overlay,LineStrength);
    [ArrowX,ArrowY,Cntrs] = ArrowheadProcess(Allx4Arrow,Ally4Arrow,...
        min_x,min_y,field_intsect,basemap,RxnId,RefOverlapData,'bigg',field_intsect_overlay);
    I = MapTrimmer(I,ArrowX,ArrowY,Cntrs,Allx4Arrow,Ally4Arrow,LineStrength,[0.8 0.992 0.8]);
    [ArrowX1,ArrowY1,Cntrs1,Director] = ArrowheadProcess(ArXo,...
        ArYo,min_x,min_y,field_intsect_overlay,basemap,[],[],'bigg',[]);
    I = MapTrimmerF(I,ArrowX1,ArrowY1,Cntrs1,ArXo,ArYo,RxnId,...
        field_intsect_overlay,LineStrength,[0.8 0.992 0.8]);
    % Compounds overlay ---------------------------------------------------
    for i1 = 1:size(cpdxy,1)
        cc1 = x_cpd(i1)-0.7*cpdxy(i1,4):x_cpd(i1)+cpdxy(i1,4)*0.7;
        cc2 = y_cpd(i1)-cpdxy(i1,4)*0.7:y_cpd(i1)+cpdxy(i1,4)*0.7;
        [x,y] = meshgrid(cc1,cc2);
        f1= ((x-x_cpd(i1)).^2+(y-y_cpd(i1)).^2) <=(cpdxy(i1,4)*0.7)^2;
        f2= ((x-x_cpd(i1)).^2+(y-y_cpd(i1)).^2) <=(cpdxy(i1,4)*0.55)^2;
        x1 = x(f1); y1= y(f1); x1=round(x1); y1=round(y1);
        x2 = x(f2); y2= y(f2); x2=round(x2); y2=round(y2);
        x1(~y1) = []; y1(~y1) = [];
        x2(~y2) = []; y2(~y2) = [];
        for i2 = 1:numel(x1) 
            I(y1(i2),x1(i2),1) = 0.502;
            I(y1(i2),x1(i2),2) = 0.502;
            I(y1(i2),x1(i2),3) = 0;
        end
        for i2 = 1:numel(x2)
            I(y2(i2),x2(i2),1) = 0;
            I(y2(i2),x2(i2),2) = 1;
            I(y2(i2),x2(i2),3) = 0.498;
        end
        clear f1 x1 x2 y1 y2 x y cc1 cc2
    end
    imshow(I,'InitialMagnification', 'fit','Parent',h1); 
    % Disabled in current settings ////////////////////////////////////////
    % Show time-series values on the corresponding image
%     TimePropUDF ('One',TextFont,TimeUnit,mTime,k,h1) 
    %//////////////////////////////////////////////////////////////////////
    arrowhead(ArrowX,ArrowY,[],0,4,[]); 
    flxDirect.D = Director; flxDirect.F = flx(:,k);
    arrowhead(ArrowX1,ArrowY1,ArrowChek,McmCollect,4,flxDirect);

    % Process pixhover settings--------------------------------------------
    Mdl = getappdata(0,'Mdl');
    Uni = load(which('UniModelKEGG.mat'));
    CRxns = getappdata(0,'CRxns'); % Map of consistent rxns
    if ~isempty(CRxns) && k == 1
        ConOverlayRxnID = OverlayRxnsIDs;
        for cnter = 1:numel(CRxns.C)
            ConOverlayRxnID(ismember(OverlayRxnsIDs,CRxns.C{cnter})) = ...
                {[CRxns.C{cnter},'(',CRxns.O{cnter},')']};
            OverlayRxnsIDs(ismember(OverlayRxnsIDs,CRxns.C{cnter})) = ...
                CRxns.O(cnter);
        end
    elseif isempty(CRxns)
        ConOverlayRxnID = OverlayRxnsIDs;
    end
    K = Uni.B2Kegg.K; K1 = Mdl.B2Kegg.K;
    B = Uni.B2Kegg.B; B1 = Mdl.B2Kegg.B;
    clear Uni
    HoverUs.K = K1; HoverUs.B = B1;
    HoverUs.R = OverlayRxnsIDs; HoverUs.P = field_intsect_overlay';
    HoverUs.CR = ConOverlayRxnID;
    HoverAll.R = rxnCds; HoverAll.P = field_intsect;
    HoverAll.K = K; HoverAll.B = B;
    pixHover(x_cpd,y_cpd,0.7*cpdxy(1,4),cpdnames,Allx,Ally,...
        HoverAll,HoverUs,abs(flx(:,k)),h1,H1)

    % Show colorbar--------------------------------------------------------
    flxstr=linspace(flxMin,flxMax,7);
    flxstr1=cell(numel(flxstr),1);
    for ir=1:numel(flxstr)
        flxstr1{ir,1}=num2str(flxstr(ir),2);
    end
    ch = colorbar('peer',h1);
    set(ch,'yticklabel',flxstr1)
    % Check for further settings ------------------------------------------
    IsWanting = input('Do you want to post-process the image? (y/n):','s');
    if isempty(IsWanting) || strcmp(IsWanting,'y')
        setappdata(0,'Img',I);
        setappdata(0,'ParentAx',h1);
        setappdata(0,'ParentFig',H1);
        PostData.field_intsect = field_intsect;PostData.Allx=Allx;
        PostData.Ally=Ally;PostData.min_x=min_x;PostData.min_y=min_y;
        PostData.RefOverlapData=RefOverlapData;PostData.Allx4Arrow=Allx4Arrow;
        PostData.Ally4Arrow=Ally4Arrow;PostData.ArX=ArX;PostData.ArY=ArY;
        PostData.OverlapData=OverlapData;PostData.field_intsect_overlay=field_intsect_overlay;
        PostData.ArXo=ArXo;PostData.ArYo=ArYo;PostData.ArrowChek=ArrowChek;
        PostData.RxnId=RxnId;PostData.ArrowX=ArrowX;PostData.ArrowY=ArrowY;
        PostData.Cntrs=Cntrs;PostData.ArrowX1=ArrowX1;PostData.ArrowY1=ArrowY1;
        PostData.Cntrs1=Cntrs1;PostData.McmCollect=McmCollect;
        PostData.cpdnames=cpdnames;PostData.cpdxy=cpdxy;PostData.x_cpd=x_cpd;
        PostData.y_cpd=y_cpd;PostData.flx=flx;PostData.rxnCds=rxnCds;
        PostData.max_x=max_x;PostData.max_y=max_y;PostData.basemap=basemap;
        PostData.OverlayRxnsIDs = OverlayRxnsIDs; PostData.k = k;
        PostData.ConOverlayRxnID = ConOverlayRxnID;
%         PostData.TimeUnit=TimeUnit;PostData.TextFont=TextFont;
%         PostData.mTime=mTime;
        setappdata(0,'PostData',PostData);
        uiwait(MapAdjuster);
    end
    % Save to image files within a new folder------------------------------
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
        SaveFlxImgs(getappdata(0,'ParentFig'),FolderPath,k,ImgRes,ImgFrmt);
    end
	%----------------------------------------------------------------------
    pause(Timevalue); % Time step between consecutive images
end
if Svid
    Img2Animation(Svid)
end

%==========================================================================
% Subfunctions ------------------------------------------------------------
%==========================================================================
function [flxOut,rxnCdsI] = ChckRptRxns (flx,rxnCds,CurrentMap,BiGGI) % NetDraw
% First, rxns are pruned based on manual assessment of pathways
rxnCdsI = rxnCds;

Fileid1 = fopen('rxnbiggmap.txt','r');
Rxnbmap = textscan(Fileid1,'%s %s %s %d');
Rtemp = Rxnbmap{1}; Btemp = Rxnbmap{2}; Mtemp = Rxnbmap{3};
Dtemp = Rxnbmap{4};
fclose(Fileid1);
AllrxnCdsI = ({}); AllBiGGI = ({}); Allflx = [];
for i1 = 1:numel(CurrentMap)
    if any(ismember(Mtemp,CurrentMap{i1}))
        ChR = Rtemp(ismember(Mtemp,CurrentMap{i1}));
        ChB = Btemp(ismember(Mtemp,CurrentMap{i1}));
        ChD = Dtemp(ismember(Mtemp,CurrentMap{i1}));
        for ct = 1:numel(ChR)
            n1_chr = find(ismember(rxnCdsI,ChR{ct}));
            n1_chb = find(ismember(BiGGI,ChB{ct}));
            if ~isempty(intersect(n1_chr,n1_chb))
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
    if size(rxnCdsI,2) > size(rxnCdsI,1)
        rxnCdsI = rxnCdsI';
    end
    if size(BiGGI,2) > size(BiGGI,2)
        BiGGI = BiGGI';
    end
    AllrxnCdsI = [AllrxnCdsI;rxnCdsI];
    AllBiGGI = [AllBiGGI;BiGGI];
    Allflx = [Allflx;flx];
end
MixBiggRxn = strcat(AllrxnCdsI,AllBiGGI);
[~, nflx] = unique(MixBiggRxn);
Allflx = Allflx(nflx,:);
AllrxnCdsI = AllrxnCdsI(nflx); AllBiGGI = AllBiGGI(nflx);
rxnCdsI = AllrxnCdsI;
BiGGI = AllBiGGI;
flx = Allflx;
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
    setappdata(0,'CurrentMap','All maps');
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
            && ~strcmp(AT1{co},'MultiK') && ~strcmp(AT1{co},'MultiFlx') ...
            && ~strcmp(AT1{co},'BiGGModel') && ~strcmp(AT1{co},'Mdl') ...
            && ~strcmp(AT1{co},'CRxns')
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
function [OutFlx,OutRxn,OutBigg] = ChckMultiRxns(BiGG4KeggDraw,RxnCds,rxns4thismap)
MultiB = getappdata(0,'MultiB'); 
MultiK = getappdata(0,'MultiK');
MultiFlx = getappdata(0,'MultiFlx');

OutFlx =[];
OutBigg = ({});
OutRxn = ({});
while numel(MultiB)
    Bigg1Loci = find(ismember(MultiB,MultiB{1}));
    KeggTemp = MultiK(Bigg1Loci);
    if all(ismember(KeggTemp,rxns4thismap)) % All child rxns are present in map
        FlxTemp = MultiFlx(Bigg1Loci(1),:);
        OutFlx = [OutFlx;FlxTemp];
        OutRxn = [OutRxn;RxnCds(ismember(BiGG4KeggDraw,MultiB{1}))];
        OutBigg = [OutBigg;MultiB{1}];
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
fclose(Fileid1);
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
fclose(Fileid);

fprintf('\nRxn consistency for %s(of %d): ','All pathways',numel(rxnIn))
ctr = 1; cnter = 1; OriginalRxns = ({}); ConsRxns = ({});
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
                W4TempCpds = Weight(ismember(cpd4weight,cpds4TempRxn));
                if numel(W4curCpds) == numel(W4TempCpds) &&...
                        isempty(setdiff(W4curCpds,W4TempCpds)) && ...
                        ~any(strcmp(TempRns{count1},{'R02630','R03232'})) % Exceptional condition 
                    rxnIn{CurLen} = TempRns{count1};
                    flxIn(CurLen,:) = TempFlx;
                    BiGG4KeggDraw{CurLen} = BiGG4KeggDraw{count};
                    CurLen = CurLen+1;
                    OriginalRxns{cnter} = rxnIn{count};
                    ConsRxns{cnter} = TempRns{count1};
                    cnter = cnter + 1;
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
if ~isempty(OriginalRxns) % A map between consisten rxns and those in model
    CRxns.O = OriginalRxns; CRxns.C = ConsRxns;
    setappdata(0,'CRxns',CRxns)
end