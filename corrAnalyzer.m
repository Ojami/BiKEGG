function corr = corrAnalyzer
% corrAnalyzer
% compares the reconciled reactions in BiGG and KEGG databases
% taken from MetaNetX repository and BiKEGG (specifically, Bigg2Kegg
% function). 
% WARNING: To execute this function, you need to include all BiGG models
% in a folder named "BiGG Models".

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Apr. 2016


clc
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth2 = fullfile(Pth,'BiGG2KEGG\*.mat');
Pth3 = fullfile(Pth,'BiGG2KEGG');
Pth22 = fullfile(Pth,'BiGG Models\*.xml');
Pth33 = fullfile(Pth,'BiGG Models');

% Read MetaNetX data
Fileid = fopen(fullfile(Pth,'MetaNetX.txt'),'r');
rawmetanet = textscan(Fileid,'%s %[^\n]');
fclose(Fileid);
metaid = rawmetanet{2}; metaref = rawmetanet{1};
metabigg = regexp(metaref,'bigg:\w*');
metabigg = ~cellfun('isempty',metabigg);
metabigg1 = metaid(metabigg);
for ct = 1:numel(metabigg1)
    metabigg11(ct) = cellstr(metabigg1{ct});
end
clear metabigg1
metabigg1 = metabigg11;
metabigg = metaref(metabigg); metabigg = strrep(metabigg,'bigg:','');
metakegg = regexp(metaref,'kegg:\w*');
metakegg = ~cellfun('isempty',metakegg);
metakegg1 = metaid(metakegg);
for ct = 1:numel(metakegg1)
    metakegg11(ct) = cellstr(metakegg1{ct});
end
clear metakegg1
metakegg1 = metakegg11;
metakegg = metaref(metakegg); metakegg = strrep(metakegg,'kegg:','');
clear metaref metaid rawmetanet metakegg11 metabigg11
netbigg = metabigg(ismember(metabigg1,metakegg1));
ct2 = 1; netkegg = ({});
for ct = 1:numel(metabigg1)
    netloci = find(ismember(metakegg1,metabigg1{ct}));
    if ~isempty(netloci)
        for ct1 = 1:numel(netloci)
            netkegg{ct2}{ct1} = metakegg{netloci(ct1)};
        end
        ct2 = ct2 + 1;
    end
end

% Read data from bigg2kegg
BMatfiles = dir(Pth2);
BMatfiles1 = dir(Pth22);
xmlname1 = struct2cell(BMatfiles1);
xmlname = xmlname1(1,:);
xmlname = strrep(xmlname,'.xml','');
ct1 = 1;
for ct = 1:numel(BMatfiles)
    if ~strcmp(BMatfiles(ct).name,'BiGG2KEGG_HMRbased-RECON1.mat') && ...
            ~strcmp(BMatfiles(ct).name,'Multirxns.mat') && ...
            ~strcmp(BMatfiles(ct).name,'UniModelKEGG.mat')
        BNames = BMatfiles(ct).name;
        btempname = BNames(1:end-8);
        load(fullfile(Pth3,BNames))
        B=B2Kegg.B;K=B2Kegg.K;
        modelname = xmlname(ismember(xmlname,btempname));
        model = readCbModel(fullfile(Pth33,[modelname{1},'.xml']));
        modelrxn = model.rxns;
        currnetkegg = netkegg(ismember(netbigg,modelrxn));
        metxn = numel(currnetkegg);
        bikegg = sum(ismember(modelrxn,B));
        corr(ct1).m = metxn; corr(ct1).bi = bikegg;
        corr(ct1).mRxns = netbigg(ismember(netbigg,modelrxn)); corr(ct1).biRxns = modelrxn(ismember(modelrxn,B));
        corr(ct1).bRxns = currnetkegg;
        names{ct1} = modelname{1};
        clear B2Kegg
        ct1 = ct1 + 1;
    end
end
