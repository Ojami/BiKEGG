function MultiRxnsLib
% MultiRxnsLib
% generates shared multi-step reactions among all models saved
% in Bigg2Kegg folder. This list is used in both function UniModel and
% MultiRxns
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Jan. 2015

Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth2 = fullfile(Pth,'BiGG2KEGG\*.mat');
Pth3 = fullfile(Pth,'BiGG2KEGG');
BMatfiles = dir(Pth2);
BNames = cell(numel(BMatfiles)-2,1);
ct1 = 1;
for ct = 1:numel(BMatfiles)
    if ~strcmp(BMatfiles(ct).name,'BiGG2KEGG_HMRbased-RECON1.mat') && ...
            ~strcmp(BMatfiles(ct).name,'Multirxns.mat')
        BNames{ct1} = BMatfiles(ct).name;
        ct1 = ct1 + 1;
    end
end
BM = ({}); NotBM = ({}); NotKM = ({});
for count =1:numel(BNames)
    load(fullfile(Pth3,BNames{count}))
    B=B2Kegg.B;K=B2Kegg.K;
    MltLoci = find(ismember(B,'MULTIR'));
    NotBM1 = B(1:MltLoci-1); NotKM1 = K(1:MltLoci-1);
    BM1 = B(MltLoci+1:end); 
    NonOverlp = ~ismember(BM1,BM);
    BM1 = BM1(NonOverlp);
    [NotBM1,NotKM1] = Rxnfinder(BM1,NotBM1,NotKM1);
    Len = numel(BM)+1:numel(BM1)+numel(BM);
    clear B2Kegg
    BM(Len) = BM1; 
    if ~isempty(NotBM1)
        Len1 = numel(NotBM)+1:numel(NotBM1)+numel(NotBM);
        NotBM(Len1) = NotBM1; NotKM(Len1) = NotKM1;
    end
end
[ParentRxn,ChildRxns] = MultiRxns(NotKM);
Multirxns.P = ParentRxn; Multirxns.C = ChildRxns;
save ([Pth3,'\Multirxns.mat'],'Multirxns')
clear
disp('Multi Reactions Library has been generated successfully!')
%--------------------------------------------------------------------------
function [NotBM,NotKM] =  Rxnfinder(BM1,NotBM1,NotKM1)
count = 1; NotBM = ({}); NotKM = ({});
while numel(BM1)
    Tfind = find(ismember(BM1,BM1{1}));
    Tfind1 = find(ismember(NotBM1,BM1{Tfind(1)}));
    for ct = 1:numel(Tfind1)
        NotBM{count} = NotBM1{Tfind1(ct)};
        NotKM{count} = NotKM1{Tfind1(ct)};
        count = count + 1;
    end
    BM1(Tfind) = [];
end