function UniModel
% UiModel
% generates one unified model containing all shared reaction among
% all BiGG models in Bigg2Kegg folder. This model is of use in cases where
% user wants not to use a special BiGG model, or a new model is in hand.
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
ct1 = 1; NewB = ({}); NewK = ({});
for ct = 1:numel(BMatfiles)
    if ~strcmp(BMatfiles(ct).name,'BiGG2KEGG_HMRbased-RECON1.mat') && ...
            ~strcmp(BMatfiles(ct).name,'Multirxns.mat')
        BNames{ct1} = BMatfiles(ct).name;
        load(fullfile(Pth3,BNames{ct1}))
        B=B2Kegg.B;K=B2Kegg.K;
        Mlt = find(ismember(B,'MULTIR'));
        if ~isempty(Mlt)
            TB = B(1:Mlt-1); TK = K(1:Mlt-1);
        end
        if isempty(NewB)
            NewB = TB; NewK = TK;
        else
            NotB = TB(~ismember(TB,NewB));
            NotK = TK(~ismember(TB,NewB));
            NewLen = numel(NewB)+1:numel(NotB)+numel(NewB);
            NewB(NewLen) = NotB; NewK(NewLen) = NotK;
        end
        clear B2Kegg
        ct1 = ct1 + 1;
    end
end
Len1 = numel(NewB);
NewB{Len1+1} = 'MULTIR';
NewK{Len1+1} = 'MULTIR';
load(fullfile(Pth3,'Multirxns.mat'))
Pr = Multirxns.P; Cr = Multirxns.C;
clear Multirxns
for c1 = 1:numel(Pr)
    TempN = NewB(ismember(NewK,Pr{c1}));
    TempC = Cr{c1};
    Len = numel(TempC);
    for c2 = 1:numel(TempN)
        Len1 = numel(NewB);
        Len2 = Len1 + 1: Len1 + Len;
        NewB(Len2) = TempN(c2);
        NewK(Len2) = TempC;
    end
end
B2Kegg.B = NewB; B2Kegg.K = NewK;
fname = fullfile(Pth3,'UniModelKEGG.mat');
save (fname,'B2Kegg')