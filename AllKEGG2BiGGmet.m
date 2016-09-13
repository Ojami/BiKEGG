function AllKEGG2BiGGmet
% AllKEGG2BiGGmet
% Extracts all KEGG metabolites from BIGG API and save them
% as a mat file (AllKEGG2BiGGmet.mat), which is used in MapAdjuster for
% displaying compound details on created customized maps by NetDraw.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016
%--------------------------------------------------------------------------

Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Fileid3 = fopen(fullfile(Pth,'bigg_models_metabolites.txt'),'r');
biggmets = textscan(Fileid3,'%s %s %[^\n]');
fclose(Fileid3);
biggmetid = biggmets{1}; biggmetdes = biggmets{3}; clear biggmets
Metkegg = ({}); Metbigg = ({});
ct3 = 1;
for ct1 = 2:numel(biggmetid)
    temploci = regexp(biggmetdes{ct1},...
        'http://identifiers.org/kegg.compound/C\S{5}','match');
    if isempty(temploci)
        continue
    end
    if iscell(temploci{1})
        temploci = temploci{1};
    end
    if isempty(temploci)
        continue
    else
        for ct2 = 1:numel(temploci)
            regtemp = regexp(temploci{ct2},'C\S{5}','match');
            Metkegg{ct3}{ct2} = regtemp{1};
            Metbigg{ct3}{ct2} = biggmetid{ct1};
        end
        ct3 = ct3 + 1;
    end
end
for ct1 = 2:numel(biggmetid) % KEGG glycan
    temploci = regexp(biggmetdes{ct1},...
        'http://identifiers.org/kegg.compound/G\S{5}','match');
    if isempty(temploci)
        continue
    end
    if iscell(temploci{1})
        temploci = temploci{1};
    end
    if isempty(temploci)
        continue
    else
        for ct2 = 1:numel(temploci)
            regtemp = regexp(temploci{ct2},'G\S{5}','match');
            Metkegg{ct3}{ct2} = regtemp{1};
            Metbigg{ct3}{ct2} = biggmetid{ct1};
        end
        ct3 = ct3 + 1;
    end
end
AllKEGG2BiGGmet.Metkegg = Metkegg; AllKEGG2BiGGmet.Metbigg = Metbigg;
save('AllKEGG2BiGGmet.mat','AllKEGG2BiGGmet')