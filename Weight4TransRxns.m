function Weight4TransRxns(Idfier)
% Weight4TransRxns
% This function generates the Weight4TranRxns.txt file containing all KEGG
% reactions in which there are 2,4 or 6 KEGG compounds, along with their
% corresponding component weights. This text file will be used in Bigg2Kegg
% function to identify potential equivalent KEGG reactions for BiGG
% transport reactions.
% Input:
%     - Idfier : A number: 2 | 4 | 6 , respectively denoting transport
%     reactions with 2, 4 and 6 components.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Jan. 2016

Fileid = fopen('cpd2rxn.txt','r');
cpd2rxn = textscan(Fileid,'%s %s');
rxnData = strrep(cpd2rxn{2},'rn:','');
cpds = strrep(cpd2rxn{1},'cpd:','');
fclose(Fileid);
% Read cpd2weight ===================
Fileid = fopen('cpd2weight.txt','r');
cpd2weight = textscan(Fileid,'%s %f');
cpd4weight = strrep(cpd2weight{1},'cpd:','');
Weight = cpd2weight{2};
fclose(Fileid);
CurLen = 1; W4CurCpds = ({}); CanRxn = ({});
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);

if (Idfier == 2)
    fprintf('Generating W4TranRxns2.txt for rxn: ')
    Pth1 = [Pth,'\','W4TranRxns2.txt'];
    Fid = fopen(Pth1,'w');
elseif (Idfier == 4)
    fprintf('Generating W4TranRxns4.txt for rxn: ')
    Pth1 = [Pth,'\','W4TranRxns4.txt'];
    Fid = fopen(Pth1,'w');
elseif (Idfier == 6)
    fprintf('Generating W4TranRxns6.txt for rxn: ')
    Pth1 = [Pth,'\','W4TranRxns6.txt'];
    Fid = fopen(Pth1,'w');
end
for count = 1:numel(rxnData)
    if count>1
        for j=0:log10(count-1)
            fprintf('\b'); 
        end
    end
    fprintf('%d', count);
    TNum = sum(ismember(rxnData,rxnData{count}));
    if TNum == Idfier
        CurCpds = cpds(ismember(rxnData,rxnData{count}));
        W4CurCpdsT = Weight(ismember(cpd4weight,CurCpds));
        if numel(W4CurCpdsT) == Idfier && ~any(ismember(CanRxn,rxnData{count})) ...
                && (2.*numel(unique(W4CurCpdsT)) == Idfier)
            W4CurCpds{CurLen} = W4CurCpdsT;
            CanRxn{CurLen} = rxnData{count};
            CurLen = CurLen + 1;
        end       
    end
end

if (Idfier == 2)
    for ct = 1:numel(W4CurCpds)
        fprintf(Fid,'%s\t%f\t%f\n',CanRxn{ct},W4CurCpds{ct});
    end
elseif (Idfier == 4)
    for ct = 1:numel(W4CurCpds)
        fprintf(Fid,'%s\t%f\t%f\t%f\t%f\n',CanRxn{ct},W4CurCpds{ct});
    end
elseif (Idfier == 6)
    for ct = 1:numel(W4CurCpds)
        fprintf(Fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\n',CanRxn{ct},W4CurCpds{ct});
    end
end
fclose(Fid);
