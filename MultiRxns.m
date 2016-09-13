function [Po,Hout1] = MultiRxns(K)
% MultiRxns
% gets a set of reactions and based on data provided in KEGG
% database for each reaction, searches for multi-step reactions and saves
% each multi-step reaction along with all of its single-step reactions
% (Child reactions).
% Inputs:
% K = A cell array of KEGG reactions identifiers. 

% Output:
% Po = A cell array of multi-step reactions found in K
% Hout1 = A cell array of decomposed single-step reactions corresponding to
%         each reaction in Po
% 
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Nov. 2015
% EDITED:
% Jan. 2016 : Offline mode was included for reducing the cpu time

NetStat = 0; % Offline mode
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth3 = fullfile(Pth,'BiGG2KEGG\Multirxns.mat');
Pth2 = fullfile(Pth,'OfflineRxnData');
if exist(Pth3,'file')
    load (Pth3)
    if exist('Multirxns','var')
        Crxns = Multirxns.C;
        Prxns = Multirxns.P;
        PoT = Prxns(ismember(Prxns,K));
        Hout1T = Crxns(ismember(Prxns,K));
        K = K(~ismember(K,Prxns));
    end
end
if isempty(K)
    Po = PoT;
    Hout1 = Hout1T;
    return
end
% Check internet connection
if NetStat
    TestLink='http://rest.kegg.jp';
    [~,Stat1]=urlread(TestLink);
    if ~Stat1
        msgbox({'Function: MultiRxns needs a stable internet connection';...
            'Seemingly you are offline!'});
        return
    end
end

ct=1;
for count = 1:numel(K)
    cts=1; setappdata(0,'cts',cts)
    if exist('Po','var') && any(strcmp(K{count},Po))
        continue % No need to check again
    end
    if NetStat
        A = urlread(['http://rest.kegg.jp/get/',K{count}]);
    else
        load([Pth2,'\',K{count},'.mat']);
        A = RnDat;
        clear RnDat
    end
    [~,c2]=regexp(A,'COMMENT\s');
    if ~isempty(c2)
        [c3,~]=regexp(A,'RPAIR');
        if isempty(c3)
            [c3,~]=regexp(A,'ENZYME');
        end
        Temp = A(c2:c3-1);
        Tip1 = regexp(Temp,'\w*-step\>');
        if ~isempty(Tip1)
            Temp2 = Temp(Tip1:end);
            c4 = regexp(Temp2,'R\d{4}');
            if ~isempty(c4)
                c5 = regexp(Temp2,'+');
                c54 = sort([c4,c5]);
                c5Loci = find(ismember(c54,c5));
                c5LociP = c5Loci+1; c5LociM = c5Loci-1;
                c5Loci1 = unique([c5LociP,c5LociM]);
                c6 = regexp(A,'part of', 'once');
                c7 = regexp(A,'second step of', 'once');
                c8 = regexp(A,'first step of', 'once');
                if ~isempty(c5) && isempty(c6) && isempty(c7) && isempty(c8)
                    Po{ct}=K{count};
                    for lop = 1:numel(c5Loci1)
                        Pol{ct}{lop} = Temp2(c54(c5Loci1(lop)):c54(c5Loci1(lop))+5);
                    end
                    Hout = SecLump(Pol{ct});
                    if ~isempty(Hout)
                        Hout1{ct} = SecCheck(Hout);
                        if ~isempty(Hout1{ct})
                            fprintf('Multi-step reaction was found:%d of %d:%s\n',...
                                count,numel(K),K{count});
                        end
                        ct=ct+1;
                    end
                end
            end
        end
        clear c2 c3 c4 c5 c6 c7 c8
    end
end
if exist('Po','var')
    Po = [Po,PoT];
else
    Po = PoT;
end
if exist('Hout1','var')
    Hout1 = [Hout1,Hout1T];
else
    Hout1 = Hout1T;
end

function Hout = SecLump(Hin)
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth2 = fullfile(Pth,'OfflineRxnData');
NetStat = 0;
ct=1;
for count = 1:length(Hin)
    if NetStat
        A = urlread(['http://rest.kegg.jp/get/',Hin{count}]);
    else
        load([Pth2,'\',Hin{count},'.mat']);
        A = RnDat;
        clear RnDat
    end
    [~,c2]=regexp(A,'COMMENT\s');
    if ~isempty(c2)
        [c3,~]=regexp(A,'RPAIR');
        if isempty(c3)
            [c3,~]=regexp(A,'ENZYME');
        end
        Temp = A(c2:c3-1);
        Tip1 = regexp(Temp,'\w*-dependent enzyme\>');
        if ~isempty(Tip1)
            Temp2 = Temp(Tip1:end);
            c4 = regexp(Temp2,'R\d{4}');
            for lop = 1:numel(c4)
                Pol{ct}{lop} = Temp2(c4(lop):c4(lop)+5);
            end
            Hout = SecLump(Pol{ct});
            ct = ct+1;
        else
            
            cts=getappdata(0,'cts');
            Hout{cts} = Hin{count};
            cts=cts+1;
            setappdata(0,'cts',cts) 
        end
    else
        cts=getappdata(0,'cts');
        Hout{cts} = Hin{count};
        cts=cts+1;
        setappdata(0,'cts',cts)
    end
end

function Hout1 = SecCheck (Hin1)
Fileid1 = fopen('rxn2map.txt','r');
rxn2map = textscan(Fileid1,'%s %s');
RawRxns = rxn2map{1};
RawMaps = rxn2map{2};
[Loci1,~] = regexp(RawMaps, 'path:rn');
Loci2 = ~cellfun('isempty', Loci1);
RawMaps(Loci2) = [];
RawRxns(Loci2) = [];
Maps = strrep(RawMaps,'path:map','');
Rxns = strrep(RawRxns,'rn:','');
MapsTemp={0};
for count = 1:length(Hin1)
    MapsTemp{count} = Maps(ismember(Rxns,Hin1{count}))';
    if isempty(MapsTemp{count})
        Hout1 = [];
        break
    end
end
if ~exist('Hout1','var')
    MapsT = [MapsTemp{:}];
    MapsTU = unique(MapsT);
    if numel(MapsT)>numel(MapsTU)
        Hout1 = Hin1;
    end
end