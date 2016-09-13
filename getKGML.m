function getKGML
% getKGML
% reads and retrieves KGML file and image data of KEGG pathways,
% and stores them in "KEGGmaps" subfolder within FluxPuya folder.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Dec. 2015

Fileid2 = fopen('KEGGmaps.txt','r');
TempMaps = textscan(Fileid2,'%s %[^\n]');
TempMaps1 = cell(size(TempMaps{1},1),1);
for count = 1:size(TempMaps{1},1)
    TempMaps1{count}=[TempMaps{1}{count},' ',TempMaps{2}{count}];
end
TempMaps1 = strrep(TempMaps1,'path:map','');
MapNme = regexp(TempMaps1,'\d{5}','match');
MapNme = [MapNme{:}];
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
ctYes = 1; ctNo = 1;
for count = 1:numel(MapNme)
    Link=['http://rest.kegg.jp/get/rn',MapNme{count},'/kgml'];
    ImgLink=['http://www.kegg.jp/kegg/pathway/map/map',MapNme{count},'.png'];
    I=imread(ImgLink);
    [Dat,Stat1]=urlread(Link);
    if ~Stat1
        fprintf('%d-Error:No data was found for map%s\n',ctNo,MapNme{count});
        ctNo = ctNo + 1;
        continue
    else
        save([Pth,'\KEGGmaps\',MapNme{count},'.mat'],'Dat','I');
        fprintf('%d-Map%s data was saved successfully to KEGGmaps folder\n'...
            ,ctYes,MapNme{count});
        ctYes = ctYes + 1;
    end
end
DateCrt = date;
save([Pth,'\KEGGmaps\VersionDate.mat'],'DateCrt')