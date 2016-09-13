function GetData4MultiRxns
% GetData4MultiRxns
% Gets offline data from KEGG database for each reaction in
% the KEGG for later use in MultiRxn function.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Jan. 2016

RawRns = urlread('http://rest.kegg.jp/list/rn');
Rns = regexp(RawRns, 'R\d{5}', 'match'); % All KEGG rxns
Pth1 = which ('Bigg2Kegg.m');
tind = find(Pth1=='\',1,'last');
Pth = Pth1(1:tind-1);
Pth = fullfile(Pth,'OfflineRxnData');
fprintf('\nDownloading rxn data (of %d): ',numel(Rns))
for ct = 1:numel(Rns)
    if ct>1
        for j=0:log10(ct-1)
            fprintf('\b'); 
        end
    end
    fprintf('%d', ct);
    RnDat = urlread(['http://rest.kegg.jp/get/',Rns{ct}]);
    save([Pth,'\',Rns{ct},'.mat'],'RnDat')
    % In case of using txt file format
%     Pth1 = [Pth,'\',Rns{ct},'.txt'];
%     Idf = fopen(Pth1,'w');
%     fprintf(Idf,RnDat);
%     fclose(Idf);
end
DateCrt = date;
save([Pth,'\','VersionDate.mat'],'DateCrt')
fprintf('\n')
disp('All reaction data has been downloaded and saved successfully!')