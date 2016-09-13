function basemap = BaseMapReader(Netstat)
% BaseMapReader
% reads KGML file of KEGG global metaboli pathway (map01100),
% extracts all required information and organizes them in a structure
% array, with three fields of r, c, and rc for respectively reactions, 
% compounds and reaction-to-compound details for subsequent use in NetDraw.
% 
% Input:
% Netstat: Operate in online (1) or offline (0) mode. In offline mode, data
% in KEGGmaps folder of BiKEGG will be used.
% 
% Output:
% basemap: Structure array with three fields of r, c, and rc for reactions,
% compounds and reaction-to-compound details respectively. 

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016
% -------------------------------------------------------------------------

if Netstat
    Dat = 'http://rest.kegg.jp/get/rn01100/kgml';
    Dat = urlread(Dat);
else
    load(which('map01100.mat'))
end
% Extract rxns details ====================================================
rxnDat = regexp(Dat,'(?<=entry) id="(\d+)" name="rn:[^"]*.*?(?=</entry>)','match');
RawRNs = regexp(rxnDat,'(?<=reaction="rn:)[^"]*','match');
Lineseg = regexp(rxnDat, '(?<=type="line" coords=")[^"]*', 'match');
rxns = ({}); ctr = 1; rxnx = ({}); rxny = ({}); rxn_id = (0); rxn_col = ({});
for i1 = 1:numel(RawRNs)
    temp_rxn = regexp(RawRNs{i1}{1}, 'R\d{5}', 'match');
    rxn_idTemp = str2double(regexp(rxnDat{i1},'(?<=id=")(\d+)[^"]*','match'));
    rxn_colTemp = unique(regexp(rxnDat{i1},'(?<=fgcolor="#)[^"]*','match'));
    temp_coord =regexp(Lineseg{i1}, '\d+', 'match');
    line_coord1 = cellfun(@str2double, temp_coord, 'UniformOutput', false);
    line_x = ({}); line_y = ({});
    for i3 = 1:numel(line_coord1) % Extract all lines for each rxn
        temp_coord1 = line_coord1{i3};
        x_temp = temp_coord1(1:2:end);
        y_temp = temp_coord1(2:2:end);
        line_x{i3} = x_temp;
        line_y{i3} = y_temp;
        clear x_temp y_temp
    end
    for i2 = 1:numel(temp_rxn)
        rxns{ctr} = temp_rxn{i2}; % All rxn details !!!!!!!!!!!!!
        rxn_id(ctr) = rxn_idTemp;
        rxn_col(ctr) = rxn_colTemp;
        rxnx{ctr} = line_x;  % All rxn details !!!!!!!!!!!
        rxny{ctr} = line_y;  % All rxn details !!!!!!!!!!!
        ctr = ctr + 1;
    end
end
clear rxnDat RawRNs Lineseg
% =========================================================================
% Extract compound details ================================================
cpdDat = regexp(Dat,'(?<=entry) id="(\d+)" name="cpd:[^"]*.*?(?=</entry>)','match');
glDat = regexp(Dat,'(?<=entry) id="(\d+)" name="gl:[^"]*.*?(?=</entry>)','match');
cpd_coord = (0); cpd_id = (0); cpds = ({});
for i1=1:numel(cpdDat)
    cpds(i1) = regexp(cpdDat{i1},'(?<=name="cpd:)[^"]*','match');
    cpd_id(i1) = str2double(regexp(cpdDat{i1},'(?<=id=")(\d+)[^"]*','match'));
    cpd_coords = regexp(cpdDat{i1},'type="circle" x="(\d+)" y="(\d+)" width="(\d+)" height="(\d+)"','tokens');
    cpd_coord(i1,1) = str2double(cpd_coords{1}{1});
    cpd_coord(i1,2) = str2double(cpd_coords{1}{2});
    cpd_coord(i1,3) = str2double(cpd_coords{1}{3});
    cpd_coord(i1,4) = str2double(cpd_coords{1}{4});
end
gl_coord = (0); gl_id = (0); gls = ({});
for i1 = 1:numel(glDat)
    gls(i1) = regexp(glDat{i1},'(?<=name="gl:)[^"]*','match');
    gl_id(i1) = str2double(regexp(glDat{i1},'(?<=id=")(\d+)[^"]*','match'));
    gl_coords = regexp(glDat{i1},'type="circle" x="(\d+)" y="(\d+)" width="(\d+)" height="(\d+)"','tokens');
    gl_coord(i1,1) = str2double(gl_coords{1}{1});
    gl_coord(i1,2) = str2double(gl_coords{1}{2});
    gl_coord(i1,3) = str2double(gl_coords{1}{3});
    gl_coord(i1,4) = str2double(gl_coords{1}{4});
end
% =========================================================================
% Extract rxn-to-compound details =========================================
% NOTE: rxn-to-compound details are collected based on rxn ID, meaning all
% rxns correlated with one certain ID are sotred together.
rxn4cmp = regexp(Dat,'(?<=<reaction)\s+id="(\d+)" (.+?/>)*?\s*(?=</reaction>)','tokens');
rxn4cmp_id = (0); rxn4cmp_subid= ({}); rxn4cmp_subname= ({});
rxn4cmp_prodid= ({}); rxn4cmp_prodname= ({}); rxn4cmp_name = ({});
for i1 = 1:numel(rxn4cmp)
    rxn4cmp_det = regexp(rxn4cmp{i1}{2},'(?<="rn:)[^"]*','match');
    temp_rxn1 = regexp(rxn4cmp_det{1}, 'R\d{5}', 'match');
    for r1 = 1:numel(temp_rxn1)
        rxn4cmp_name{i1}{r1} = temp_rxn1{r1};
    end
    rxndet_sub = regexp(rxn4cmp{i1}{2},'(?<=<substrate)\s+id="(\d+)"\s+name="cpd:C(\d+)"[^/]*','tokens');
    rxndet_sub1 = regexp(rxn4cmp{i1}{2},'(?<=<substrate)\s+id="(\d+)"\s+name="gl:G(\d+)"[^/]*','tokens');
    rxndet_prod = regexp(rxn4cmp{i1}{2},'(?<=<product)\s+id="(\d+)"\s+name="cpd:C(\d+)"[^/]*','tokens');
    rxndet_prod1 = regexp(rxn4cmp{i1}{2},'(?<=<product)\s+id="(\d+)"\s+name="gl:G(\d+)"[^/]*','tokens');
    for i2 = 1:numel(rxndet_sub)
        rxn4cmp_subid{i1}(i2) = str2double(rxndet_sub{i2}(1));
        rxn4cmp_subname{i1}{i2} = strcat('C',rxndet_sub{i2}{2});
    end
    for i2 = numel(rxndet_sub)+1:numel(rxndet_sub1)+numel(rxndet_sub)
        rxn4cmp_subid{i1}(i2) = str2double(rxndet_sub1{i2-numel(rxndet_sub)}(1));
        rxn4cmp_subname{i1}{i2} = strcat('G',rxndet_sub1{i2-numel(rxndet_sub)}{2});
    end
    for i3 = 1:numel(rxndet_prod)       
        rxn4cmp_prodid{i1}(i3) = str2double(rxndet_prod{i3}(1));
        rxn4cmp_prodname{i1}{i3} = strcat('C',rxndet_prod{i3}{2});       
    end
    for i3 = 1+numel(rxndet_prod):numel(rxndet_prod1)+numel(rxndet_prod)      
        rxn4cmp_prodid{i1}(i3) = str2double(rxndet_prod1{i3-numel(rxndet_prod)}(1));
        rxn4cmp_prodname{i1}{i3} = strcat('G',rxndet_prod1{i3-numel(rxndet_prod)}{2});
    end
    rxn4cmp_id(i1) = str2double(rxn4cmp{i1}{1});
end
basemap.r.rxn = rxns;
basemap.r.rxnid = rxn_id;
basemap.r.col = rxn_col;
basemap.r.rxnx = rxnx;
basemap.r.rxny = rxny;
basemap.c.cpd = cpds;
basemap.c.cpdid = cpd_id;
basemap.c.cpdxy = cpd_coord;
basemap.c.gl = gls;
basemap.c.glid = gl_id;
basemap.c.glxy = gl_coord;
basemap.rc.rxn = rxn4cmp_name;
basemap.rc.rxnid = rxn4cmp_id;
basemap.rc.sub = rxn4cmp_subname;
basemap.rc.subid = rxn4cmp_subid;
basemap.rc.prod = rxn4cmp_prodname;
basemap.rc.prodid = rxn4cmp_prodid;

% Modify basemap: There are some discrepancies between KGML and basemap
% Error in R08940
basemap.rc.sub{basemap.rc.rxnid==18} = 'C00019';
basemap.rc.subid{basemap.rc.rxnid==18} = 4254;
% Error in R04097
basemap.rc.prod{basemap.rc.rxnid==718} = 'C15975';
basemap.rc.prodid{basemap.rc.rxnid==718} = 2997;