function pixFit(x,y,rxn)
% pixFit 
% uses discrete position coordinates from BaseMapReader and interpolates
% them to generate spatially continuous points for each reaction. This
% function uses the already existing Pixdata.mat file.
% 
% Inputs:
% x,y: Original discrete coordinates from BaseMapReader.
% rxn: Reaction identifier in the form of [KEGG reaction ID,Entry ID provided
% in the global map's KGML file], e.g.: 'R000933116' is a unique reaction
% identifier for reaction ID R00093 and entry ID 3116.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016

if size(x,2) > size(x,1)
    x = x';
end
if size(y,2) > size(y,1)
    y = y';
end
xy = [x,y];
xy = unique(xy,'rows');
x = xy(:,1); y = xy(:,2);

if ~all(diff(x)==0) % Not a vertical line
    pixxy = interparc(2000,x,y,'pchip'); % Copyright (c) 2012, John D'Errico
    pixx = pixxy(:,1);
    pixy = pixxy(:,2);
    pixy = round(pixy);
    pixx = round(pixx);
    if size(pixx,2) > size(pixx,1)
        pixx = pixx';
    end
    if size(pixy,2) > size(pixy,1)
        pixy = pixy';
    end
    pixxy = [pixx,pixy];
    pixxy = unique(pixxy,'rows');
    pixx = pixxy(:,1); pixy = pixxy(:,2);   
else
    pixy = min(y):max(y);
    pixx = ones(1,numel(pixy)).*x(1);
    if size(pixx,2) > size(pixx,1)
        pixx = pixx';
    end
    if size(pixy,2) > size(pixy,1)
        pixy = pixy';
    end
    pixxy = [pixx,pixy];
    pixxy = unique(pixxy,'rows');
    pixx = pixxy(:,1); pixy = pixxy(:,2);
end

path_name = which('Pixdata.mat');
if isempty(path_name)
    msgbox({'Pixdata.mat cannot be found!';...
                'Make sure the file is present in BiKEGG folder'},'Error','error');
    return
end
load(path_name)

if isfield(pix_x,rxn)
    fprintf('%s has been replaced with new pixels!\n',rxn)
    temp_pixx = pix_x.(rxn);
    temp_pixy = pix_y.(rxn);
    if size(pixx,2) > size(pixx,1)
        pixx = pixx';
    end
    if size(pixy,2) > size(pixy,1)
        pixy = pixy';
    end
    if size(temp_pixx,2) > size(temp_pixx,1)
        temp_pixx = temp_pixx';
    end
    if size(temp_pixy,2) > size(temp_pixy,1)
        temp_pixy = temp_pixy';
    end
    pixx = [pixx;temp_pixx];
    pixy = [pixy;temp_pixy];
    pixxy = [pixx,pixy];
    pixxy = unique(pixxy,'rows');
    pixx = pixxy(:,1); pixy = pixxy(:,2);
    pix_x.(rxn) = pixx;
    pix_y.(rxn) = pixy;
else
    pix_x.(rxn) = pixx;
    pix_y.(rxn) = pixy;
end
delete(path_name)
save(path_name,'pix_x','pix_y')

