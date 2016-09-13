function I = flx2col(BoxCol,flx,Cm,I,k,flxMax,flxMin,m,idx,idw,Tempx,Thresh)
% flx2col
% Performs color mapping for caculated flux values and normalizes them based on
% maximum and minimum flux values in the input data.
% This function uses following inputs:
% 1- BoxCol: RGB code entered by user for background color of each box
% 2- flx : Contains flux values in form of:
%       2-1 A 2D matrix, in which each column represents flux values for
%           each state of one certain species/model (i.e. dynamic state 
%           of a certain model extracted from COBRA Toolbox). In this
%           situation flx may contain various columns, indicating different
%           circumstances of a certain species.
%       2-2 A 2D matrix, in which each column represents flux valuse for
%           each state of two different species/models (similarly, from
%           COBRA Toolbox). In this situation, if a dynamic comparison
%           between two species/models is performed, the number of flx
%           columns must be even. 
%       2-3 A column vector. This situation is similar to case 1-1, but in
%           static state.
% 3- Cm : A colormap based on MATLAB built-in colormaps.
% 4- I : The output of imread, containing RGB values of KEGG map pixels.
% 5- k, m : Column and row indices, respectively.
% 6- idx : Coordinates of each enzyme box on a KEGG map image, extracted
%          from corresponding KGML file
% 7- idw : Width of each box corresponding to idx, extracted from KGML
%          file
% 8- Tempx : Indicates the height of each box which should be colored
% 9- Thresh: Threshold for each flux value, which lower than this threshold
%            the box will be left blank, suggesting that the flx value is
%            near zero.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Feb. 2015
% Modified: Lines 70-78 for zero carrying fluxes
% =========================================================================

% Normalization of color indices in Cm, Red for highest flux and
% Blue for lowest flux. 
flx = abs(flx);
% flxMax=max(abs(flx(:,k)));
% flxMin=min(abs(flx(:,k)));
% if flxMax==flxMin % One row (i.e. for 1 rxn)
%     flxMax=max(abs(flx(m,:)));
%     flxMin=min(abs(flx(m,:)));
% end
    flxSl=(size(Cm,1)-1)./(flxMax-flxMin);
    if isinf(flxSl) % in case of flxMax == flxMin
        flxSl = 1e15.*(size(Cm,1)-1);
    end
    Mcm=ceil((flx(m,k)-flxMax).*flxSl+size(Cm,1));
    Tempy=(round(idx(m)+1):round(idx(m)+idw(m)-1));
if flxMax > 0   
    [x1,y1]=find(I(Tempx,Tempy,1)==0 & I(Tempx,Tempy,2)==0 & I(Tempx,Tempy,3)==0);
    Tempy1=round(idx(m)+1):round(idx(m)+(flx(m,k)*idw(m))./flxMax-1);
    Tempy2=round(idx(m)+(flx(m,k)*idw(m))./flxMax)-1:(idx(m)+idw(m)-1);
    [x2,y2]=find(I(Tempx,Tempy1,1)==0 & I(Tempx,Tempy1,2)==0 & I(Tempx,Tempy1,3)==0);
    I(Tempx,Tempy1,1)=Cm(Mcm,1);
    I(Tempx,Tempy1,2)=Cm(Mcm,2);
    I(Tempx,Tempy1,3)=Cm(Mcm,3);
    if isempty(BoxCol) % Default color
        I(Tempx,Tempy2,1)=169;%Dark gray for background of each box
        I(Tempx,Tempy2,2)=169;
        I(Tempx,Tempy2,3)=169;
    else % User defined color
        I(Tempx,Tempy2,1)=BoxCol(1);
        I(Tempx,Tempy2,2)=BoxCol(2);
        I(Tempx,Tempy2,3)=BoxCol(3);
    end
end
if abs(flx(m,k))< Thresh % Threshold for considering a flux = 0
    if isempty(BoxCol) % Default color
        for g1=1:length(x1)
            I(Tempx(x1(g1)),Tempy(y1(g1)),:)=169;
        end
    else
        I(Tempx,Tempy,1)=BoxCol(1);
        I(Tempx,Tempy,2)=BoxCol(2);
        I(Tempx,Tempy,3)=BoxCol(3);
    end
else
    for g1=1:length(x2)
        I(Tempx(x2(g1)),Tempy1(y2(g1)),1)=Cm(Mcm,1);
        I(Tempx(x2(g1)),Tempy1(y2(g1)),2)=Cm(Mcm,2);
        I(Tempx(x2(g1)),Tempy1(y2(g1)),3)=Cm(Mcm,3);
    end
end
end