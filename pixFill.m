function [I,pixx,pixy,pixx1,pixy1] = pixFill(I,pixx,pixy,LineStrength,ColSpec)
% pixFill
% is a subfunction of NetDraw for generating pixel data points
% based on the specified line width for reactions. 
% 
% Inputs:
% I: An array containing the customized metabolic map.
% pixx, pixy: Reaction coordinates.
% ColSpec: Color code of the specified reaction.
% 
% Outputs:
% I: Customized metabolic map after adjustments.
% pixx,pixy,pixx1,pixy1: Boundary pixel data based on the specified
% reactions' width.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016

if ~std(pixx)
    pixx1 = pixx + LineStrength/2;
    pixx = pixx - LineStrength/2;
    pixy1 = pixy;
elseif ~std(pixy)
    pixy1 = pixy + LineStrength/2;
    pixy = pixy - LineStrength/2;
    pixx1 = pixx;
else
    pixx1 = pixx + LineStrength/2;
    pixx = pixx - LineStrength/2;
    pixy1 = pixy + LineStrength/2;
    pixy = pixy - LineStrength/2;
end
pixx1 = round(pixx1); pixx = round(pixx);
pixy1 = round(pixy1); pixy = round(pixy);
for i5 = 1:numel(pixx)
    I(pixy(i5):pixy1(i5),pixx(i5),1) = ColSpec(1);
    I(pixy(i5),pixx(i5):pixx1(i5),1) = ColSpec(1); 
    I(pixy(i5):pixy1(i5),pixx(i5),2) = ColSpec(2);
    I(pixy(i5),pixx(i5):pixx1(i5),2) = ColSpec(2); 
    I(pixy(i5):pixy1(i5),pixx(i5),3) = ColSpec(3);
    I(pixy(i5),pixx(i5):pixx1(i5),3) = ColSpec(3); 
end