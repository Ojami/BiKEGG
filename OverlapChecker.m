function [I,ArX,ArY,ArrowChek,DoubleData] = OverlapChecker(I,ArX,ArY,OverlapData,...
    flx,McmCollect,field_intsect_overlay,LineStrength)
% OverlapChecker
% is a subfunction of NetDraw and identifies overlapping reactions.
% 
% Inputs:
% I: An array containing the customized metabolic map.
% ArX,ArY: Coordinates for flux carrying reactions.
% OverlapData: Data for identifying overlapping reactions.
% flx: Flux rates.
% McmCollect: Color codes corresponding to flux rates.
% field_intsect_overlay: Reaction identifiers for flux carrying reactions.
% LineStrength: Reaction's line width.
% 
% Outputs:
% I: Customized metabolic map after modification for overlapping reactions.
% ArX,ArY: Coordinates for flux carrying reactions after modification for
% overlapping reactions.
% ArrowChek: Data for arrowhead function.
% DoubleData: Data for ArrowheadProcess function.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016

[~,Idx1,~] = unique(OverlapData);
Chck1 = setdiff(1:size(OverlapData,1),Idx1);
ArrowChek = []; ct = 1; RxnId = ({}); RxnX = ({}); RxnY = ({});
RxnX1 = RxnX; RxnY1 = RxnY;
for m1 = 1:numel(Chck1)
    WhrC = find(ismember(OverlapData,OverlapData{Chck1(m1)}));
    if std(flx(WhrC))<1e-4 % Same flux rates
        continue
    end
    pixx = ArX{WhrC(1)}; pixy = ArY{WhrC(1)};
    ColSpec = McmCollect(WhrC,:);
    % set start from the most top/left part of line
    if ~std(pixx)
            pixx = pixx - LineStrength*numel(WhrC)/2;
    elseif ~std(pixy)
            pixy = pixy - LineStrength*numel(WhrC)/2;
    else
            pixx = pixx - LineStrength*numel(WhrC)/2;
            pixy = pixy - LineStrength*numel(WhrC)/2;
    end
%     LineSegments = linspace(0,LineStrength,numel(WhrC)+1);
%     LineSegments = sort(LineSegments,'descend');
    LineSegments = [];
    LineSegments(1) = 0; LineSegments(2) = LineStrength;
    for i1 = 3:numel(WhrC)*2
        if mod(i1,2)
            LineSegments(i1) = LineSegments(i1-1) + 1;
        else
            LineSegments(i1) = LineSegments(i1-1) + LineStrength;
        end
    end
    
    for m2 = 1:numel(WhrC)
        if ~std(pixx)
            pixx1 = pixx + LineSegments(2*m2-1);
            pixx1A = pixx1 + LineStrength/2;
            pixx2 = pixx + LineSegments(2*m2);
            pixy1 = pixy;
            pixy1A = pixy;
            pixy2 = pixy;
        elseif ~std(pixy)
            pixy1 = pixy + LineSegments(2*m2-1);
            pixy1A = pixy1 + LineStrength/2;
            pixy2 = pixy + LineSegments(2*m2);
            pixx1 = pixx;
            pixx1A = pixx;
            pixx2 = pixx;
        else
            pixx1 = pixx + LineSegments(2*m2-1);
            pixx1A = pixx1 + LineStrength/2;
            pixx2 = pixx + LineSegments(2*m2);
            pixy1 = pixy + LineSegments(2*m2-1);
            pixy1A = pixy1 + LineStrength/2;
            pixy2 = pixy + LineSegments(2*m2);
        end
        pixx1 = round(pixx1); pixx2 = round(pixx2);
        pixy1 = round(pixy1); pixy2 = round(pixy2);
        ArX{WhrC(m2)} = []; ArX{WhrC(m2)} = pixx1A; RxnId{ct} = field_intsect_overlay{WhrC(m2)};
        ArY{WhrC(m2)} = []; ArY{WhrC(m2)} = pixy1A; 
        RxnX{ct} =  pixx1; RxnY{ct} = pixy1; RxnX1{ct} = pixx2;RxnY1{ct} = pixy2;
        ArrowChek(ct,1) = WhrC(m2); ArrowChek(ct,2) = numel(WhrC);
        ct = ct + 1;
        for i5 = 1:numel(pixx)
            I(pixy1(i5):pixy2(i5),pixx1(i5),1) = ColSpec(m2,1);
            I(pixy1(i5),pixx1(i5):pixx2(i5),1) = ColSpec(m2,1); 
            I(pixy1(i5):pixy2(i5),pixx1(i5),2) = ColSpec(m2,2);
            I(pixy1(i5),pixx1(i5):pixx2(i5),2) = ColSpec(m2,2); 
            I(pixy1(i5):pixy2(i5),pixx1(i5),3) = ColSpec(m2,3);
            I(pixy1(i5),pixx1(i5):pixx2(i5),3) = ColSpec(m2,3); 
        end
    end
end
DoubleData.X = RxnX; DoubleData.Y = RxnY; DoubleData.Id = RxnId;
DoubleData.X1 = RxnX1; DoubleData.Y1 = RxnY1;