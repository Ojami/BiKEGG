function I = MapTrimmer(I,ArrowX,ArrowY,Cntrs,Allx4Arrow,Ally4Arrow,LineStrength,ColSpec)
% MapTrimmer
% is a subfunction of NetDraw for drawing removing extra pixels
% after attaching arrowheads to reaction ends. 
% 
% Inputs:
% I: An array containing the customized metabolic map.
% ArrowX,ArrowY,Cntrs: Outputs from ArrowheadProcess.
% Allx4Arrow,Ally4Arrow: All reactions' coordinates.
% LineStrength: Reaction width.
% ColSpec: Color code corresponding to that of the reactions.
% 
% Output:
% I: Modified image without extra pixels.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016

for i1 = 1:numel(ArrowX)
    if isempty(ArrowX{i1})
        continue
    end
    Xin = ArrowX{i1}; Yin = ArrowY{i1}; Cents = Cntrs{i1};
    pixx = Allx4Arrow{i1};
    pixy = Ally4Arrow{i1};
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
    
    for cnt = 1:size(Xin,1)
        X = Xin(cnt,:); Y = Yin(cnt,:); Cent = Cents(cnt,:);
        DistToCpd = hypot( Allx4Arrow{i1} - Cent(1), Ally4Arrow{i1} - Cent(2));
        CpdDistToPoint1 = hypot( X(1) - Cent(1), Y(1) - Cent(2));
%         CpdDistToPoint2 = hypot( X(2) - Cent(1), Y(2) - Cent(2));
        Ind = find(DistToCpd < 0.8*CpdDistToPoint1);% & DistToCpd>CpdDistToPoint2);        
        pixx1T = pixx1(Ind); pixxT = pixx(Ind);        
        pixy1T = pixy1(Ind); pixyT = pixy(Ind);
        pixx1T = round(pixx1T); pixxT = round(pixxT);
        pixy1T = round(pixy1T); pixyT = round(pixyT);
        for i5 = 1:numel(pixxT)
            I(pixyT(i5):pixy1T(i5),pixxT(i5),1) = ColSpec(1);
            I(pixyT(i5),pixxT(i5):pixx1T(i5),1) = ColSpec(1); 
            I(pixyT(i5):pixy1T(i5),pixxT(i5),2) = ColSpec(2);
            I(pixyT(i5),pixxT(i5):pixx1T(i5),2) = ColSpec(2); 
            I(pixyT(i5):pixy1T(i5),pixxT(i5),3) = ColSpec(3);
            I(pixyT(i5),pixxT(i5):pixx1T(i5),3) = ColSpec(3); 
        end
    end
end