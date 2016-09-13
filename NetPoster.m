function NetPoster(LineStrength,OverlapC,ColormapC,BackCol,CpdCol,inactiveCol)
% NetPoster
% is a subfunction of MapAdjuster for post-processing the
% customized metabolic maps created by NetDraw, and possesses a structure
% very similar to that of the NetDraw.
% 
% Inputs:
% LineStrength: User defined of reactions' line width.
% OverlapC: Check if user wants to remove (0) or retain (1) overlapping
% reactions on the created map.
% ColormapC: Check if user wants to remove (1) or apply (0) color mapping.
% BackCol: Background color.
% CpdCol: Compounds' color.
% inactiveCol: Inactive reactions' color.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016


PostData = getappdata(0,'PostData');
h1 = getappdata(0,'ParentAx');
H1 = getappdata(0,'ParentFig');
delete(h1); delete(H1);
% Create new canvas -------------------------------------------------------
H1 = figure;
h1 = axes;
set(H1,'position',get(0,'screensize'));
set(h1,'position',[0,0,1,1],'xtick',[],'ytick',[]);
set(H1,'toolbar','figure');
set(H1,'menubar','figure');
setappdata(0,'ParentAx',h1); setappdata(0,'ParentFig',H1);
% Set Background color ----------------------------------------------------
if isempty(BackCol) || numel(BackCol) == 1
    BackColr = [0.8,0.992,0.8];
else
    BackColr = BackCol;
end
I = ones(PostData.max_y-PostData.min_y+50,...
    PostData.max_x-PostData.min_x+50,3);
I(:,:,1) = I(:,:,1).*BackColr(1);
I(:,:,2) = I(:,:,2).*BackColr(2); 
I(:,:,3) = I(:,:,3).*BackColr(3);
% Set inactive rxns -------------------------------------------------------
Allx4Arrow = PostData.Allx4Arrow; Ally4Arrow = PostData.Ally4Arrow;
if isempty(inactiveCol) || numel(inactiveCol) == 1
    Flaginactive = [0 0 0];
else
    Flaginactive = inactiveCol;
end
for i1 = 1:numel(Allx4Arrow)
    I = pixFill(I,Allx4Arrow{i1},Ally4Arrow{i1},LineStrength,Flaginactive);
end
% Flux overlay ------------------------------------------------------------
field_intsect_overlay = PostData.field_intsect_overlay;
field_intsect = PostData.field_intsect;
Cm=jet(numel(field_intsect_overlay)); flx = PostData.flx;
flx_abs = abs(flx);
for k=PostData.k
    McmCollect = zeros(size(flx,1),3); 
    if ~ColormapC
        for m=1:size(flx,1)
            % Normalize colors and insertion into enzyme boxes
            flxMax=max(abs(flx(:)));
            flxMin=min(abs(flx(:)));
            if flxMax==flxMin % One row (i.e. for 1 rxn)
                flxMax=max(abs(flx(:)));
                flxMin=min(abs(flx(:)));
            end
            flxSl=(size(Cm,1)-1)./(flxMax-flxMin);
            if isinf(flxSl) % in case of flxMax == flxMin
                flxSl = 1e15.*(size(Cm,1)-1);
            end
            Mcm=ceil((flx_abs(m,k)-flxMax).*flxSl+size(Cm,1));

            I = pixFill(I,PostData.ArX{m},PostData.ArY{m},LineStrength,...
                [Cm(Mcm,1) Cm(Mcm,2) Cm(Mcm,3)]);
            McmCollect(m,:) = [Cm(Mcm,1) Cm(Mcm,2) Cm(Mcm,3)];               
            clear pixx pixy
        end
    else % No color for flux carrying rxns
        for m=1:size(flx,1)
            I = pixFill(I,PostData.ArX{m},PostData.ArY{m},LineStrength,...
                Flaginactive);
            McmCollect(m,:) = Flaginactive;            
        end
    end
    ArX = PostData.ArX; ArY = PostData.ArY; OverlapData = PostData.OverlapData;
    min_x = PostData.min_x; min_y = PostData.min_y; basemap = PostData.basemap;
    RefOverlapData = PostData.RefOverlapData;
    if OverlapC
        [I,ArXo,ArYo,ArrowChek,RxnId] = OverlapChecker(I,ArX,ArY,OverlapData,...
        flx(:,k),McmCollect,field_intsect_overlay,LineStrength);
        [ArrowX,ArrowY,Cntrs] = ArrowheadProcess(Allx4Arrow,Ally4Arrow,...
            min_x,min_y,field_intsect,basemap,RxnId,RefOverlapData,'bigg',field_intsect_overlay);
        I = MapTrimmer(I,ArrowX,ArrowY,Cntrs,Allx4Arrow,Ally4Arrow,LineStrength,BackColr);
        [ArrowX1,ArrowY1,Cntrs1,Director] = ArrowheadProcess(ArXo,ArYo,min_x,min_y,...
            field_intsect_overlay,basemap,[],[],'bigg',[]);
        I = MapTrimmerF(I,ArrowX1,ArrowY1,Cntrs1,ArXo,ArYo,RxnId,field_intsect_overlay,LineStrength,BackColr);
    else
        [ArrowX,ArrowY,Cntrs] = ArrowheadProcess(Allx4Arrow,Ally4Arrow,...
            min_x,min_y,field_intsect,basemap,[],[],'bigg',[]);
        I = MapTrimmer(I,ArrowX,ArrowY,Cntrs,Allx4Arrow,Ally4Arrow,LineStrength,BackColr);
        [ArrowX1,ArrowY1,~,Director] = ArrowheadProcess(ArX,ArY,min_x,min_y,...
            field_intsect_overlay,basemap,[],[],'bigg',[]);
    end
    % Compounds overlay ---------------------------------------------------
    if isempty(CpdCol) || numel(CpdCol) == 1
        CpdColr = [0,1,0.498];
    else
        CpdColr = CpdCol;
    end
    for i1 = 1:size(PostData.cpdxy,1)
        cc1 = PostData.x_cpd(i1)-PostData.cpdxy(i1,4)*0.7:PostData.x_cpd(i1)+PostData.cpdxy(i1,4)*0.7;
        cc2 = PostData.y_cpd(i1)-PostData.cpdxy(i1,4)*0.7:PostData.y_cpd(i1)+PostData.cpdxy(i1,4)*0.7;
        [x,y] = meshgrid(cc1,cc2);
        f1= ((x-PostData.x_cpd(i1)).^2+(y-PostData.y_cpd(i1)).^2) <=(PostData.cpdxy(i1,4)*0.7)^2 ;
        f2= ((x-PostData.x_cpd(i1)).^2+(y-PostData.y_cpd(i1)).^2) <=(PostData.cpdxy(i1,4)*0.55)^2 ;
        x1 = x(f1); y1= y(f1);x1=round(x1); y1=round(y1);
        x2 = x(f2); y2= y(f2); x2=round(x2); y2=round(y2);
        x1(~y1) = []; y1(~y1) = [];
        x2(~y2) = []; y2(~y2) = [];
        for i2 = 1:numel(x1)
            I(y1(i2),x1(i2),1) = 0.502;
            I(y1(i2),x1(i2),2) = 0.502;
            I(y1(i2),x1(i2),3) = 0;
        end
        for i2 = 1:numel(x2)
            I(y2(i2),x2(i2),1) = CpdColr(1);
            I(y2(i2),x2(i2),2) = CpdColr(2);
            I(y2(i2),x2(i2),3) = CpdColr(3);
        end
        clear f1 x1 x2 cc1 cc2
    end
    % Show map ------------------------------------------------------------
    imshow(I,'InitialMagnification', 'fit','Parent',h1); 
    % Show time-series values on the corresponding image
%     TimePropUDF ('One',PostData.TextFont,PostData.TimeUnit,PostData.mTime,k,h1) 
    arrowhead(ArrowX,ArrowY,[],Flaginactive,LineStrength/2,[]);
    flxDirect.D = Director; flxDirect.F = flx(:,k);
    if OverlapC
        arrowhead(ArrowX1,ArrowY1,ArrowChek,McmCollect,LineStrength/2,flxDirect);
    else
        arrowhead(ArrowX1,ArrowY1,[],McmCollect,LineStrength/2,flxDirect);
    end
    
    Mdl = getappdata(0,'Mdl');
    Uni = load(which('UniModelKEGG.mat'));
    K = Uni.B2Kegg.K; K1 = Mdl.B2Kegg.K;
    B = Uni.B2Kegg.B; B1 = Mdl.B2Kegg.B;
    clear Uni
    HoverUs.K = K1; HoverUs.B = B1;
    HoverUs.R = PostData.OverlayRxnsIDs; HoverUs.P = field_intsect_overlay';
    HoverUs.CR = PostData.ConOverlayRxnID; 
    HoverAll.R = PostData.rxnCds; HoverAll.P = field_intsect;
    HoverAll.K = K; HoverAll.B = B;
    pixHover(PostData.x_cpd,PostData.y_cpd,0.7*PostData.cpdxy(1,4),PostData.cpdnames...
        ,PostData.Allx,PostData.Ally,...
        HoverAll,HoverUs,abs(flx(:,k)),h1,H1)
    % Show colorbar--------------------------------------------------------
    if ~ColormapC
        flxstr=linspace(flxMin,flxMax,7);
        flxstr1=cell(numel(flxstr),1);
        for ir=1:numel(flxstr)
            flxstr1{ir,1}=num2str(flxstr(ir),2);
        end
        ch = colorbar('peer',h1);
        set(ch,'yticklabel',flxstr1)
    end
end
