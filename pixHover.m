function pixHover(x,y,rad,cnames,Allx,Ally,AllHover,HoverUs,flx,h1,H1)
% pixHover
% a subfunction of NetDraw for showing information about
% reactions/compounds upon hovering the cursor on the customized map 
% created by NetDraw
% 
% Inputs:
% x, y: Compounds' coordinates.
% Allx, Ally: Reactions' coordinates.
% AllHover: BiGG/KEGG identifiers for all reactions.
% HoverUs: BiGG/KEGG identifiers for flux carrying reactions.
% flx: Flux rates.
% h1: Handle to map's axes.
% H1: Handle to map's figure.

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% July 2016

 Thresh = 15; % Threshold identifying the hovering area
 axes(h1)
 set(H1,'WindowButtonMotionFcn', @hoverCallback);
 textHdl = text('Color', 'black', 'VerticalAlign',...
     'Bottom','Interpreter','none','Parent',gca);
 set(textHdl,'BackgroundColor',[173,255,47]./255)
 set(textHdl,'EdgeColor',[34,139,34]./255)
 set(textHdl,'Parent',h1)
 textHdl1 = text('Color', 'black', 'VerticalAlign',...
     'Bottom','Interpreter','none','Parent',gca);
 set(textHdl1,'BackgroundColor',[144,238,144]./255)
 set(textHdl1,'EdgeColor',[34,139,34]./255)
 set(textHdl1,'Parent',h1)

 function hoverCallback(src, evt)
    mousePoint = get(h1, 'CurrentPoint');
    mouseX = mousePoint(1,1);
    mouseY = mousePoint(1,2);
    distancesToMouse = hypot(x - mouseX, y - mouseY);
    [~, ind] = min(abs(distancesToMouse));
    if abs(mouseX - x(ind)) < Thresh && abs(mouseY - y(ind)) < Thresh
        set(textHdl, 'String', ['ID = ', cnames(ind)]);
        set(textHdl, 'Position', [x(ind) + Thresh, y(ind) + Thresh])
        axes(h1)
        viscircles(h1,[x(ind) y(ind)],rad, 'EdgeColor',[1,0.8667,0.1765]);           
    else
        set(textHdl, 'String', '')
        hg1 = findobj('type', 'line');
        set(hg1, 'Visible','off');
    end
    inFlag = 0; cum_str = [];
    for i1 =1:numel(Allx)
        SX = Allx{i1};
        SY = Ally{i1};
        distancesToMouse1 = hypot(SX - mouseX, SY - mouseY);
        [~, ind1] = min(abs(distancesToMouse1));
        if abs(mouseX - SX(ind1)) < Thresh/4 && abs(mouseY - SY(ind1)) < Thresh/4
            inFlag = [i1,ind1];
            cum_str = [cum_str,AllHover.P(i1)];
        end
    end
    if inFlag(1)
        [N1,~] = (ismember(HoverUs.P,cum_str));
        AllKEGG = AllHover.R(ismember(AllHover.P,cum_str));
        [~,NA] = intersect(AllHover.K,AllKEGG);
        AllBiGG = AllHover.B(NA);
        ThisKEGGH = HoverUs.R(N1); ThisKEGG = HoverUs.CR(N1);
        [~,NW] = intersect(HoverUs.K,ThisKEGGH);
        ThisBiGG = HoverUs.B(NW);
        if isempty(ThisBiGG); ThisBiGG = {'None'}; end;
        if isempty(AllBiGG); AllBiGG = {'None'}; end;
        if ~isempty(ThisKEGG) % Flux carrying rxns ////////////////////
            if numel(ThisKEGG) > 1
                ForThisFlux = strrep(cellstr(num2str(flx(N1))),' ','');
                if size(ForThisFlux,1) > size(ForThisFlux,2)
                    ForThisFlux = ForThisFlux';
                end
                if size(ThisKEGG,1) > size(ThisKEGG,2)
                    ThisKEGG = ThisKEGG';
                end
                if size(ThisBiGG,1) > size(ThisBiGG,2)
                    ThisBiGG = ThisBiGG';
                end
                set(textHdl1, 'String', {['KEGG: ', strjoin(ThisKEGG,',')],...
                    ['BiGG: ', strjoin(ThisBiGG,',')],...
                    ['Rate: ',strjoin(ForThisFlux,',')]});
            elseif numel(ThisKEGG) == 1 && numel(ThisBiGG)>1
                if size(ThisBiGG,1) > size(ThisBiGG,2)
                    ThisBiGG = ThisBiGG';
                end
                ForThisFlux = num2str(flx(N1));
                set(textHdl1, 'String', {strjoin(['KEGG: ', ThisKEGG]),...
                    ['BiGG: ', strjoin(ThisBiGG,',')],...
                    ['Rate: ',ForThisFlux]});
            else
                ForThisFlux = num2str(flx(N1));
                set(textHdl1, 'String', {strjoin(['KEGG: ', ThisKEGG]),...
                    ['BiGG: ', ThisBiGG{1}],...
                    ['Rate: ',ForThisFlux]});
            end
        else % Inactive rxns //////////////////////////////////////////
            if numel(AllKEGG) > 1
                if size(AllBiGG,1) > size(AllBiGG,2)
                    AllBiGG = AllBiGG';
                end
                set(textHdl1, 'String', {['KEGG: ', strjoin(AllKEGG,',')],...
                    ['BiGG: ', strjoin(AllBiGG,',')]});
            elseif numel(AllKEGG) == 1 && numel(AllBiGG)>1
                if size(AllBiGG,1) > size(AllBiGG,2)
                    AllBiGG = AllBiGG';
                end
                set(textHdl1, 'String', {['KEGG: ', AllKEGG{1}],...
                    ['BiGG: ', strjoin(AllBiGG,',')]});
            else
                set(textHdl1, 'String', {['KEGG: ', AllKEGG{1}],...
                    ['BiGG: ', AllBiGG{1}]});
            end
        end
        set(textHdl1, 'Position', [Allx{inFlag(1)}(inFlag(2)) + Thresh, Ally{inFlag(1)}(inFlag(2)) + Thresh])
    else
        set(textHdl1, 'String', '')    
    end

end
end