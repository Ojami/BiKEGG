function SaveFlxImgs (hax,FolderPath,k,ImgRes,ImgFrmt)
% SavFlxImgs
% employs following input to export created maps as image files:
% 1- hax- Handle to current axes
% 2- FolderPath- Target folder for images to be saved
% 3- k- A counter, used for image numbering
% 4- ImgRes- Output image resolution 
% 5- ImgFrmt- Output image format

% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Mar. 2015

FrmDat={'PNG';'JPEG';'BMP';'TIFF';'PDF'};
ResDat={'300','350','400','450','500','550','600','1200'};
FrmSp={'-dpng';'-djpeg';'-dbmp';'-dtiff';'-dpdf'};
if isempty(ImgRes) || isempty(ImgFrmt)
    ImgRes=1;
    ImgFrmt=1;
end
Res1=ResDat{ImgRes};
Frm1=FrmDat{ImgFrmt};
Res=['-r',Res1];
Frmt=FrmSp{ImgFrmt};
ImgNm=['Image ',num2str(k),'.',lower(Frm1)];
FileNm=fullfile(FolderPath,ImgNm);
set(gcf,'PaperPositionMode','auto')
print (hax,Res,Frmt,FileNm);