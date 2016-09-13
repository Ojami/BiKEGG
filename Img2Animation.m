function Img2Animation (Idf)
% Img2Animation
% gets consecutive image files, and based on Idf value
% generates a Video or GIF file as output. Frames are sorted based on image
% numbers, which performed by SaveFlxImgs function
%
% O. Jamialahmadi
% TMU, Chem. Eng. Dept., Biotech. Group 
% Sep. 2015

if isempty(Idf) || Idf > 2 || Idf <= 0
    waitfor(msgbox('Unknown identifier (1 | 2)','Error','error'));
    return
end
% Check if image formats are consistent
ImagesFolder=uigetdir;
TempF=dir(ImagesFolder);
TempF1=TempF(3).name;
T1=strfind(TempF1,'.');
TempF2=TempF1(T1+1:end);
if size(TempF,1) < 4
     msgbox(['There is only one image ',...
     'in this directory!'],'Warning','warn');
     return
end
TempG1=TempF(4).name;
G1=strfind(TempF1,'.');
TempG2=TempG1(G1+1:end);
if ~strcmp(TempF2,TempG2)
     msgbox(['File types in the chosen ',...
         'folder are not consistent!'],'Warning','warn');
     return
end

ImgFiles = dir(strcat(ImagesFolder,'\*',TempF2));

S = [ImgFiles(:).datenum]; 
[~,S] = sort(S);
ImgFilesS = ImgFiles(S);

if Idf == 1 % Video file
    VideoFile=strcat(ImagesFolder,'\Output video');
    writerObj = VideoWriter(VideoFile);
%     Input dialog box for framerate
    prompt='Insert the desired framerate value:';
    name='Frame rate: ';
    numlines=1;
    defaultanswer={'1'};
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    fps=str2double(answer(1)); 
    writerObj.FrameRate = fps;
    open(writerObj);
    for t= 1:length(ImgFilesS)
         Frame=imread(strcat(ImagesFolder,'\',ImgFilesS(t).name));
         writeVideo(writerObj,im2frame(Frame));
    end
    close(writerObj);
elseif Idf == 2 % GIF file
    GIFfilename = strcat(ImagesFolder,'\OutputAnimatedImg.gif');
    for t= 1:length(ImgFilesS)
        Frame=imread(strcat(ImagesFolder,'\',ImgFilesS(t).name));
        [imind,cm] = rgb2ind(Frame,256);
        if t == 1;
            imwrite(imind,cm,GIFfilename,'gif','LoopCount',Inf,'DelayTime',1);
        else
           imwrite(imind,cm,GIFfilename,'gif','WriteMode','append','DelayTime',1);
        end
    end
end