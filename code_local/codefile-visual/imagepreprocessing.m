clc;clear; close all;
% read the ndpi file through the file; 
% Here users need to create their own file.
projectpath='..\ProjectImgFile\';
ndpipath=[projectpath,'imgsrc\']; 
dirname = [ndpipath,'*.jpg'];
dirfile= dir(fullfile(dirname));
num = length(dirfile);
% create the file project file, H, E file;
% the following is abosulte path;
mkdir(projectpath,'E');
mkdir(projectpath,'H');
%mkdir(projectpath,'origin');
mkdir(projectpath,'imgsrclabel');
mkdir(projectpath,'imgsrclabelrgb');
imgsrc = [projectpath,'imgsrc\']; 
imgsrc1 = [projectpath,'imgsrclabel\'];
imgsrc2 = [projectpath,'imgsrclabelrgb\'];
Hstr0 = [projectpath,'H\']; 
Estr0 = [projectpath,'E\']; 

%file = imread([filepath dirfile(k).name]);%
for k=1:num   
    sequence=num2str(k);      
    %ndpifile=[ndpipath,'p',sequence,'.jpg'];
    ndpifile=[ndpipath,sequence,'.jpg'];
     %read ndpi file %imread('I','PixelRegion',{[begin end ],[begin end]});
    I1=imread(ndpifile);
     %I1=imread(ndpifile,'PixelRegion',{[10003,11002],[10003,11002]});
    
    [Inorm1,H1,E1] = normalizeStaining(I1);
    %H1 and E1 are rgb;rgb to gray then gray to black; 
    Hgray=rgb2gray(H1);
    HBW=imbinarize(Hgray);
    Egray=rgb2gray(E1);
    EBW=imbinarize(Egray);
   
    %originfile=[Originstr0,sequence,'.png'];
    %imwrite(I1,originfile);
    
    Hfile=[Hstr0,sequence,'.png'];
    imwrite(HBW,Hfile); 
    
    Efile=[Estr0,sequence,'.png'];
    imwrite(EBW,Efile);
    
    original=lableimg(HBW,EBW,2);
    cmap=[0 0 1; 1 0 0];
    labelRGB =label2rgb(original,cmap,[1 1 1]);

    labelpic_path=[imgsrc1,sequence,'.png'];
    labelpicrgb_path=[imgsrc2,sequence,'.png'];
    imwrite(original,labelpic_path);
    imwrite(labelRGB,labelpicrgb_path);
   
end


