clc;clear; close all;
% read the ndpi file through the file; 
% Here users need to create their own file.
ndpipath= 'C:\Users\Qiang\Desktop\ndpi\';
dirfile= dir(fullfile('C:\Users\Qiang\Desktop\ndpi\*.ndpi'));
num = length(dirfile);
% create the file project file, H, E file;
% the following is abosulte path;
mkdir('C:\Users\Qiang\Desktop\project1','E');
mkdir('C:\Users\Qiang\Desktop\project1','H');
mkdir('C:\Users\Qiang\Desktop\project1','origin');
mkdir('C:\Users\Qiang\Desktop\project1','originlabel');
mkdir('C:\Users\Qiang\Desktop\project1','originlabelrgb');
Originstr0 = 'C:\Users\Qiang\Desktop\project1\origin\';
Originstr1 = 'C:\Users\Qiang\Desktop\project1\originlabel\';
Originstr2 = 'C:\Users\Qiang\Desktop\project1\originlabelrgb\';
Hstr0 = 'C:\Users\Qiang\Desktop\project1\H\';
Estr0 = 'C:\Users\Qiang\Desktop\project1\E\';

%file = imread([filepath dirfile(k).name]);%
for k=1:num   
    sequence=num2str(k);      
    ndpifile=[ndpipath,'p',sequence,'.ndpi'];
     %read ndpi file %imread('I','PixelRegion',{[begin end ],[begin end]});
    I1=imread(ndpifile,'PixelRegion',{[10003,11002],[10003,11002]});
    
    [Inorm1,H1,E1] = normalizeStaining(I1);
    %H1 and E1 are rgb;rgb to gray then gray to black; 
    Hgray=rgb2gray(H1);
    HBW=imbinarize(Hgray);
    Egray=rgb2gray(E1);
    EBW=imbinarize(Egray);
   
    originfile=[Originstr0,sequence,'.png'];
    imwrite(I1,originfile);
    
    Hfile=[Hstr0,sequence,'.png'];
    imwrite(HBW,Hfile); 
    
    Efile=[Estr0,sequence,'.png'];
    imwrite(EBW,Efile);
    
    original=lablepicture(HBW,EBW,150);
    %original(original==255)=0;
    labelRGB = label2rgb(original);

    labelpic_path=[Originstr1,sequence,'.png'];
    labelpicrgb_path=[Originstr2,sequence,'.png'];
    imwrite(original,labelpic_path);
    imwrite(labelRGB,labelpicrgb_path);
   
end
