%The output: ratio is liver fibrosis ratio in the biopsy sample
%             CSimg is tube cross section image is saved in the setting savepath
%The input:
%           projectpath is origin image file path；
%           x0,y0,z0 is penetrarion coordinate position
%           l is the tube length; r is tube radius
%           theta,phi: tube insert angle
%           CSR:  tube cross section rotation angle; cross-section planes can rotate in different rotation angles
%                  CSR=0 means cross section parallel with x-y plane;
%                  CSR=90 prependicular with x-y plane;
%           nc: number of cross-section images perpendicular to tube middle axis
% % tube size:
% l=1000;r=25;
% %penetration position
% x0=1000;y0=1000;z0=25;
% %penetration angle
% theta=90;phi=90;
% CSR = 0; 
% % load the image to build tissue voxel model
function [ratio,CSimg]=tissuemodel(projectpath,x0,y0,z0,l,r,theta,phi,CSR,nc,Scale)
if isempty(projectpath)
     error('The first inputs of projectpath are needed.')  %https://www.mathworks.com/matlabcentral/fileexchange/27056-set-default-values 
end
ndpipath=[projectpath,'imgsrc\'];
imgsrc =ndpipath;
imgsrc1 = [projectpath,'imgsrclabel\'];
imgsrc2 =[projectpath,'imgsrclabelrgb\']; 
mkdir(projectpath,'output');
savepath=[projectpath,'output\'];


%XX,YY,ZZ tissue length, width and height
% Get the 3-D size information from the slides
path =imgsrc;
imgname= [path,'1.jpg']; 
Isize=imread(imgname);
imgnamen=[path,'*.jpg'];
dirfile= dir(fullfile(imgnamen));
[YY,XX,~]=size(Isize);
ZZ = length(dirfile);
NUM=ZZ;


% set defalut value
switch nargin
     case 0
        x0= [];
        y0= [];
        z0= [];
     case 1
       l= [];
       r= [];
     case 2
       theta= [];
       phi= [];
     case 3      
       CSR= [];
     case 4
       nc= [];
     case 5
       imgsrc=[];
       imgsrc1=[];
       imgsrc2=[];
       savepath=[];
end
if isempty(x0)&&isempty(y0)&&isempty(z0)
    x0 =XX;
    y0=YY;
    z0=Scale*ZZ-10;
end
if isempty(l)&&isempty(r)
  l = XX;
  r = 40;
end
if isempty(theta)&&isempty(phi)
    theta=60;
    phi=60;
end
if isempty(CSR)
    CSR=0;
end
if isempty(nc)
    nc=2;
end



%XX,YY,ZZ tissue length, width and height
labelimgmodel(1:YY,1:XX,1:ZZ) = uint8(zeros(YY,XX,ZZ));
imgmodel(1:YY,1:XX,1:3,1:ZZ) = uint8(zeros(YY,XX,3,ZZ));

%RGB value of the image VR,VG,VB
VR(1:YY,1:XX,1:NUM) = double(zeros(YY,XX,NUM));
VG(1:YY,1:XX,1:NUM) = double(zeros(YY,XX,NUM));
VB(1:YY,1:XX,1:NUM) = double(zeros(YY,XX,NUM));

for k=1:ZZ  
% store the labeled pixels
labelpic=imgsrc1+string(k)+'.png';
labelpic = imread(labelpic);
labelpic=double(labelpic);
labelimgmodel(1:YY,1:XX,k)=labelpic;

Oname=imgsrc+string(k)+'.jpg';
Oname = imread(Oname);
Oname=double(Oname);

%load labeled origin image to build tissue
%Oname1=imgsrc2+string(k)+'.png';
%Oname1 = imread(Oname1);
%Oname1=double(Oname1);
imgmodel(1:YY,1:XX,1:3,k)=Oname;
% read original rgb value
VR(1:YY,1:XX,k) = Oname(:,:,1);
VG(1:YY,1:XX,k) = Oname(:,:,2);
VB(1:YY,1:XX,k) = Oname(:,:,3);
end
%aim for computing the ratio
%aim1 for interpolation for multiple layer
%aim2 for  rotation cross-section plane interpolation
%aim3 for  vertical cross-section plane interpolation 
[aim,aim1,aim2,aim3]=tubevoxel(XX,YY,ZZ,x0,y0,z0,l,r,theta,phi,CSR,nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%display tissue model if running in the hpc this part can be deleted
num1 =[1:1:ZZ];
k1=1;
[x,y] = meshgrid(1:XX,1:YY); %!!!!y,x rotate
for k2=1:ZZ
  z = k1*ones(YY,XX);
  b = surf(x,y,z,uint8(imgmodel(:,:,1:3,num1(k2))));
  grid off
  set(b,'linestyle','none');
  view(3);hold on;%title('3-D tissue model');
  k1=k1+Scale;
end
 daspect([10,10,10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %color interpolation for multiple layer 
% for ccs=0:(ZZ-1)
% aa0=2*r+1;
% aa1=aa0*(l+1);
% aim1x=reshape(aim1((aa1*ccs+1):(aa1*(ccs+1)),1),aa0,[]);
% aim1y=reshape(aim1((aa1*ccs+1):(aa1*(ccs+1)),2),aa0,[]);
% aim1z=reshape(aim1((aa1*ccs+1):(aa1*(ccs+1)),3),aa0,[]);
% %aim1z = aim1z/Scale; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% %interpolation method   
% VRq = interp3(VR,aim1x,aim1y,aim1z,'makima',255);
% VGq = interp3(VG,aim1x,aim1y,aim1z,'makima',255);
% VBq = interp3(VB,aim1x,aim1y,aim1z,'makima',255);
% % R,G,B channels
% CSimg(:,:,1)=VRq;
% CSimg(:,:,2)=VGq;
% CSimg(:,:,3)=VBq;
% % saved the image
% CSimg=uint8(CSimg);
% sequence=num2str(ccs+1); 
% CSimgname=[savepath,'CSimg',sequence,'.png'];
% imwrite(CSimg,CSimgname);
% clearvars CSimg;
% end
% 
% % central cross-section color interpolation  
% aim2x=reshape(aim2(:,1),2*r+1,[]);
% aim2y=reshape(aim2(:,2),2*r+1,[]);
% aim2z=reshape(aim2(:,3),2*r+1,[]);
% %interpolation method   
% VRq2 = interp3(VR,aim2x,aim2y,aim2z,'makima',255);
% VGq2 = interp3(VG,aim2x,aim2y,aim2z,'makima',255);
% VBq2 = interp3(VB,aim2x,aim2y,aim2z,'makima',255);
% CSMimg(:,:,1)=VRq2;
% CSMimg(:,:,2)=VGq2;
% CSMimg(:,:,3)=VBq2;
% %save image
% CSMimg=uint8(CSMimg);
% CSimgname=[savepath,'CSMimg','.png'];
% imwrite(CSMimg,CSimgname);
% clearvars CSMimg;
% % figure;
% % imshow(CSimg1); %Uint8 can show the cross-section
% % hold on;title('Tube cross-section');
% 
% %color interpolation of cross-section images perpendicular to tube middle axis
% for vi=0:(nc-1)
%     vaa=2*r+1;
%     vaa2=vaa*(ZZ+1); 
% aim3x=reshape(aim3(vaa2*vi+1:(vi+1)*vaa2,1),vaa,[]); %(vi*aa+1):(2*vi*aa+1)
% aim3y=reshape(aim3(vaa2*vi+1:(vi+1)*vaa2,2),vaa,[]);
% aim3z=reshape(aim3(vaa2*vi+1:(vi+1)*vaa2,3),vaa,[]);
% VRq3 = interp3(VR,aim3x,aim3y,aim3z,'makima',255);
% VGq3 = interp3(VG,aim3x,aim3y,aim3z,'makima',255);
% VBq3 = interp3(VB,aim3x,aim3y,aim3z,'makima',255);
% VCSimg(:,:,1)=VRq3;
% VCSimg(:,:,2)=VGq3;
% VCSimg(:,:,3)=VBq3;
% 
% % saved the image
% VCSimg=uint8(VCSimg);
% sequence1=num2str(vi+1); 
% VCSimgname=[savepath,'VCSimg',sequence1,'.png'];
% imwrite(VCSimg,VCSimgname);
% clearvars VCSimg;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% [nn,~]=size(aim);
% for jj1=1:nn
% O(jj1)=labelimgmodel(aim(jj1,2),aim(jj1,1),aim(jj1,3));
% end
% 
% NumH0= sum(sum(O==1));% H1 is blue voxel   26   51  77
% NumE0= sum(sum(O==2));% E2 is red voxel
% NumB0= sum(sum(O==0));  %0
% 
% ratio1=NumH0/(NumH0+NumE0);
% ratio2=NumE0/(NumH0+NumE0);
% ratio3=NumH0/(NumH0+NumE0+NumB0);
% ratio4=NumE0/(NumH0+NumE0+NumB0);
% ratio5=NumB0/(NumH0+NumE0+NumB0);
% ratio=[ratio1,ratio2,ratio3,ratio4,ratio5];
% fprintf('\n The liver fibrosis ratio blue is %0.4f \n', ratio1);
% fprintf('\n The liver fibrosis ratio red is %0.4f \n', ratio2);
% fprintf('\n The liver fibrosis ratio blue (include background) is %0.4f \n', ratio3);
% fprintf('\n The liver fibrosis ratio red (include background) is %0.4f \n', ratio4);
% fprintf('\n The liver fibrosis ratio background (include background) is %0.4f \n', ratio5);
% 
% formatSpec = 'ratio is %4.2f \n';
% A=fprintf(formatSpec,ratio);
% dlmwritename=[savepath,'ratio.txt'];
% dlmwrite(dlmwritename,ratio);
end