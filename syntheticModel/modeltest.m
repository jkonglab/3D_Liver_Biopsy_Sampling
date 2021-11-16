clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build blocks tissue size 
XX=300;
YY=300;
ZZ=600;
%Tube size:l is length; d is diameter 
l=600;d=50;
%and penetraion position and angle
x0=200;y0=100;z0=600;
theta=0;phi=0;

crossectR = 60; % crossection rotation angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image=zeros(XX,YY,3); %initialize
image(:,1:100,1)=0.5;   %Red (dark red)
image(1:150,1:100,3)=1;
image(:,101:200,1)=1;   %Red (maximum value)
image(1:150,:,2)=1;%rand(100,300);  %Green
image(150:XX,200:YY,3)=0.5;

image2=zeros(XX,YY,3); %initialize
image2(:,1:100,1)=0.3;   %Red (dark red)
image2(1:150,1:100,3)=2;
image2(:,101:200,1)=1;   %Red (maximum value)
image2(1:150,:,2)=0.5;%rand(100,300);  %Green
image2(150:XX,200:YY,3)=1;
rgbimage2 = imrotate(image2,270);
rgbimage = imrotate(image,270);


img1label=rgb2gray(rgbimage);
img1label=img1label*10000;
img1label=round(img1label);
img1label(img1label==2989)=3000;
img1label2=rgb2gray(rgbimage2);
img1label2=img1label2*10000;
img1label2=round(img1label2);

NUM=ZZ;
orgbpicgather(1:XX,1:YY,1:3,1:NUM) = double(zeros(XX,YY,3,NUM));
opicgather(1:XX,1:YY,1:NUM) = double(zeros(XX,YY,NUM));


%RGB value of the pictures VR,VG,VB
VR(1:XX,1:YY,1:NUM) = double(zeros(XX,YY,NUM));
VG(1:XX,1:YY,1:NUM) = double(zeros(XX,YY,NUM));
VB(1:XX,1:YY,1:NUM) = double(zeros(XX,YY,NUM));
% build the model
NN1=NUM/2;
 for k=1:NN1  %
 orgbpicgather(1:XX,1:YY,1:3,k) = image;%
 opicgather(1:XX,1:YY,k)=img1label;
 
 VR(1:XX,1:YY,k) = image(:,:,1);% load the RGB value
 VG(1:XX,1:YY,k) = image(:,:,2);
 VB(1:XX,1:YY,k) = image(:,:,3);
 end
  for k1=(NN1+1):600  %
orgbpicgather(1:XX,1:YY,1:3,k1) =image2;%
opicgather(1:XX,1:YY,k1)=img1label2;

VR(1:XX,1:YY,k1) = image2(:,:,1);%
VG(1:XX,1:YY,k1) = image2(:,:,2);
VB(1:XX,1:YY,k1) = image2(:,:,3);
  end
% display the model
num1 =[1:1:NUM];
k1=0;
[x,y] = meshgrid(1:XX,1:YY);%!!!!y,x rotate
for k2=1:NUM
    z = k1*ones(XX,YY);
    b = surf(x,y,z,double(orgbpicgather(:,:,1:3,num1(k2)))); 
%  colormap(gray);
  grid off
  set(b,'linestyle','none');
  view(3);hold on;
  k1=k1+1;
end
 daspect([10,10,10]);
%aim is tube discrete point, aim1 is discrete point in the tube cross-section 
[aim,aim1]=tubevoxel(XX,YY,ZZ,x0,y0,z0,l,theta,phi,d,crossectR);%tubevoxel619
[nn,mm]=size(aim);
aim1x=reshape(aim1(:,1),d-1,[]);
aim1y=reshape(aim1(:,2),d-1,[]);
aim1z=reshape(aim1(:,3),d-1,[]);
% image interpolation throgh R,G,B channels
VRq = interp3(VR,aim1x,aim1y,aim1z,'linear',255);
VGq = interp3(VG,aim1x,aim1y,aim1z,'linear',255);
VBq = interp3(VB,aim1x,aim1y,aim1z,'linear',255);
origion3(:,:,1)=VRq;
origion3(:,:,2)=VGq;
origion3(:,:,3)=VBq;
figure
imshow(origion3);

for jj1=1:nn
O(jj1)=opicgather(aim(jj1,1),aim(jj1,2),aim(jj1,3));
end

Num11_1=sum(sum(O==8505));% 8505
Num11_2=sum(sum(O==8860));% 8860
Num11_3=sum(sum(O==5870));% 5870
Num12_1=sum(sum(O==1495));% 1495
Num12_2=sum(sum(O==3000));% 3000
Num12_3=sum(sum(O==570));% 570

Num21_1=sum(sum(O==6112));% 6122
Num21_2=sum(sum(O==5925));% 
Num21_3=sum(sum(O==2935));% 5870
Num22_1=sum(sum(O==897));% 897
Num22_2=sum(sum(O==2989));% 2989
Num22_3=sum(sum(O==1140));% 1140
overall=Num21_1+Num21_2+Num21_3+Num22_1+Num22_2+Num22_3+Num11_1+Num11_2+Num11_3+Num12_1+Num12_2+Num12_3;
ratio11_1=Num12_1/overall;
ratio11_2=Num12_2/overall;
ratio11_3=Num12_3/overall;
ratio12_1=Num11_1/overall;
ratio12_2=Num11_2/overall;
ratio12_3=Num11_3/overall;

% ratio21_1=Num21_1/overall;
% ratio21_2=Num21_2/overall;
% ratio21_3=Num21_3/overall;
% ratio22_1=Num22_1/overall;
% ratio22_2=Num22_2/overall;
% ratio22_3=Num22_3/overall;
ratio21_1=Num22_1/overall;                      
ratio21_2=Num22_2/overall;                 
ratio21_3= Num22_3/overall;
ratio22_1= Num21_1/overall;
ratio22_2= Num21_2/overall;
ratio22_3= Num21_3/overall;


figure
imshow(image2);
figure
imshow(image);
ratio=[ratio11_1,ratio11_2,ratio11_3,ratio12_1,ratio12_2,ratio12_3;
    ratio21_1,ratio21_2,ratio21_3,ratio22_1,ratio22_2,ratio22_3]
