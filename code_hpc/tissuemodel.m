%The output: ratio is liver fibrosis ratio in the biopsy sample
%             CSimg is tube cross section image is saved in the setting savepath
%The input:
%           projectpath is origin image file path；
%           x0,y0,z0 is penetrarion coordinate position
%           l is the tube length; r is tube radius
%           theta,phi: tube insert angle
% clc;clear all;
% l=1650; r=160; x0=1770; y0=973;  z0=110; theta=89.9;phi=0;

function[ratio]=tissuemodel(projectpath,probeset,x0,y0,z0,l,r,theta,phi) %,CSR
if isempty(projectpath)
     error('The first inputs of projectpath are needed.')  %https://www.mathworks.com/matlabcentral/fileexchange/27056-set-default-values 
end
ndpipath=[projectpath,'imgsrc/']; 
imgsrc =ndpipath;
imgsrc1 = [projectpath,'imgsrclabel/'];

% write the image into the output file
outputname=strcat('output_',probeset);
mkdir(projectpath,outputname); 
savepath=strcat(projectpath,outputname,'/');

%XX,YY,ZZ tissue length, width and height
% Get the 3-D size information from the slides
path =imgsrc1;
imgname= [path,'1.png']; 
Isize=imread(imgname);
imgnamen=[path,'*.png'];
dirfile= dir(fullfile(imgnamen));
[YY,XX,~]=size(Isize); clearvars Isize; 
ZZ = length(dirfile);
NUM=ZZ;

ScaleC = 5/0.44371478013932641;

%XX,YY,ZZ tissue length, width and height
labelimgmodel(1:YY,1:XX,1:ZZ) = uint8(zeros(YY,XX,ZZ));

%RGB value of the image VR,VG,VB
VR(1:YY,1:XX,1:NUM) = double(zeros(YY,XX,NUM));
VG(1:YY,1:XX,1:NUM) = double(zeros(YY,XX,NUM));
VB(1:YY,1:XX,1:NUM) = double(zeros(YY,XX,NUM));

fprintf('\n load image data \n');
for k=1:ZZ  
% store the labeled pixels
labelimg=imgsrc1+string(k)+'.png';
labelimg = imread(labelimg);
labelimg=double(labelimg);
labelimgmodel(1:YY,1:XX,k)=labelimg;
clearvars labelimg; 
img=imgsrc+string(k)+'.jpg';
img = imread(img);
img=double(img);

% read original rgb value
VR(1:YY,1:XX,k) = img(:,:,1); VR(1,1,1)=255;
VG(1:YY,1:XX,k) = img(:,:,2); VG(1,1,1)=255;
VB(1:YY,1:XX,k) = img(:,:,3); VB(1,1,1)=255;

clearvars img; 
end
fprintf('\n image data loaded  \n');
%aim1 is coordinates for interpolation for multiple layer

%interpolate 3D data volume dx=dy=0.44371478013932641, dz=5
% [X,Y,Z] = meshgrid(1:XX, 1:YY, [ScaleC:ScaleC:(ZZ)*ScaleC]);
x1=x0-l*sind(theta)*cosd(phi);
y1=y0-l*sind(theta)*sind(phi);
z1=z0-l*cosd(theta);

fprintf('\n cross-section images perpendicular to tube middle axis  \n');
% cross-section images perpendicular to tube middle axis

x=0:1:l; % l represents probe length
y=-r:1:r; % r represents radius
small = min(ZZ*ScaleC, 2*r); 
zstart=-round(small/2);  
zend =round(small/2);
%NUM_PLANE = length(z);
% initial the number of voxel: Total_Background=0; Total_Red=0;Total_Blue=0;
Total_Background=0; Total_Red=0;Total_Blue=0;
%for CSR=(0.1:-0.02:-0.1)
%for CSR=(0.25:-0.05:-0.25)
for CSR=(2.5:-0.5:-2.5)
for zz1=0:0
%for zz1=zstart:(zend) 
%generate the each layers discrete points coordinates
[xx,yy,zz] = meshgrid(x,y,zz1);
R=sqrt(yy.^2+zz.^2);
t=true(size(xx,1),size(xx,2),size(xx,3)); % t for the tube internal points index
t(R<=r)=true; % obtain the tube point index t
aim(:,1)= xx(:); clearvars xx;
aim(:,2)= yy(:); clearvars yy;
aim(:,3)= zz(:); clearvars zz;
%-----------------------------------------------------------------------
% plane rotation when CSR=0 it paralle with the x-y plane initially; 
% CSR=90 means it prependicular to x-y plane initially;

aim(:,2) = aim(:,2)*cosd(CSR);
aim(:,3) = aim(:,2)*sind(CSR);
%-----------------------------------------------------------------------


htheta=90-theta;
htheta=-htheta;
hphi=-phi;

Rz=[cosd(hphi),sind(hphi),0;-sind(hphi),cosd(hphi),0;0,0,1];
Ry=[cosd(htheta),0,sind(htheta);0,1,0;-sind(htheta),0,cosd(htheta)];
% rotation
aim=Ry*aim';
aim=Rz*(aim);
clearvars hphi htheta;
p2=[x1,y1,z1];

siz=size(aim');
Vp2=repmat(p2,siz(1),1);
aim=aim'+Vp2;

%aim1(:,3)=aim(:,3)/ScaleC;
clearvars Vp2 Vp12 Vp22 siz siz1 siz2;
aim1x=reshape(aim(:,1), 2*r+1, l+1, []);
aim1y=reshape(aim(:,2), 2*r+1, l+1, []);
aim1z=reshape(aim(:,3), 2*r+1, l+1, []);

%creat the image file space: height:2*r+1,length:l+1; 
CSimg = uint8(zeros((2*r+1),(l+1),3));
    for row1=1:(2*r+1)
    %https://en.wikipedia.org/wiki/Trilinear_interpolation
    % aim1xmin indicates the lattice point below point x and aim1xmax,
    % similarly for y and z
        aim1xmin=floor(aim1x(row1,:));
        aim1xmax=ceil(aim1x(row1,:));
    %On a periodic and cubic lattice, let xd be the differences between each
    %of x,y, z and the smaller coordinate related
        xd=aim1x(row1,:)-aim1xmin;
        ixd=1-xd;

        aim1ymin=floor(aim1y(row1,:));
        aim1ymax=ceil(aim1y(row1,:));
        yd=aim1y(row1,:)-aim1ymin;
        iyd=1-yd; 

        Zaim1z=aim1z(row1,:)/ScaleC;
        aim1zmin=floor(Zaim1z(1,:)); %for qury pixel on z axis below points
        aim1zmax=ceil(Zaim1z(1,:)); %for qury pixel on z axis above points

        Z_min=aim1zmin*ScaleC;  %pixel position along the z below the points
        zd=(aim1z(row1,:)-Z_min)/ScaleC; %distance of pixel position along the z below the points
        izd=1-zd; % Distance of pixel position along the z above the points


    %construct the eight points around the query points
    %eight points: C000,C100,C101,C001,C011,C111,C110,C010
        C000=[aim1ymax;aim1xmin;aim1zmin];%C000(1,:)==aim1ymax;
        C000( :,C000(1,:)<1) = 1; C000(:, C000(1,:)>YY) =1;
        C000( :,C000(2,:)<1) = 1; C000(:, C000(2,:)>XX) = 1;
        C000( :,C000(3,:)<1) = 1; C000(:,C000(3,:)>ZZ) = 1;
        C000ind = sub2ind(size(VR), C000(1,:), C000(2,:), C000(3,:)); clearvars C000;
        VRC000=VR(C000ind );
        VGC000=VG(C000ind );
        VBC000=VB(C000ind ); clearvars C000ind;

        C100=[aim1ymax;aim1xmax;aim1zmin];
        C100( :,C100(1,:)<1) = 1; C100(:,C100(1,:)>YY) =1;
        C100( :,C100(2,:)<1) = 1; C100(:,C100(2,:)>XX) = 1;
        C100( :,C100(3,:)<1) = 1; C100(:,C100(3,:)>ZZ) = 1;
        C100ind = sub2ind(size(VR), C100(1,:), C100(2,:), C100(3,:)); clearvars C100;
        VRC100=VR(C100ind);
        VGC100=VG(C100ind);
        VBC100=VB(C100ind); clearvars C100ind;

        VRC00= VRC000.*ixd+VRC100.*xd;  clearvars VRC000 VRC100;
        VGC00= VGC000.*ixd+VGC100.*xd;  clearvars VGC000 VGC100;
        VBC00= VBC000.*ixd+VBC100.*xd;  clearvars VBC000 VBC100;

        C010=[aim1ymin;aim1xmin;aim1zmin];
        C010( :,C010(1,:)<1) = 1; C010(:, C010(1,:)>YY) =1;
        C010( :,C010(2,:)<1) = 1; C010(:, C010(2,:)>XX) =1;
        C010( :,C010(3,:)<1) = 1; C010(:, C010(3,:)>ZZ) =1;
        C010ind = sub2ind(size(VR), C010(1,:), C010(2,:), C010(3,:)); clearvars C010;
        VRC010=VR(C010ind); VGC010=VG(C010ind); VBC010=VB(C010ind); clearvars C010ind;

        C110=[aim1ymin;aim1xmax;aim1zmin];
        C110( :,C110(1,:)<1) = 1; C110(:, C110(1,:)>YY) =1;
        C110( :,C110(2,:)<1) = 1; C110(:, C110(2,:)>XX) = 1;
        C110( :,C110(3,:)<1) = 1; C110(:,C110(3,:)>ZZ) = 1;
        C110ind = sub2ind(size(VR), C110(1,:), C110(2,:), C110(3,:)); clearvars C110;
        VRC110=VR(C110ind); VGC110=VR(C110ind); VBC110=VR(C110ind); clearvars C110ind;
        VRC10= VRC010.*ixd+VRC110.*xd; clearvars VRC010 VRC110;
        VGC10= VGC010.*ixd+VGC110.*xd; clearvars VGC010 VGC110; 
        VBC10= VBC010.*ixd+VBC110.*xd; clearvars VBC010 VBC110;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VRC0= VRC00.*iyd+VRC10.*yd;  clearvars VRC10 VRC00;
        VGC0= VGC00.*iyd+VGC10.*yd;  clearvars VGC10 VGC00;
        VBC0= VBC00.*iyd+VBC10.*yd;  clearvars VBC10 VBC00;

        C001=[aim1ymax;aim1xmin;aim1zmax];%C000(1,:)=aim1ymax;
        C001( :,C001(1,:)<1) = 1; C001(:, C001(1,:)>YY) =1;
        C001( :,C001(2,:)<1) = 1; C001(:, C001(2,:)>XX) = 1;
        C001( :,C001(3,:)<1) = 1; C001(:,C001(3,:)>ZZ) = 1;
        C001ind = sub2ind(size(VR), C001(1,:), C001(2,:), C001(3,:)); clearvars C001;
        VRC001=VR(C001ind); VGC001=VG(C001ind); VBC001=VB(C001ind);  clearvars C001ind;
        C101=[aim1ymax;aim1xmax;aim1zmax];
        C101( :,C101(1,:)<1) = 1; C101(:, C101(1,:)>YY) =1;
        C101( :,C101(2,:)<1) = 1; C101(:, C101(2,:)>XX) = 1;
        C101( :,C101(3,:)<1) = 1; C101(:, C101(3,:)>ZZ) = 1;
        C101ind = sub2ind(size(VR), C101(1,:), C101(2,:), C101(3,:)); clearvars C101;
        VRC101=VR(C101ind); VGC101=VG(C101ind); VBC101=VB(C101ind);  clearvars C101ind;

        VRC01= VRC001.*ixd+VRC101.*xd; clearvars VRC001 VRC101;
        VGC01= VGC001.*ixd+VGC101.*xd; clearvars VGC001 VGC101;
        VBC01= VBC001.*ixd+VBC101.*xd; clearvars VBC001 VBC101;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        C011=[aim1ymin;aim1xmin;aim1zmax];
        C011( :,C011(1,:)<1) = 1; C011(:, C011(1,:)>YY) =1;
        C011( :,C011(2,:)<1) = 1; C011(:, C011(2,:)>XX) = 1;
        C011( :,C011(3,:)<1) = 1; C011(:, C011(3,:)>ZZ) = 1;
        C011ind = sub2ind(size(VR), C011(1,:), C011(2,:), C011(3,:)); 
        VRC011=VR(C011ind); VGC011=VG(C011ind); VBC011=VB(C011ind); clearvars C011ind;

        C111=[aim1ymin;aim1xmax;aim1zmax];
        C111( :,C111(1,:)<1) = 1; C111(:, C111(1,:)>YY) =1;
        C111( :,C111(2,:)<1) = 1; C111(:, C111(2,:)>XX) = 1;
        C111( :,C111(3,:)<1) = 1; C111(:, C111(3,:)>ZZ) =1;
        C111ind = sub2ind(size(VR), C111(1,:), C111(2,:), C111(3,:));
        VRC111=VR(C111ind);  VGC111=VG(C111ind); VBC111=VB(C111ind);  clearvars C111ind;
        VRC11= VRC011.*ixd+VRC111.*xd;  clearvars VRC011 VRC111;
        VGC11= VGC011.*ixd+VGC111.*xd;  clearvars VGC011 VGC111;
        VBC11= VBC011.*ixd+VBC111.*xd;  clearvars VBC011 VBC111;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        VRC1= VRC01.*iyd+VRC11.*yd;  clearvars VRC11 VRC01; 
        VGC1= VGC01.*iyd+VGC11.*yd;  clearvars VGC11 VGC01;
        VBC1= VBC01.*iyd+VBC11.*yd;  clearvars VBC11 VBC01;
        VRC=VRC0.*izd +VRC1.*zd; clearvars VRC0 VRC1;
        VGC=VGC0.*izd +VGC1.*zd; clearvars VGC0 VGC1;
        VBC=VBC0.*izd +VBC1.*zd; clearvars VBC0 VBC1;

        CSimg(row1,:,1)=VRC;  clearvars VRC;
        CSimg(row1,:,2)=VGC;  clearvars VGC;
        CSimg(row1,:,3)=VBC;  clearvars VBC;
        clearvars ixd xd yd iyd zd izd;

    end
CSimg=uint8(CSimg);
namelayer=zz1+57;
sequence=num2str(namelayer);
STRCSR=num2str(CSR);
CSimgname=[savepath,'CSimg',sequence,STRCSR,'.jpg'];%,CSR
imwrite(CSimg,CSimgname);
clearvars CSimg;


% %t is the discrete points in the tube
% coord = round([aim1y(t), aim1x(t), aim1z(t)/ScaleC]); % t indices for tube
% %clearvars aim1x aim1y aim1z t;
% coord( coord(:,1)<1, :) = []; coord( coord(:,1)>YY, :) = [];
% coord( coord(:,2)<1, :) = []; coord( coord(:,2)>XX, :) = [];
% coord( coord(:,3)<1, :) = []; coord( coord(:,3)>ZZ, :) = [];
% ind = sub2ind(size(labelimgmodel), coord(:,1), coord(:,2), coord(:,3));
% 
% labels = labelimgmodel(ind);
% Blue_voxels= sum(labels(:)==1);% H1 is blue voxel 
% Total_Blue= Blue_voxels+ Total_Blue;  
% Red_voxels= sum(labels(:)==2);% E2 is red voxel   Total_Red
% Total_Red = Total_Red +Red_voxels;
% Background_voxels= sum(labels(:)==0);  %0  Background_voxels Total_Background
% Total_Background =Total_Background+ Background_voxels;
end

end
clearvars coord;
clearvars labelimgmodel ind;

clearvars X Y Z;
clearvars VR VG VB;

clearvars labels;
%================================
ratio1=Total_Blue/(Total_Blue+Total_Red);
ratio2=Total_Red/(Total_Blue+Total_Red);
ratio3=Total_Blue/(Total_Blue+Total_Red+Total_Background);
ratio4=Total_Red/(Total_Blue+Total_Red+Total_Background);
ratio5=Total_Background/(Total_Blue+Total_Red+Total_Background);
ratio=[ratio1,ratio2,ratio3,ratio4,ratio5];
fprintf('\n The liver fibrosis ratio blue is %0.4f \n', ratio1);
fprintf('\n The liver fibrosis ratio red is %0.4f \n', ratio2);
fprintf('\n The liver fibrosis ratio blue (include background) is %0.4f \n', ratio3);
fprintf('\n The liver fibrosis ratio red (include background) is %0.4f \n', ratio4);
fprintf('\n The liver fibrosis ratio background (include background) is %0.4f \n', ratio5);

dlmwritename=[savepath,'ratio.txt'];
dlmwrite(dlmwritename,ratio);






end