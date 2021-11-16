% clc;clear all;
% XX=500;YY=600;ZZ=10;
% theta=90;phi=90;
% %penetration position
% x0=400;y0=300;z0=8;
% l=500;
% r=7;
% ZZ=10;

% aim for computing the ratio
%aim1 for interpolation for multiple layer
%aim2 for  rotation cross-section plane interpolation
%aim3 for  vertical cross-section plane interpolation
function [aim,aim1,aim2,aim3]=tubevoxel(XX,YY,ZZ,x0,y0,z0,l,r,theta,phi,CSR,nc) 
x1=x0-l*sind(theta)*cosd(phi);
y1=y0-l*sind(theta)*sind(phi);
z1=z0-l*cosd(theta);
% cross-section images perpendicular to tube middle axis
NC=[0:( l/nc):(l-1)];
x=0:1:l;
y=-r:1:r;
%z=-r:1:r;
z=-round(ZZ/2):1:round(ZZ/2); 
[x,y,z] = meshgrid(x,y,z);
% contain the Voxel in the tube
% R=sqrt(y.^2+z.^2);
t=true(size(x,1),size(x,2),size(x,3));
% t(R<=r)=true;
ind = find(t);
aim(:,1)=x((ind))+1;
aim(:,2)=y((ind))+1;
aim(:,3)=z((ind))+1;


%keep the x-y cross section z=0
id2 = aim(:,3)==0;
aim2= aim(id2,:);
aim2(:,1) = aim2(:,1);
aim2(:,2) = aim2(:,2)*cosd(CSR);
aim2(:,3) = aim2(:,2)*sind(CSR);

t3=false(size(x,1),size(x,2),size(x,3));
[~,nnc] =size(NC);
for nnnc = 1:nnc
    t3(x==NC(nnnc))=true;
end
ind3 = find(t3);
aim3(:,1)=x((ind3))+1;
aim3(:,2)=y((ind3))+1;
aim3(:,3)=z((ind3))+1;
aim3=sortrows(aim3);
clearvars x y z ind t 


htheta=90-theta;
htheta=-htheta;
phi=-phi;

Rz=[cosd(phi),sind(phi),0;-sind(phi),cosd(phi),0;0,0,1];
Ry=[cosd(htheta),0,sind(htheta);0,1,0;-sind(htheta),0,cosd(htheta)];
% rotation
aim=Ry*aim';
aim=Rz*(aim);
aim=aim';

aim2=Ry*aim2';
aim2=Rz*(aim2);
aim2=aim2';

aim3=Ry*aim3';
aim3=Rz*(aim3);
aim3=aim3';


%scatter3(aim(:,1),aim(:,2),aim(:,3));

aim=aim';
aim2=aim2';
aim3=aim3';
% add the displacement for the tube

p2=[x1,y1,z1];

siz=size(aim');
Vp2=repmat(p2,siz(1),1);
aim=aim'+Vp2;
aim1=aim; % aim1 for interpolation  aim for calculate percentage

siz1=size(aim2');
Vp12=repmat(p2,siz1(1),1);
aim2=aim2'+Vp12;

siz2=size(aim3');
Vp22=repmat(p2,siz2(1),1);
aim3=aim3'+Vp22;

clearvars Vp2 Vp12 Vp22 siz siz1 siz2


aim=round(aim);
id =aim(:,1)<=0 ;
aim(id,:)=[];
id= aim(:,1)>XX;
aim(id,:)=[];
id= aim(:,2)<=0;
aim(id,:)=[];
id= aim(:,2)>YY;
aim(id,:)=[];
id = aim(:,3)<=0;
aim(id,:)=[];
id= aim(:,3)>ZZ;
aim(id,:)=[];

scatter3(aim1(:,1),aim1(:,2),aim1(:,3) );
hold on;
scatter3(aim2(:,1),aim2(:,2),aim2(:,3) );
hold on;
scatter3(aim3(:,1),aim3(:,2),aim3(:,3) );
hold on;

end


