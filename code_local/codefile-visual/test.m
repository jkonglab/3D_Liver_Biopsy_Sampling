    clc;clear;close all;
%The output: ratio file contains liver fibrosis ratio in the biopsy sample
%             CSimg is tube cross section image;cross-section images
%             perpendicular to tube middle axis;
%The input parameter:
%           projectpath is origin image file path；
%           x0,y0,z0 is penetrarion coordinate position;
%           l is the tube length; r is tube radius;
%           scale is the adjacent distance between each WSI;
%           theta,phi: tube insert angle
%           CSR:  tube cross section rotation angle; cross-section planes can rotate in different rotation angles
%                  CSR=0 means cross section parallel with x-y plane;
%                  CSR=90 prependicular with x-y plane;
%           nc: number of cross-section images perpendicular to tube middle axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tube size:
l=540;r=40;
%penetration position
x0=400;y0=500;z0=100; Scale = 5/0.44371478013932641;
% load the image to build tissue voxel model
projectpath= '..\ProjectImgFile\';
%penetration angle
theta=80;phi=45;
% crossection rotation angle
CSR =0; nc=4;
%tissuemodel(projectpath,[],[],[],[],[],theta,phi,CSR,nc,Scale);%
tissuemodel(projectpath,x0,y0,z0,l,r,theta,phi,CSR,nc);
%tissuemodel([],x0,y0,z0,l,r,theta,phi,CSR,nc);
%https://ww2.mathworks.cn/help/matlab/ref/inputparser.html
% https://www.mathworks.com/matlabcentral/answers/314783-how-to-assign-default-values-to-function-inputs
alpha(0.4);

