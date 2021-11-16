clc;clear;
%project path
projectpath= './ProjectImgFile/';
%projectpath= 'C:\Users\Qiang\Desktop\0822\';
probeset = input('Choose one probe : ','s');
%CSR=0;
%penetration angle
switch probeset
    case 'a'
        l=16500; r=1600; x0=17700; y0=9730; z0=110/2; theta=90;phi=0;%a
    case 'b'   
        l=20520; r=1600; x0=20900; y0=12810; z0=100/2; theta=90;phi=0;%b
    case 'c'
        l=24460;r=1600; x0=24820; y0=16430; z0=110/2; theta=90;phi=0;%c
    case 'd'
        l=27860;r=1600; x0=28195; y0=18975; z0=110/2; theta=90;phi=0;%d
    case 'e'
        l=22500;r=1600; x0=23245; y0=21493; z0=110/2; theta=90;phi=0;%e
    case 'f'
        l=21150;r=1600; x0=2140 ; y0=22818; z0=110/2; theta=90;phi=100;%f
    case 'g'
        l=19000;r=1600; x0=2140 ;  y0=22818; z0=110/2; theta=90;phi=135;%g  
    case 'h'
        l=23500;r=1600; x0=4720 ; y0=24553; z0=110/2; theta=90;phi=90;%h
    case 'i'
        l=18900;r=1600; x0=4720 ; y0=24553; z0=110/2; theta=90;phi=135;%i
    case 'j'
        l=22900;r=1600; x0=6490 ; y0=25108; z0=110/2; theta=90;phi=90;%j
    case 'k'
        l=21000;r=1600; x0=9540 ; y0=25678; z0=110/2; theta=90;phi=90;%k
    case 'l'
        l=21213;r=1600; x0=9540 ; y0=26100; z0=110/2; theta=90;phi=135;%l
    case 'm'
        l=18500;r=1600;x0=12440; y0=25600; z0=110/2; theta=90;phi=90;%m
    case 'n'
        l=20500;r=1600 ;x0=15400; y0=25900; z0=110/2; theta=90;phi=45;%n
    case 'o'
        l=17000;r=1600 ;x0=15400; y0=25900; z0=100/2; theta=90;phi=90; %o
    case 'p'
        l=23000;r=1600 ;x0=18085; y0=23578; z0=100/2; theta=90;phi=45; % p
    case 'q'
        l=15110;r=1600 ;x0=18910; y0=24440; z0=110/2; theta=90;phi=90; %q
    case 'r'
        l=25000;r=1600 ;x0=19960; y0=22828; z0=110/2; theta=90;phi=45; %r 
    case 's'
        l=28000;r=1600 ;x0=23110; y0=21568; z0=110/2; theta=90;phi=45; %s
end
%l=1000;r=100 ;x0=1000; y0=1000; z0=100; theta=90;phi=45;
tissuemodel(projectpath,probeset,x0,y0,z0,l,r,theta,phi); %,CSR
%the results of the liver fibrosis ratio and probe cross-section image are
%in the project output file
