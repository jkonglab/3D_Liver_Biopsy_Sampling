
function original = lablepicture(A,B,thershold)
B=uint8(B);
A=uint8(A);
% C is backgroud part
C=A&B;
%D is non background part
D=~C;
D=uint8(D);

A1=A.*D;
B1=B.*D;

A1=logical(A1)*thershold;
B1=logical(B1);
%1 125 255
%original=A1+B1+C*255;
original=A1+B1;
original=uint8(original);
%original(original==0)=1;

%imshow(original)
clear Hstr1 Hstr1;
end