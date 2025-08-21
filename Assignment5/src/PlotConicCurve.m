% Function to plot a conic curve 'Omega' in an image 'img' painted 
% with color 'Paint' (3 dimensional vector with RGB values between 0 and
% 255)
% 
% [img]=PlotConicCurve(Omega,img,Paint);
% 
% Outputs:
% img - color image with thhe painted conic
%
% ---------
% Joao P. Barreto / University of Coimbra. 
% This software is exclusively meant for academic use. 
% Send your feedback to jpbar@deec.uc.pt



function [img]=PlotConicCurve(C,img,color);

 [m,n]=size(C);
 if m==n & m==3
  C=[C(1,1) C(1,2) C(2,2) C(1,3) C(2,3) C(3,3)];
 end;
 a=C(1);b=C(2);c=C(3);d=C(4);e=C(5);f=C(6);
 
 [m,n,k]=size(img);
 Xn=1:1:n;Yn=1:1:m;
 Ya=[]; Yb=[];
 if c==0
  for i=1:1:n
   Ya=[Ya;Xn(i) -(a*Xn(i)^2+2*d*Xn(i)+f)/(2*(b*Xn(i)+e))];
  end;
 else
  for i=1:1:n
   aux=(b^2-a*c)*Xn(i)^2+2*(b*e-c*d)*Xn(i)+e^2-c*f;
   if aux >= 0
    Ya=[Ya; Xn(i) ((-(b*Xn(i)+e)+sqrt(aux))/c)];
    Yb=[Yb; Xn(i) ((-(b*Xn(i)+e)-sqrt(aux))/c)]; 
   end;
  end;
 end;
 Xa=[]; Xb=[];
 if a==0
  for i=1:1:m
   Xa=[Xa;-(c*Yn(i)^2+2*e*Yn(i)+f)/(2*(b*Yn(i)+d)) Yn(i)];
  end;
 else
  for i=1:1:m
   aux=(b^2-a*c)*Yn(i)^2+2*(b*d-a*e)*Yn(i)+d^2-a*f;
   if aux >= 0
    Xa=[Xa; ((-(b*Yn(i)+d)+sqrt(aux))/a) Yn(i)];
    Xb=[Xb; ((-(b*Yn(i)+d)-sqrt(aux))/a) Yn(i)]; 
   end;
  end;
 end;
 
 points=[];
 [ma,na]=size(Ya);
 for i=1:1:ma
  if round(Ya(i,2))>0 & round(Ya(i,2))<m
   points=[points; Ya(i,:)];
  end;
 end;
 [mb,nb]=size(Yb);
 for i=1:1:mb
  if round(Yb(i,2))>0 & round(Yb(i,2))<m
   points=[points; Yb(i,:)];
  end;
 end;
 [ma,na]=size(Xa);
 for i=1:1:ma
  if round(Xa(i,1))>0 & round(Xa(i,1))<n
   points=[points; Xa(i,:)];
  end;
 end;
 [mb,nb]=size(Xb);
 for i=1:1:mb
  if round(Xb(i,1))>0 & round(Xb(i,1))<n
   points=[points; Xb(i,:)];
  end;
 end;

 [m,n]=size(points);
 for i=1:1:m
  img(round(points(i,2)),round(points(i,1)),1)=color(1);
  img(round(points(i,2)),round(points(i,1)),2)=color(2); 
  img(round(points(i,2)),round(points(i,1)),3)=color(3);
 end;
 
 
 
