% Function to randomly generate conic curves 
% and noisy points over the curve (assignment 2 of CMEDI)
%
% [Points,Omega_gt,img]=CurveGenerator(noise_std,n_points,conic_type)
%
% Inputs:
% noise_std = Standard deviation of the aditive gaussian noise
%             added in X amd Y to each point
% n_points = number of sample points to be generated
% conic_type = 0 any type 
%              1 only elipses
%              2 only hyperbols
%
% Outputs:
% Points - 2xn_points matrix with the noisy point coordinates
% Omega_gt - 3x3 matrix representing the Conic Curve (Ground Truth)
% img - 640x480 image with the Conic Curve
%
% ---------
% Joao P. Barreto / University of Coimbra. 
% This software is exclusively meant for academic use. 
% Send your feedback to jpbar@deec.uc.pt



function [Points,Omega_gt,img]=CurveGenerator(noise_std,n_points,conic_type)


qsi=0.55;
FOV=85*pi/180;
N=640;M=480;

% Auxiliay Computation
H_c=[N/(2*tan(FOV)) 0 N/2;0 M/(2*tan(FOV)) M/2;0 0 1];
cangle=pi/2-asin(qsi);

% Compute Normal
theta=2*pi*rand;
phi=pi/2*rand;
if conic_type==1
 while phi<=cangle
  phi=pi/2*rand;
 end;
elseif conic_type==2
 while phi>=(cangle) 
  phi=pi/2*rand;
 end;
end;
n=[cos(phi)*sin(theta) cos(phi)*cos(theta) sin(phi)];
% Compute Conic
Omega_gt=transpose(inv(H_c))*[n(1)^2*(1-qsi^2)-n(3)^2*qsi^2 n(1)*n(2)*(1-qsi^2) n(1)*n(3);
		    n(1)*n(2)*(1-qsi^2) n(2)^2*(1-qsi^2)-n(3)^2*qsi^2 n(2)*n(3); 
            n(1)*n(3) n(2)*n(3) n(3)^2]*inv(H_c);
Omega_gt=Omega_gt*Omega_gt(3,3)^-1;
%Computing Points
img=255*ones(M,N,3);
[Points,img]=PlotConicCurveN(Omega_gt,img,[0 255 0]);
n_points=min(max(size(Points)),n_points);
I=ceil(max(size(Points))*rand(1,n_points));
aux=find(I==0);
I(aux)=1;
Points=transpose(Points(I,:))+noise_std*(randn(2,n_points));

function [points,img]=PlotConicCurveN(C,img,color);

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
 


