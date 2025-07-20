% Function comparing conic Omega1 with conic Omega2. The output is the
% mean distance between the principal points
% 
% result=CompareConic(Omega1,Omega2)
% 
% ---------
% Joao P. Barreto / University of Coimbra. 
% This software is exclusively meant for academic use. 
% Send your feedback to jpbar@deec.uc.pt


function result=CompareConic(Omega1,Omega2)

%Compute Principal Points of Omega1
conic=[Omega1(1,1) Omega1(1,2) Omega1(2,2)];
lambda_1=0.5*(-conic(1)+conic(3)+sqrt(conic(1)^2-2*conic(1)*conic(3)+conic(3)^2+4*conic(2)^2))/conic(2);
lambda_2=0.5*(-conic(1)+conic(3)-sqrt(conic(1)^2-2*conic(1)*conic(3)+conic(3)^2+4*conic(2)^2))/conic(2);
axis_1=transpose(Omega1(1,:))+lambda_1*transpose(Omega1(2,:));
axis_2=transpose(Omega1(1,:))+lambda_2*transpose(Omega1(2,:));
P_1=inter_line_conic(Omega1,axis_1);
P_2=inter_line_conic(Omega1,axis_2);
if isreal(P_1)
 norma_1=norm(P_1(1:2,1)-P_1(1:2,2));
else
 norma_1=Inf;
end;
if isreal(P_2)
 norma_2=norm(P_2(1:2,1)-P_2(1:2,2));
else
 norma_2=Inf;
end;
if norma_1>norma_2
 axis_M=axis_1; P_M_1=P_1;
 axis_m=axis_2; P_m_1=P_2;
else
 axis_M=axis_2; P_M_1=P_2;
 axis_m=axis_1; P_m_1=P_1;
end;
%Compute Principal Points for Omega2
conic=[Omega2(1,1) Omega2(1,2) Omega2(2,2)];
lambda_1=0.5*(-conic(1)+conic(3)+sqrt(conic(1)^2-2*conic(1)*conic(3)+conic(3)^2+4*conic(2)^2))/conic(2);
lambda_2=0.5*(-conic(1)+conic(3)-sqrt(conic(1)^2-2*conic(1)*conic(3)+conic(3)^2+4*conic(2)^2))/conic(2);
axis_1=transpose(Omega2(1,:))+lambda_1*transpose(Omega2(2,:));
axis_2=transpose(Omega2(1,:))+lambda_2*transpose(Omega2(2,:));
P_1=inter_line_conic(Omega2,axis_1);
P_2=inter_line_conic(Omega2,axis_2);
if isreal(P_1)
 norma_1=norm(P_1(1:2,1)-P_1(1:2,2));
else
 norma_1=Inf;
end;
if isreal(P_2)
 norma_2=norm(P_2(1:2,1)-P_2(1:2,2));
else
 norma_2=Inf;
end;
if norma_1>norma_2
 axis_M=axis_1; P_M_2=P_1;
 axis_m=axis_2; P_m_2=P_2;
else
 axis_M=axis_2; P_M_2=P_2;
 axis_m=axis_1; P_m_2=P_1;
end;
% Compute Error
distances=[];
if isreal(P_M_1) & isreal(P_M_2)
 aux=P_M_1-P_M_2;
 aux_a=sqrt(diag(transpose(aux)*aux));
 aux=P_M_1-P_M_2*[0 1;1 0];
 aux_b=sqrt(diag(transpose(aux)*aux));
 if ([1 1]*aux_a) > ([1 1]*aux_b)
  distances=[distances transpose(aux_b)];
 else
  distances=[distances transpose(aux_a)];
 end;
end;
aux=P_m_1-P_m_2;
aux_a=sqrt(diag(transpose(aux)*aux));
aux=P_m_1-P_m_2*[0 1;1 0];
aux_b=sqrt(diag(transpose(aux)*aux));
if ([1 1]*aux_a) > ([1 1]*aux_b)
 distances=[distances transpose(aux_b)];
else
 distances=[distances transpose(aux_a)];
end;
result=mean(distances);







function result=inter_line_conic(Oi,r);

[m,n]=size(Oi);

if m==n & m==3
 Omega=Oi;
 Oi=[Omega(1,1) Omega(1,2) Omega(2,2) Omega(1,3) Omega(2,3) Omega(3,3)];
 Omega_conj=[Oi(6)*Oi(3)-Oi(5)^2 Oi(5)*Oi(4)-Oi(2)*Oi(6) Oi(2)*Oi(5)-Oi(3)*Oi(4);Oi(5)*Oi(4)-Oi(2)*Oi(6) Oi(1)*Oi(6)-Oi(4)^2 Oi(2)*Oi(4)-Oi(1)*Oi(5);Oi(2)*Oi(5)-Oi(3)*Oi(4) Oi(2)*Oi(4)-Oi(1)*Oi(5) Oi(1)*Oi(3)-Oi(2)^2];
else
 Omega=[Oi(1) Oi(2) Oi(4);Oi(2) Oi(3) Oi(5);Oi(4) Oi(5) Oi(6)];
 Omega_conj=[Oi(6)*Oi(3)-Oi(5)^2 Oi(5)*Oi(4)-Oi(2)*Oi(6) Oi(2)*Oi(5)-Oi(3)*Oi(4);Oi(5)*Oi(4)-Oi(2)*Oi(6) Oi(1)*Oi(6)-Oi(4)^2 Oi(2)*Oi(4)-Oi(1)*Oi(5);Oi(2)*Oi(5)-Oi(3)*Oi(4) Oi(2)*Oi(4)-Oi(1)*Oi(5) Oi(1)*Oi(3)-Oi(2)^2];
end; 
r=[r(1);r(2);r(3)];

I=eye(3,3);
Is=[0 0 -1; 0 0 -1; 1 1 0];
P1=(sqrt(-transpose(r)*Omega_conj*r)*I+skew_symetric_v(r)*Omega)*Is*r;
P2=(-sqrt(-transpose(r)*Omega_conj*r)*I+skew_symetric_v(r)*Omega)*Is*r;
P1=P1*(1/P1(3)); P2=P2*(1/P2(3));
result=[P1 P2];


function result=skew_symetric_v(a);
result=[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];


