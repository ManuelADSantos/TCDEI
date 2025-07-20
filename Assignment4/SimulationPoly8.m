% Function to simulate readings of o polynomial of 8th degree
% as function of a scalar parameter - theta (assignment 4 of TCDEI)
%
% [y,theta,sol]=SimulationPoly8(m,NoiseStd,OutlierRate,draw);
%
% Inputs:
% m - Number of data points to generate
% NoiseStd - Noise standard deviation (in meters)
% OutlierRate - Percentage of outliers in the data set (from 0 to 1)
% draw - if draw then plots are produced for data visualization
%
% Outputs:
% y - set of readings of the dependent variable
% theta - sampling values of the independent variable
% sol - ground truth coefficients [a,b,c,d,e,f,g,h]
%
% ---------
% Nuno Gonçalves / University of Coimbra. 
% This software is exclusively meant for academic use. 
% Send your feedback to nunogon@deec.uc.pt

function [y,theta,sol]=SimulationPoly8(m,NoiseStd,OutlierRate,draw)

lower_bound=-2.5;
upper_bound=1.25;

a=-0.2;
b=-0.5;
c=0.4;
d=-0.1;
e=-2;
f=3.0;
g=1.75;
h=-2.5;

sol=[a;b;c;d;e;f;g;h];

theta=lower_bound+(upper_bound-lower_bound)*rand(m,1);
y=a*theta.^8+b*theta.^7+c*theta.^6+d*theta.^5+e*theta.^4+f*theta.^3+g*theta.^2+h*theta;

if draw
    figure;
    plot(theta,y,'*');
    hold on;
end

% Data Points
y=y+NoiseStd*randn(m,1);
if draw
    plot(theta,y,'.');
end

% Outliers
n=round(OutlierRate*m);
if n<=0
    return;
end

go=1;
while go && n
 order=round(m*rand(1,n));
 aux=find(order==0);
 if isempty(aux)
  go=0;
 end
end
y(order)=y(order)+25*randn(n,1);
if draw
    plot(theta,y,'o');
end
