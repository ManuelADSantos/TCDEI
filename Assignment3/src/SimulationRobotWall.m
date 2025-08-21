 % Function to simulate sonar readings of a robot 
% moving towards a wall (assignment 1 of CMEDI)
%
% [d,time,sol]=SimulationRobotWall(NoiseStd,OutlierRate,ExpDuration,TimeStep);
%
% Inputs:
% NoiseStd - Noise standard deviation (in meters)
% OutlierRate - Percentage of outliers in the data set (from 0 to 1)
% ExpDuration - Time duration of the experiment (in seconds)
% TimeStep - Sampling Period (in seconds)
%
% Outputs:
% d - set of sonar readings (in meters)
% time - sampling instants (in seconds)
% sol - ground truth [initial position, initial velocity, acceleration]
%
% ---------
% Joao P. Barreto / University of Coimbra. 
% This software is exclusively meant for academic use. 
% Send your feedback to jpbar@deec.uc.pt

function [d,time,sol]=SimulationRobotWall(NoiseStd,OutlierRate,ExpDuration,TimeStep);

 go=1;
 while go
  r=rand(1,3);
  v=r(1)*(-1)+(1-r(1))*(-4);
  a=r(2)*(-0.8*v/ExpDuration)+(1-r(2))*(-v/ExpDuration);
  aux=-0.5*a*ExpDuration^2-v*ExpDuration;
  if aux<40
   s=r(3)*aux+(1-r(3))*40;
   go=0;
  end;
 end;
 sol=[s; v; a];
 time=transpose(0:TimeStep:ExpDuration-TimeStep);
 m=max(size(time));
 
 % Data Points
 gt=s+v*time+0.5*a*time.^2;
 d=gt+NoiseStd*randn(m,1);
 % Outliers
 n=round(OutlierRate*m);
 go=1;
 while go & n
  order=round(m*rand(1,n));
  aux=find(order==0);
  if isempty(aux)
   go=0;
  end;
 end;
 d(order)=gt(order)+25*randn(n,1);
 
 
      
