clear
close all
clc

% Page 40 of "Parameter Estimation and Inverse Problems"

G = ones(10,3);
t = 1:10;
G(:,2) = t;
G(:,3) = -(1/2).*t.^2;

cov_ml2 = (8.^2)*inv(G'*G);

intervals = 1.96.*(diag(sqrt(cov_ml2)))