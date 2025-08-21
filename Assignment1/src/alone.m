clear
close all
clc

% Load given values
load('./test1.mat')

d=d(:,1);

% Plot given values
figure(1);
plot(d)
title('Readings of the first experiment'); grid on;
xlabel('Time [s]');ylabel('Distance to wall [m]')

% Rationale : 50 runs that gather 100 values each
% d = 100 X 50
% => G = 100 X 3
% => m =  3  X 50

m = 100; % 100 observations
n = 3; % 3 model parameters

% Create Ti, to insert in G
time = 0.1:0.1:10;

% Create matrix G
G = ones(size(d,1),3);
G(:,2) = time(:);
G(:,3) = (time(:).^2)/2;

% Calculate mL2 values
m_l2 = (G'*G)\G'*d;
% Is equivalent to m_l2 = inv(G'*G)*G'*d;

% Show mL2
figure(2);
for i = 1:3
    plot(m_l2(i,:),'.')
    hold on;
end
legend('d0 [m]','v0 [m/s]','a0 [m/s^2]','Location','best')
title('Estimated motion parameters');grid on;
xlabel('Test ID');


% Show estimated d values
figure(3);
d_est = G*m_l2;
plot(d_est)
ylim([-5 25])
hold on;
plot(d)
title('Esimated reading values');grid on;
xlabel('Time [s]');ylabel('Distance to wall [m]')

% Proof that the noise follows a standard normal distribution
figure(4);
histogram(normalize(d-d_est,'range',[0 1]),20);
title('Noise distribution');grid on;
xlabel('Value (normalised)');ylabel('Number of Ocurrences')

% Calculate p-values
figure(5);
res = d-d_est;
res_2 = res.^2;
xobs = sum(res_2,1);
p = 1-cdf('Chisquare', xobs, m-n);
histogram(p, 10);
title('P-test values distribution');
hold on;grid on;
% pd = makedist('Uniform','lower',0,'upper',1);
% plot(pd);
xlabel('Value');ylabel('Number of Ocurrences')


% Compute theoretical value of confidence interval
% cov_ml2 = cov(m_l2);
cov_ml2 = 1.*inv(G'*G);
intervals = 1.96*(diag(sqrt(cov_ml2)));

% Show confidence intervals
m_true = mean(m_l2,2);
m_limits = [m_true+intervals,m_true-intervals];
