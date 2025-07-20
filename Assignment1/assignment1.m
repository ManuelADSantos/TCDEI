clear
close all
clc

% Load given values
load('./test1.mat')

% Plot given values
figure
plot(d)
title('Readings of the experiment'); grid on;
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

% =============== A) ===============
% Calculate mL2 values
m_l2 = (G'*G)\G'*d;
% Is equivalent to m_l2 = inv(G'*G)*G'*d;

% Show mL2
figure()
colors = ["#0072BD","#D95319","#EDB120"];
legs = ["d0 [m]","v0 [m/s]","a0 [m/s^2]"];
lims = [[20 23];[-4.5 -3.2];[0.20 0.55]];
for j=1:3
    subplot(3,1,j)
    plot(m_l2(j,:),'o','MarkerEdgeColor','k','MarkerFaceColor',colors(j))
    grid on; hold on
    xlabel('Repetition');legend(legs(j),'Location','eastoutside');
    xlim([0 51]);ylim(lims(j,:))
end

% Show estimated d values
figure
d_est = G*m_l2;
subplot(2,1,2)
plot(d_est)
ylim([-5 25])
title('Estimated reading values');grid on;
xlabel('Time [s]');ylabel('Distance to wall [m]')
subplot(2,1,1)
plot(d)
ylim([-5 25])
title('Obtained reading values');grid on;
xlabel('Time [s]');ylabel('Distance to wall [m]')

% Calculate p-values
figure
res = d-d_est;
res_2 = res.^2;
xobs = sum(res_2,1);
p = 1-cdf('Chisquare', xobs, m-n);
histogram(p, 10);
title('P-values distribution');
hold on;grid on;
xlabel('Value');ylabel('Number of Ocurrences')

% =============== B) ===============
N=50; % sample size
a=0; % lower boundary
b=1; % higher boundary

nbins = 10; % number of bins
edges = linspace(a,b,nbins+1); % edges of the bins
E = N/nbins*ones(nbins,1); % expected value (equal for uniform dist)

figure
histogram(p,10)
hold on; grid on;
histogram(linspace(a,b,50),10,'FaceAlpha',0.2)
title('Observed p-values distribution VS uniform distribution');
xlabel('Value');ylabel('Number of Ocurrences')
legend('obtained distribution','expected distribution','location','best')

[h,p_result,stats] = chi2gof(p,'Expected',E,'Edges',edges);

% =============== C) ===============
% Compute theoretical value of confidence interval
cov_ml2 = 1.*inv(G'*G);
intervals = 1.96*(diag(sqrt(cov_ml2)));

% =============== D) ===============
% Show confidence intervals
aux_intervals = ones(3,50);
aux_intervals = aux_intervals(1:3,:).*intervals(:);
m_l2_plus = m_l2 + aux_intervals;
m_l2_minus = m_l2 - aux_intervals;
m_true = ones(3,50);
m_true = m_true(1:3,:).*mean(m_l2,2);
figure
for j=1:3
    subplot(3,1,j)
    plot(m_l2(j,:),'o','MarkerEdgeColor','k','MarkerFaceColor',colors(j))
    grid on; hold on
    % plot(m_l2(j,:),'Color',colors(j))
    plot(m_true(j,:),'Color',colors(j))
    plot(m_l2_plus(j,:),'_k')
    plot(m_l2_minus(j,:),'_k')
    for i=1:size(m_l2,2)
        plot([i,i],[m_l2_minus(j,i),m_l2_plus(j,i)],'--k')
    end
    xlabel('Repetition');legend(legs(j),'m_{true}','Location','eastoutside');
    xlim([0 51])
end
sgtitle('Estimated motion parameters with confidence intervals');
