% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear
close all
clc

% =============== A) ===============
% Load given values
load('./test2.mat')

% Plot given values
figure
plot(d,'*')
title('Readings of the second experiment'); grid on;
xlabel('Time [s]');ylabel('Distance to wall [m]')

m = 100; % 100 observations
n = 3; % 3 model parameters


% =============== B) ===============
% Create Ti, to insert in G
time = 0.1:0.1:10;

% Create matrix G
G = ones(size(d,1),3);
G(:,2) = time(:);
G(:,3) = (time(:).^2)/2;

% Calculate mL2 values
m_l2 = pinv(G'*G)*G'*d;

% Show estimated d values
figure
d_est = G*m_l2;
plot(d,'*')
xlabel('Time [s]');ylabel('Distance to wall [m]')
hold on; grid on;
plot(d_est);
title('Obtained VS Estimated');
legend('Obtained d values','Estimated d values')

% Calculate residues
residuals = d-d_est;
figure; grid on; hold on;
plot(residuals);
title('Residuals');

% Quantile-quantile test for residuals
figure
qqplot(sort(residuals)); grid on;


% =============== C) ===============
% Calcular desvio padrão estimado
dof = m-n;
s = vecnorm(residuals)/sqrt(dof);
s_i = time.*s;

% Calculate parameters with weighted expressions
W = diag(1./s_i);
GW = W*G;
dW = W*d;
m_l2_W = pinv(GW'*GW)*GW'*dW;
d_est_W = G*m_l2_W;
residuals_W = dW - GW*m_l2_W; % same as (d-G*m_l2_W)_i./sigma_i'

% Plot estimated d values vs original ones
figure; hold on; grid on;
plot(d,'*');plot(d_est); plot(d_est_W);
xlabel('Time [s]');ylabel('Distance to wall [m]')
title('Obtained VS Estimated VS Estimated Weighted');
legend('Obtained d values','Estimated non weighted d values','Estimated weighted d values')

% Quantile-quantile test for weighted residues
figure
qqplot(sort(residuals_W)); grid on;

% Residues before and after weighting
figure;hold on; grid on;
plot(residuals); plot(residuals_W)
title('Residuals comparison')
legend('Before weighted','After weighted','Location','northwest')


% =============== D) ===============
% res_2 = residuals_W.^2;
% xobs = sum(res_2,1);
% p = 1-cdf('Chisquare', xobs, m-n);


% =============== E) ===============
% Calculate covariance matrix with s approximation
C = (s.^2).*inv(GW'*GW);

% Calculate intervals
intervals = 1.96*(diag(sqrt(C)));

% Show confidence intervals
m_l2_plus = m_l2_W + intervals;
m_l2_minus = m_l2_W - intervals;
figure
colors = ["#0072BD","#D95319","#EDB120"];
legs = ["d0 [m]","v0 [m/s]","a0 [m/s^2]"];
lims = [[28 33];[-6 -4];[0 1]];
for j=1:3
    subplot(3,1,j);
    plot(m_l2_W(j),'o','MarkerEdgeColor','k','MarkerFaceColor',colors(j));
    grid on; hold on;
    plot(m_l2_plus(j),'_k'); plot(m_l2_minus(j),'_k');
    plot([1,1],[m_l2_minus(j),m_l2_plus(j)],'--k');
    legend(legs(j),'Location','eastoutside');
    xlim([0 2]); ylim(lims(j,:));
end
sgtitle('Estimated motion parameters with confidence intervals');
