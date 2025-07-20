% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear
close all
clc

rng(42);

m = 100; % 100 observations
n = 3; % 3 model parameters

% Create Ti, to insert in G
time = 0.1:0.1:10;
% Create matrix G
G = ones(m,3);
G(:,2) = time(:);
G(:,3) = (time(:).^2)/2;


% ============ Part 2 ============
outliers=[0.1, 0.3, 0.5];
noise= 0:5;
runs = 100;

idx_noise = 3; % noise = 2
idx_out = 1;   % outliers = 0.1 

[d_0, timestamp, sol] = SimulationRobotWall(noise(idx_noise),outliers(idx_out),10,0.1);

% Normal Least Squares
m_0 = pinv(G'*G)*G'*d_0;
%         res_0 = d_0 - G*m_0;

% Create 100 synthetic datasets
% d_q = zeros(m,runs);
% m_L1 = zeros(n,runs);
% d_L1 = zeros(m,runs);
% res_L1 = zeros(m, runs);     

thresh = 1; max_iter = 200;
for i = 1:runs
    d_q(:,i) = (G*m_0) + noise(idx_noise) .* randn(m,1);
    
    m_L1_begin = pinv(G'*G)*G'*d_q(:,i);
    res_q = d_0 - d_q(:,1);
    iter = 1;
    R = diag(1./abs(res_q));
    while iter < max_iter
        m_L1_now = pinv(G.'*R*G)*(G.'*R*d_q(:,1));
        res_L1_now = d_q(:,1) - G*m_L1_now;

        rms_reg(iter) = rms(res_L1_now);

        if rms(res_L1_now)<thresh
%             m_L1(:,i) = m_L1_now;
            break
        end

        iter = iter + 1;
        R = diag(1./abs(res_L1_now));
    end
    m_L1(:,i) = m_L1_now;


end
mean(m_L1,2)

% m_L1_begin = pinv(G'*G)*G'*d_q;
% res_q = d_0 - d_q;


% m_L1(:,i) = robustfit(G(:,2:3),d_q(:,i));
% d_L1(:,i) = G*m_L1(:,i);
% res_L1(:,i) = d_q(:,i) - d_L1(:,i);


%RMS
% sqrt(mean(sum(res.^2))      


