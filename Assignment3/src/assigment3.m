% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear
close all
clc

% Load given values
load('./test3.mat')

% Plot given values
figure
plot(d,'*'); grid on;
title('Readings of the third experiment');
xlabel('Time [s]');ylabel('Distance to wall [m]');


% ============ Part 1 a) ============ 
m = 100; % 100 observations
n = 3; % 3 model parameters

% Create Ti, to insert in G
time = 0.1:0.1:10;

% Create matrix G
G = ones(size(d,1),3);
G(:,2) = time(:);
G(:,3) = (time(:).^2)/2;

% Calculate mL2 values
m_LS = pinv(G'*G)*G'*d;

% Show estimated d values
d_LS = G*m_LS;
figure
plot(d,'*'); grid on; hold on; plot(d_LS); 
title('Least Squares Solution');
xlabel('Time [s]');ylabel('Distance to wall [m]');
legend('Data','LS solution');


% ============ Part 1 b) ============
res_LS = d-d_LS;
%% 
dof = m-n;
s_LS = vecnorm(res_LS)/sqrt(dof);
Cov_LS = (s_LS.^2).*inv(G'*G);

% close all

% Aplicar Bootstrap
numBoots = 10000;
x_size = length(time);

rng(42);

W = diag(ones(1,100).*s_LS);

% perform the bootstraps
modelfit = zeros(numBoots,x_size);
params = zeros(numBoots,n);
for boot=1:numBoots
  % prepare data indices
  ix = ceil(x_size*rand(1,x_size));
  % construct regressor matrix
  G_boot = W*[ones(x_size,1) time(ix)' (time(ix)'.^2)/2];
  % estimate parameters
  h = pinv(G_boot'*G_boot)*G_boot'*(W*d(ix));
  % evaluate the model
  modelfit(boot,:) = time*h(1) + h(2);
  % record the parameters
  params(boot,:) = h;
end

m_boot = mean(params)';
d_boot = G*m_boot;
res_boot = d-d_boot;

% plot(res_boot,'r+'); grid on; hold on;

aux = ones(numBoots,n);
for i=1:numBoots
    aux(i,:)=m_boot;
end

A = params - aux;
Cov_boot = (A'*A)/numBoots;


% ============ Part 2 ============
outliers=[0.1, 0.3, 0.5];
noise= 0:5;
runs = 100;

% m_all_LS is 3(n = size of m) x 6(noise) x 3(outliers)
m_all_LS = zeros(n,size(noise,2),size(outliers,2));
m_all_L1 = m_all_LS;
m_all_RANSAC = m_all_LS;

for idx_noise = 1:size(noise,2)
    for idx_out = 1:size(outliers,2)
        [d_0, timestamp, sol] = SimulationRobotWall(noise(idx_noise),outliers(idx_out),10,0.1);
        
        % Normal Least Squares
        m_0 = pinv(G'*G)*G'*d_0;
        m_all_LS(:,idx_noise,idx_out) = m_0;
        res_0 = d_0 - G*m_0;
        
        % Robust M-estimation using the L1 Norm
        % Create 100 synthetic datasets
        d_q = zeros(m,runs);
        m_L1 = zeros(n,runs);
        for i = 1:runs
            d_q(:,i) = (G*m_0) + noise(idx_noise) .* randn(m,1);
            m_L1(:,i) = robustfit(G(:,2:3),d_q(:,i));
        end
        m_all_L1(:,idx_noise,idx_out) = mean(m_L1,2);

        % RANSAC
        p = 0.9999;
        num_points = 10;
        e = outliers(idx_out);
        IT = floor( (log(1- p)) / (log(1-(1-e).^(num_points))));
        threshold = 10;

        for cont=1:IT
            ind = randsample(length(d_0),num_points);
            flag = 0;
            data_hinliner = d_0(ind,:);
            Gr=G(ind,:);
            while(1)
                comp = isoutlier(data_hinliner);
                data_hinliner = randsample(d_0,num_points);
                if(comp == zeros(1,num_points))
                    break;
                end
            end

            m_RANSAC = pinv(Gr'*Gr)*Gr'*data_hinliner;
            d_RANSAC = G*m_RANSAC;
            
            for j=1:size(d_RANSAC,1)
                if(norm(d_RANSAC(j)- d_0(ind))>threshold)
                    flag=1;
                end
            end

            if (flag==0)
                m_all_RANSAC(:,idx_noise,idx_out) = m_RANSAC;
                break;
            end
        end

    end
end

for i = 1:3
    figure 
    plot([m_all_LS(:,:,i),sol])
    legend('noise = 0','noise = 1','noise = 2','noise = 3','noise = 4','noise = 5','ground truth')
    hold on
    grid on
end
% plot(sol,'LineWidth',3)


