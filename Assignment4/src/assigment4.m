% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear
close all
clc

rng(42);

% number of points: 10, 16, 32
% noise standard deviations: 0.015, 0.15, 0.3
% y = a.θ^8 + b.θ^7 + c.θ^6 + d.θ^5 + e.θ^4 + f.θ^3 + g.θ^2 + h.θ

numPoints = 16;
noiseStd = 0.15;
outliersPerc = 0;


% =============== 1) ===============
[d,theta,sol] = SimulationPoly8(numPoints,noiseStd,outliersPerc,false);

G = [theta.^8,theta.^7,theta.^6,theta.^5,theta.^4,theta.^3,theta.^2,theta];
p = rank(G);
[m,n] = size(G); 

m_LS = pinv(G'*G)*G'*d;
d_ls = G*m_LS;

figure
plot(G*sol); grid on; hold on;
plot(d_ls); legend('Ground Truth', 'LS solution','Location','north');


% =============== 2) ===============
[U,S,V] = svd(G);

Vp = V(:,1:p);
Up = U(:,1:p);
Sp = S(1:p,1:p);

G_dag = Vp*inv(Sp)*Up';
m_dag = G_dag*d;
d_dag = G*m_dag;

figure
plot(G*sol); grid on; hold on; 
plot(d); plot(d_ls,'o'); plot(d_dag,'x');
legend('Ground Truth','Data','LS solution','SVD solution');
xlim([0 m+1]); ylim([min([min(d),min(d_ls),min(d_dag)])-1 max([max(d),max(d_ls),max(d_dag)])+1]);


% =============== 3) ===============
dof = m-n;
s = vecnorm((d-d_dag))/sqrt(dof);
% Calculate covariance matrix with s approximation
% Cov = (s.^2).*inv(G'*G);

Cov = zeros(n,n);
for i=1:p
    Cov = Cov + V(:,i)*V(:,i)'/S(i,i);
end
Cov = (s.^2).*Cov;

% Calculate intervals
% inverse cumulative distribution
icd = tinv(0.975,dof);
intervals = icd*(diag(sqrt(Cov)));
% Como eh a distribuiçao de student, para obter 5% retiramos 2.5% de cada lado (eh simetrica)
m_dag_plus = m_dag + intervals;
m_dag_minus = m_dag - intervals;

figure
errorbar(m_dag,intervals,'*'); 
grid on; hold on; plot(sol,'+');
xlim([0 9]); title("Confidence Intervals");
legend('Confidence Interval of m dagger','Ground Truth')


% =============== 4) ===============
sing_values = diag(S);
k = min(m,n);
cond_number = sing_values(1)./sing_values(k);


% =============== 5) ===============
picard = zeros(1,n);
p_TSVD = 0;
done = false;
for i = 1:n
    picard(i) = (U(:,i)'*d)/sing_values(i);
    if abs(picard(i)) <= 1 && ~done
        p_TSVD = i;
    else
        done = true;
    end
end

figure
plot(abs(picard)); grid on; hold on;
plot(abs(picard),'*'); plot(ones(1,n),'--k');
title("Discrete Picard Condition analysis");

G_TSVD = Vp(:,1:p_TSVD)*inv(Sp(1:p_TSVD,1:p_TSVD))*Up(:,1:p_TSVD)';
m_TSVD = G_TSVD*d;
d_TSVD = G*m_TSVD; 


% =============== 6) ===============
Rm = Vp*Vp';
figure
heatmap(Rm); title('SVD Model Resolution')

Rm_trunc = Vp(:,1:p_TSVD)*Vp(:,1:p_TSVD)';
figure 
heatmap(Rm_trunc); title('Truncated SVD Model Resolution')


%% =============== 7) ===============
% clear
% close all
% clc

[d_7,theta_7,sol_7] = SimulationPoly8(32,0.15,0.1,false);

G_7 = [theta_7.^8,theta_7.^7,theta_7.^6,theta_7.^5,theta_7.^4,theta_7.^3,theta_7.^2,theta_7];
p_7 = rank(G_7);

[U_7,S_7,V_7] = svd(G_7);
Vp_7 = V_7(:,1:p_7);
Up_7 = U_7(:,1:p_7);
Sp_7 = S_7(1:p_7,1:p_7);

G_dag_7 = Vp_7*inv(Sp_7)*Up_7';
m_dag_7 = G_dag_7*d_7;
d_dag_7 = G_7*m_dag_7;

picard_7 = zeros(1,size(G_7,2));
p_TSVD_7 = 0;
done_7 = false;
sing_values_7 = diag(S_7);
for i = 1:size(G_7,2)
    picard_7(i) = (U_7(:,i)'*d_7)/sing_values_7(i);
    if abs(picard_7(i)) <= 1 && ~done_7
        p_TSVD_7 = i;
    else
        done_7 = true;
    end
end

figure
plot(abs(picard_7)); grid on; hold on;
plot(abs(picard_7),'*'); plot(ones(1,size(G_7,2)),'--k');
title("Discrete Picard Condition analysis");

G_TSVD_7 = Vp_7(:,1:p_TSVD_7)*inv(Sp_7(1:p_TSVD_7,1:p_TSVD_7))*Up_7(:,1:p_TSVD_7)';
m_TSVD_7 = G_TSVD_7*d_7;
d_TSVD_7 = G_7*m_TSVD_7;

figure 
plot(G_7*sol_7); hold on; grid on;
plot(d_7);plot(d_dag_7);plot(G_7*m_TSVD_7)
legend('Ground Truth', 'Data', 'Generalised Inverse Solution', 'Truncated Solution')