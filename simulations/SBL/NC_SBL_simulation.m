%% ===== NC-SBL + CRB 蒙特卡洛仿真 =====
clear; clc;

%% 强制重启并行池
delete(gcp('nocreate'));
parpool('local');

%% 1. 参数配置
M          = 8;
N          = 800;
K_true     = 2;
grid_theta = -90:1:90;    % 粗网格（算法自行细化）
SNR_list   = -20:2:10;
Monte      = 200;

nSNR     = length(SNR_list);
rmse_sbl = zeros(1, nSNR);
crb_val  = zeros(1, nSNR);

%% 2. SNR循环
for s = 1:nSNR
    SNR_dB  = SNR_list(s);
    snr_lin = 10^(SNR_dB/10);

    rmse_mc = zeros(1, Monte);
    crb_mc  = zeros(1, Monte);
    success_mc = zeros(1, Monte); 
    parfor mc = 1:Monte

        %每次Monte Carlo随机生成非整数角度
        theta_true = generate_angles(K_true, -60, 60, 10);

        A = exp(1j*pi*(0:M-1)' * sind(theta_true));
        S = 2*randi([0,1], K_true, N) - 1;

        sig_pow = mean(abs(A*S).^2, 'all');
        noi_pow = sig_pow / snr_lin;
        noise   = sqrt(noi_pow/2) * (randn(M,N) + 1j*randn(M,N));
        Y       = A*S + noise;

        est_theta = NC_SBL(Y, grid_theta);

        true_s = sort(theta_true);
        if ~isempty(est_theta)
            % 对所有检测峰，找与真实角度最优匹配的K_true个
            if length(est_theta) >= K_true
                est_s = sort(est_theta(1:K_true));
            else
                est_s = sort(est_theta);
                % 不足K_true个，用最近网格点补齐
                est_s(end+1:K_true) = grid_theta(1);  % 惩罚补齐
            end
            true_s  = sort(theta_true);
            min_mse = Inf;
            p_all   = perms(1:K_true);
            for pi_ = 1:size(p_all,1)
                mse_p = mean((est_s(p_all(pi_,:)) - true_s).^2);
                if mse_p < min_mse, min_mse = mse_p; end
            end
            rmse_mc(mc)    = sqrt(min_mse);
            success_mc(mc) = 1;
        else
            rmse_mc(mc)    = mean(diff(grid_theta))/2;
            success_mc(mc) = 0;
        end

        crb_mc(mc) = sqrt(mean(compute_NC_CRB(theta_true, N, noi_pow, M)));
    end

    rmse_sbl(s) = mean(rmse_mc);   %  平均RMSE
    crb_val(s)  = mean(crb_mc);

   fprintf('SNR=%3d dB | RMSE=%.4f° | sqrt(CRB)=%.4f° | 成功率=%.1f%%\n', ...
    SNR_dB, rmse_sbl(s), crb_val(s), 100*mean(success_mc));
end
%% 3. 绘图
figure('Color','w','Position',[200 200 800 500]);
semilogy(SNR_list, rmse_sbl, 'b-s','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
hold on;
semilogy(SNR_list, crb_val, 'r--','LineWidth',2);
grid on;
xlabel('SNR (dB)','FontSize',12);
ylabel('RMSE (°)','FontSize',12);
title(sprintf('NC-SBL DOA估计 vs NC-CRB  (M=%d, N=%d, K=%d)', M, N, K_true));
legend('NC-SBL RMSE','NC-CRB 理论下界','Location','northeast');
set(gca,'FontSize',11);


function crb_vec = compute_NC_CRB(theta_deg, N, sigma2, M)
    K         = length(theta_deg);
    theta_rad = theta_deg * pi/180;
    m         = (0:M-1)';
    A = zeros(M,K); D = zeros(M,K);
    for k = 1:K
        a_k    = exp(1j*pi*m*sin(theta_rad(k)));
        A(:,k) = a_k;
        D(:,k) = 1j*pi*m*cos(theta_rad(k)) * (pi/180) .* a_k;
    end

    A_e = [A; conj(A)];
    D_e = [D; conj(D)];

    R  = A_e*A_e' + sigma2*eye(2*M);  % 2M×2M
    Ri = inv(R);

    FIM = zeros(K,K);
    for i = 1:K
        dRi = D_e(:,i)*A_e(:,i)' + A_e(:,i)*D_e(:,i)';
        for j = 1:K
            dRj = D_e(:,j)*A_e(:,j)' + A_e(:,j)*D_e(:,j)';
            FIM(i,j) = 2*N*real(trace(Ri*dRi*Ri*dRj));
        end
    end
    crb_vec = diag(inv(FIM))';
end