%% NC-L1-SVD DOA 估计仿真 — RMSE vs SNR
clc; clear;

%% ===== 参数设定 =====
N   = 8;
M   = 5;
K   = 200;
SNR = 10;

x_n = (0:N-1)';

%% ===== 单次演示信号生成 =====
theta_true = sort(generate_angles(M, -60, 60, 10));
S          = (randn(M, K) > 0) * 2 - 1;
A          = exp(1j * pi * x_n * sin(deg2rad(theta_true)));
As_k       = A * S;
sigma2     = mean(abs(As_k(:)).^2) / (10^(SNR/10));
V          = sqrt(sigma2/2) * (randn(N,K) + 1j*randn(N,K));
Y          = As_k + V;

grid_step  = 0.1;
grid_theta = -90:grid_step:90;
A_dict     = exp(1j * pi * x_n * sin(deg2rad(grid_theta)));
L          = length(grid_theta);

fprintf('真实角度: %s 度\n', num2str(theta_true));

%% ===== CVX 单次校准 lambda =====
fprintf('正在用CVX校准lambda（仅此一次）...\n');

Y_aug_c = [Y; conj(Y)];
Psi_c   = [A_dict; conj(A_dict)];
[U_c, S_c, ~] = svd(Y_aug_c, 'econ');
s_c     = diag(S_c);
Y_sv_c  = U_c(:,1:M) * S_c(1:M,1:M);

noise_sv_c   = s_c(M+1:end);
sigma2_c     = sum(noise_sv_c.^2) / (length(noise_sv_c) * K);
epsilon_c    = 0.8 * sqrt(sigma2_c * M * K);

cvx_begin quiet
    variable S_cvx(L, M) complex;
    variable r_cvx(L) nonnegative;
    minimize(sum(r_cvx))
    subject to
        norm(Y_sv_c - Psi_c * S_cvx, 'fro') <= epsilon_c;
        for i = 1:L
            norm(S_cvx(i,:), 2) <= r_cvx(i);
        end
cvx_end

S_cvx_val  = double(S_cvx);
resid_cvx  = norm(Y_sv_c - Psi_c * S_cvx_val, 'fro');
sparse_sum = sum(sqrt(sum(abs(S_cvx_val).^2, 2)));
if sparse_sum > 1e-8
    lambda_ref = resid_cvx / (sparse_sum + 1e-8) * 10;
else
    lambda_ref = epsilon_c / (sqrt(M * K) + 1e-6) * 0.3;
end
fprintf('CVX校准完成 → epsilon=%.4f, 残差=%.4f, lambda_ref=%.6f\n', ...
    epsilon_c, resid_cvx, lambda_ref);

%% ===== 公共参数 =====
snr_range       = -20:2:20;
monte_carlo_num = 300;
FAIL_THRESH     = 5;
theta_fixed     = sort(generate_angles(M, -60, 60, 10));

if isempty(gcp('nocreate')); parpool; end

%% ===== RMSE vs SNR =====
rmse_snr = nan(length(snr_range), 1);
crb_snr  = nan(length(snr_range), 1);

fprintf('\n实验：RMSE vs SNR\n');
for s = 1:length(snr_range)
    cur_SNR    = snr_range(s);
    crb_vec    = compute_NC_CRB_internal(theta_fixed, N, cur_SNR, K);
    crb_snr(s) = sqrt(mean(crb_vec));
    all_sq     = nan(monte_carlo_num, M);
    lam        = lambda_ref;

    parfor mc = 1:monte_carlo_num
        S_mc      = (randn(M,K)>0)*2-1;
        A_mc      = exp(1j*pi*(0:N-1)'*sin(deg2rad(theta_fixed)));
        As_mc     = A_mc*S_mc;
        sigma2_mc = mean(abs(As_mc(:)).^2)/(10^(cur_SNR/10));
        V_mc      = sqrt(sigma2_mc/2)*(randn(N,K)+1j*randn(N,K));
        Y_mc      = As_mc + V_mc;

        [~, est_mc] = NC_L1_SVD(Y_mc, A_dict, grid_theta, M, sigma2_mc, lam);
        est_s  = sort(est_mc(:));
        true_s = sort(theta_fixed(:));
        if length(est_s)==M && max(abs(est_s-true_s))<FAIL_THRESH
            all_sq(mc,:) = (est_s - true_s).^2;
        end
    end

    valid = all_sq(~isnan(all_sq(:)));
    if ~isempty(valid)
        rmse_snr(s) = sqrt(mean(valid));
    end
    fprintf('  SNR=%+3ddB  RMSE=%.4f°  CRB=%.4f°\n', cur_SNR, rmse_snr(s), crb_snr(s));
end

%% ===== 绘图 =====
c_rmse = [0.15 0.45 0.78];
c_crb  = [0.85 0.33 0.10];

figure('Color','w','Name','RMSE vs SNR');
ax = axes;

semilogy(snr_range, rmse_snr, '-o', 'Color',c_rmse, 'LineWidth',2, ...
    'MarkerFaceColor',c_rmse, 'MarkerSize',6, 'DisplayName','NC-L1-SVD RMSE');
hold on;
semilogy(snr_range, crb_snr, '--^', 'Color',c_crb, 'LineWidth',2, ...
    'MarkerFaceColor',c_crb, 'MarkerSize',6, 'DisplayName','NC-CRB 理论下界');

% 数值标注
for i = 1:length(snr_range)
    if ~isnan(rmse_snr(i))
        text(ax, snr_range(i), rmse_snr(i)*1.20, sprintf('%.3f', rmse_snr(i)), ...
            'FontSize',6.5, 'Color',c_rmse, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    end
    if ~isnan(crb_snr(i))
        text(ax, snr_range(i), crb_snr(i)*0.78, sprintf('%.3f', crb_snr(i)), ...
            'FontSize',6.5, 'Color',c_crb, ...
            'HorizontalAlignment','center', 'VerticalAlignment','top');
    end
end

grid on; box on;
xlabel('信噪比 (dB)', 'FontSize',11);
ylabel('RMSE (°)', 'FontSize',11);
title('RMSE vs SNR', 'FontSize',12);
legend('Location','northeast', 'FontSize',10);
text(0.04, 0.18, ...
    sprintf('NC-L1-SVD\nK=%d, N=%d, M=%d\n非相干信号', K, N, M), ...
    'Units','normalized', 'FontSize',9, ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[0.95 0.97 1.0], 'EdgeColor',[0.4 0.4 0.8], 'LineWidth',1.0);

fprintf('\n仿真完成。\n');