%% ================== NC-ESPRIT Monte Carlo + NC-CRB ==================
clear; clc; close all;

%% ---------- 固定参数 ----------
M      = 8;
K      = 5;
N_snap = 50;
N_mc   = 1000;
SNR_dB = -20:2:20;
n_snr  = length(SNR_dB);

%% ---------- 固定角度 ----------
while true
    theta_true = sort(-60 + 120 * rand(1, K));
    if min(diff(theta_true)) >= 5; break; end
end
theta_rad = deg2rad(theta_true);
fprintf('固定角度: '); fprintf('%.2f°  ', theta_true); fprintf('\n\n');

%% ---------- 阵列流形 ----------
m = (0:M-1)';
A = exp(1j * pi * m * sin(theta_rad));   % M×K

%% ---------- 预分配 ----------
rmse_nc = zeros(1, n_snr);
crb_nc  = zeros(1, n_snr);

%% ---------- 主循环 ----------
for si = 1:n_snr
    snr    = SNR_dB(si);
    sigma2 = 10^(-snr/10);

    %% CRB
    crb_each   = compute_NC_CRB(theta_true, M, N_snap, sigma2);
    crb_nc(si) = sqrt(mean(crb_each.^2));

    %% Monte Carlo
    err2 = zeros(N_mc, K);
    n_fail = 0;

    for mc = 1:N_mc
        S     = randn(K, N_snap);
        noise = sqrt(sigma2/2) * (randn(M,N_snap) + 1j*randn(M,N_snap));
        Y     = A * S + noise;

        try
            est = NC_ESPRIT(Y, K);           % 1×K，已排序
            err2(mc,:) = (est(:)' - sort(theta_true)).^2;
        catch
            n_fail = n_fail + 1;
            err2(mc,:) = NaN;
        end
    end

    valid = sum(~isnan(err2(:,1)));
    if valid > 0
        rmse_nc(si) = sqrt(mean(err2(~isnan(err2(:,1)),:), 'all'));
    else
        rmse_nc(si) = NaN;
    end

    fprintf('SNR=%+4ddB | RMSE=%.4f° | CRB=%.4f° | 失败=%d/%d\n', ...
        snr, rmse_nc(si), crb_nc(si), n_fail, N_mc);
end

%% ---------- 绘图 ----------
figure('Color','w','Position',[100 100 750 500]);

semilogy(SNR_dB, rmse_nc, 'b-o', ...
    'LineWidth',2,'MarkerSize',6,'MarkerFaceColor','b','DisplayName','NC-ESPRIT RMSE');
hold on;
semilogy(SNR_dB, crb_nc, 'r--s', ...
    'LineWidth',2,'MarkerSize',6,'MarkerFaceColor','r','DisplayName','NC-CRB');

grid on;
xlabel('SNR (dB)','FontSize',12);
ylabel('RMSE / CRB (°)','FontSize',12);
legend('Location','southwest','FontSize',11);
title({sprintf('NC-ESPRIT Monte Carlo vs NC-CRB'), ...
       sprintf('M=%d  |  K=%d  |  N=%d 快拍  |  %d 次Monte Carlo', ...
               M, K, N_snap, N_mc)}, 'FontSize',11);
xlim([SNR_dB(1) SNR_dB(end)]);
ylim([1e-3 1e2]);