%% NC-SBL 相干源 vs 非相干源 对比仿真
clear; clc; close all;

if isempty(gcp('nocreate'))
    parpool('local', 4);
end

%% ================== 参数==================
M_arr      = 8;
K          = 3;                    % 信源数
theta_true = sort(generate_angles(K,-60,60,8));        %随机3个信号角度
grid_theta = -90:1:90;             %粗网格，算法自身细化

n_MC       = 200;                 %蒙特卡洛仿真次数
N_snap     = 200;                 %快拍数
rho        = 0.9;                  % 相干系数

SNR_vec    = -10:2:20;
n_snr      = length(SNR_vec);

rmse_incoh = nan(1, n_snr);
rmse_coh   = nan(1, n_snr);
rate_incoh = zeros(1, n_snr);     % 成功检测率
rate_coh   = zeros(1, n_snr);

%% ================== 主循环 ==================
for si = 1:n_snr
    snr_db = SNR_vec(si);
    sigma2 = 10^(-snr_db/10);

    err_incoh  = nan(1, n_MC);
    err_coh    = nan(1, n_MC);
    ok_incoh   = false(1, n_MC);
    ok_coh     = false(1, n_MC);

    parfor mc = 1:n_MC
        A     = exp(1j*pi*(0:M_arr-1).' * sind(theta_true));  % M×K
        noise = sqrt(sigma2/2) * (randn(M_arr,N_snap) + 1j*randn(M_arr,N_snap));

        %% 非相干：K路独立BPSK
        S_incoh = sign(randn(K, N_snap));
        Y_incoh = A * S_incoh + noise;

        %% 相干：前两路相干，第三路独立
        s1    = sign(randn(1, N_snap));
        S_coh = [s1; rho*s1 + sqrt(1-rho^2)*sign(randn(1,N_snap)); sign(randn(1,N_snap))];
        Y_coh = A * S_coh + noise;   % 同一噪声，公平对比

        %% 非相干估计
        try
            est = NC_SBL(Y_incoh, grid_theta);
            if length(est) == K
                err_incoh(mc) = mean((sort(theta_true) - sort(est)).^2);
                ok_incoh(mc)  = true;
            end
        catch; end

        %% 相干估计
        try
            est = NC_SBL(Y_coh, grid_theta);
            if length(est) == K
                err_coh(mc) = mean((sort(theta_true) - sort(est)).^2);
                ok_coh(mc)  = true;
            end
        catch; end
    end

    % 成功率
    rate_incoh(si) = mean(ok_incoh);
    rate_coh(si)   = mean(ok_coh);

    % RMSE（只对成功的MC取均值，分母用成功次数，才是真实误差）
    if any(ok_incoh)
        rmse_incoh(si) = sqrt(mean(err_incoh(ok_incoh)));
    end
    if any(ok_coh)
        rmse_coh(si) = sqrt(mean(err_coh(ok_coh)));
    end

    fprintf('SNR=%+3ddB  非相干成功率=%.0f%%  相干成功率=%.0f%%\n', ...
        snr_db, rate_incoh(si)*100, rate_coh(si)*100);
end

%% ================== (2) RMSE vs 快拍数 ==================
SNR_snap        = 10;
sigma2_sn       = 10^(-SNR_snap/10);
Nsnap_vec       = [1,2,5,10,20,50,100,200,500,1000];
rmse_incoh_snap = nan(1, length(Nsnap_vec));
rmse_coh_snap   = nan(1, length(Nsnap_vec));

for ni = 1:length(Nsnap_vec)
    N_s         = Nsnap_vec(ni);
    err_incoh_s = nan(1, n_MC);
    err_coh_s   = nan(1, n_MC);
    ok_incoh_s  = false(1, n_MC);
    ok_coh_s    = false(1, n_MC);

    parfor mc = 1:n_MC
        A     = exp(1j*pi*(0:M_arr-1).' * sind(theta_true));
        noise = sqrt(sigma2_sn/2) * (randn(M_arr,N_s) + 1j*randn(M_arr,N_s));

        S_incoh = sign(randn(K, N_s));
        Y_incoh = A * S_incoh + noise;

        s1    = sign(randn(1, N_s));
        S_coh = [s1; rho*s1 + sqrt(1-rho^2)*sign(randn(1,N_s)); sign(randn(1,N_s))];
        Y_coh = A * S_coh + noise;

        try
            est = NC_SBL(Y_incoh, grid_theta);
            if length(est) == K
                err_incoh_s(mc) = mean((sort(theta_true) - sort(est)).^2);
                ok_incoh_s(mc)  = true;
            end
        catch; end

        try
            est = NC_SBL(Y_coh, grid_theta);
            if length(est) == K
                err_coh_s(mc) = mean((sort(theta_true) - sort(est)).^2);
                ok_coh_s(mc)  = true;
            end
        catch; end
    end

    if any(ok_incoh_s), rmse_incoh_snap(ni) = sqrt(mean(err_incoh_s(ok_incoh_s))); end
    if any(ok_coh_s),   rmse_coh_snap(ni)   = sqrt(mean(err_coh_s(ok_coh_s)));     end
    fprintf('N=%4d  非相干RMSE=%.3f°  相干RMSE=%.3f°\n', N_s, ...
        rmse_incoh_snap(ni), rmse_coh_snap(ni));
end

%% ================== 绘图辅助函数 ==================
function label_points(xdata, ydata, color, ax)
    for i = 1:length(xdata)
        if ~isnan(ydata(i))
            text(ax, xdata(i), ydata(i)*1.08, ...
                sprintf('%.3f', ydata(i)), ...
                'FontSize', 6.5, 'Color', color, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'bottom');
        end
    end
end

%% ================== 绘图 ==================
c_blue = [0.15 0.45 0.78];
c_red  = [0.80 0.10 0.10];

%% 图1：RMSE vs SNR（semilogy，含数值标注）
figure('Color','w','Name','RMSE vs SNR');
ax1 = axes;
semilogy(SNR_vec, rmse_incoh, '-o', 'Color',c_blue, 'LineWidth',2, ...
    'MarkerFaceColor',c_blue, 'MarkerSize',7); hold on;
semilogy(SNR_vec, rmse_coh,   '--^', 'Color',c_red,  'LineWidth',2, ...
    'MarkerFaceColor',c_red,  'MarkerSize',7);
label_points(SNR_vec, rmse_incoh, c_blue, ax1);
label_points(SNR_vec, rmse_coh,   c_red,  ax1);
grid on; box on;
ylim([1e-2, 1e2]); xlim([SNR_vec(1), SNR_vec(end)]);
xlabel('SNR (dB)', 'FontSize',11);
ylabel('RMSE (°)', 'FontSize',11);
title('相干源信号 vs 非相干源信号', 'FontSize',12);
legend('非相干源', sprintf('相干源(\\rho = %.1f)', rho), ...
    'Location','northeast', 'FontSize',10);
text(0.04, 0.18, ...
    sprintf('NC-SBL\nN=%d, M=%d, K=%d', N_snap, M_arr, K), ...
    'Units','normalized','FontSize',9, ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[1.0 0.97 0.97], 'EdgeColor',[0.8 0.3 0.3], ...
    'LineWidth',1.0);

%% 图2：检测成功率 vs SNR（含数值标注）
figure('Color','w','Name','检测成功率 vs SNR');
ax2 = axes;
plot(SNR_vec, rate_incoh*100, '-o', 'Color',c_blue, 'LineWidth',2, ...
    'MarkerFaceColor',c_blue, 'MarkerSize',7); hold on;
plot(SNR_vec, rate_coh*100,   '--^', 'Color',c_red,  'LineWidth',2, ...
    'MarkerFaceColor',c_red,  'MarkerSize',7);
label_points(SNR_vec, rate_incoh*100, c_blue, ax2);
label_points(SNR_vec, rate_coh*100,   c_red,  ax2);
grid on; box on;
ylim([0, 115]); xlim([SNR_vec(1), SNR_vec(end)]);
xlabel('SNR (dB)', 'FontSize',11);
ylabel('正确检测率 (%)', 'FontSize',11);
title('检测成功率 vs SNR', 'FontSize',12);
legend('非相干 BPSK', sprintf('相干 BPSK (\\rho=%.1f)', rho), ...
    'Location','southeast', 'FontSize',10);
yline(90, 'k--', '90%阈值', 'LabelHorizontalAlignment','left', 'FontSize',9);
text(0.04, 0.18, ...
    sprintf('NC-SBL\nN=%d, M=%d, K=%d', N_snap, M_arr, K), ...
    'Units','normalized','FontSize',9, ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[1.0 0.97 0.97], 'EdgeColor',[0.8 0.3 0.3], ...
    'LineWidth',1.0);

%% 图3：RMSE vs 快拍数（semilogx 线性Y轴，含数值标注）
figure('Color','w','Name','RMSE vs 快拍数');
ax3 = axes;
semilogx(Nsnap_vec, rmse_incoh_snap, '-o', 'Color',c_blue, 'LineWidth',2, ...
    'MarkerFaceColor',c_blue, 'MarkerSize',7); hold on;
semilogx(Nsnap_vec, rmse_coh_snap,   '--^', 'Color',c_red,  'LineWidth',2, ...
    'MarkerFaceColor',c_red,  'MarkerSize',7);
label_points(Nsnap_vec, rmse_incoh_snap, c_blue, ax3);
label_points(Nsnap_vec, rmse_coh_snap,   c_red,  ax3);
grid on; box on;
xlim([Nsnap_vec(1), Nsnap_vec(end)]);
xlabel('快拍数', 'FontSize',11);
ylabel('RMSE (°)', 'FontSize',11);
title('RMSE vs 快拍数', 'FontSize',12);
legend('非相干 BPSK', sprintf('相干 BPSK (\\rho=%.1f)', rho), ...
    'Location','northeast', 'FontSize',10);
text(0.04, 0.18, ...
    sprintf('NC-SBL\nSNR=%ddB, M=%d, K=%d', SNR_snap, M_arr, K), ...
    'Units','normalized','FontSize',9, ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[1.0 0.97 0.93], 'EdgeColor',[0.8 0.5 0.2], ...
    'LineWidth',1.0);

fprintf('\n仿真完成。\n');