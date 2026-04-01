clc; clear;
%% ===== 参数设定 =====
N   = 8;
M   = 5;
K   = 100;
SNR = 10;

x_n = (0:N-1)';

%% ===== 信号生成 =====
theta_true = sort(generate_angles(M, -60, 60, 10));
S          = (randn(M, K) > 0) * 2 - 1;
A          = exp(1j * pi * x_n * sin(deg2rad(theta_true)));
As_k       = A * S;
sigma2     = mean(abs(As_k(:)).^2) / (10^(SNR/10));
V          = sqrt(sigma2/2) * (randn(N,K) + 1j*randn(N,K));
Y          = As_k + V;

%% ===== 角度搜索网格 =====
grid_step  = 0.1;
grid_theta = -90:grid_step:90;

%% ===== 调用 NC-MUSIC =====
fprintf('真实角度: %s 度\n', num2str(theta_true));
[P_db, est_theta] = NC_MUSIC(Y, grid_theta, M);
fprintf('估计角度: %s 度\n', num2str(est_theta));
fprintf('估计误差: %s 度\n', num2str(est_theta - sort(theta_true)));

%% ===== 单次估计绘图 =====
figure('Color','w');
plot(grid_theta, P_db, 'Color',[0.15 0.45 0.78], 'LineWidth', 1.5);
hold on;

% 找每个估计角度对应的谱峰值
est_sorted  = sort(est_theta(:)');
true_sorted = sort(theta_true(:)');

for k = 1:M
    % 估计角度对应的谱值
    [~, idx_est] = min(abs(grid_theta - est_sorted(k)));
    peak_val = P_db(idx_est);

    % 红圈画在峰值处
    plot(est_sorted(k), peak_val, 'ro', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor','r');

    % 估计角度标注（峰值上方）
    text(est_sorted(k), peak_val + 1.5, sprintf('%.1f°', est_sorted(k)), ...
        'HorizontalAlignment','center', 'FontSize', 8, ...
        'Color',[0.15 0.45 0.78], 'FontWeight','bold');

    % 真实角度标注（红圈下方）
    text(est_sorted(k), peak_val - 3.5, sprintf('真实:%.1f°', true_sorted(k)), ...
        'HorizontalAlignment','center', 'FontSize', 7.5, 'Color','r');
end

grid on;
ylim([-60, 5]);
xlabel('入射角 (°)', 'FontSize', 12);
ylabel('归一化空间谱 (dB)', 'FontSize', 12);
legend('NC-MUSIC 空间谱', '估计峰值', 'Location','southeast', 'FontSize',10);

% 参数文字框
text(0.02, 0.97, ...
    sprintf('阵元数 N = %d\n信源数 M = %d\n快拍数 K = %d\n信噪比 SNR = %d dB', N, M, K, SNR), ...
    'Units','normalized', 'FontSize', 9, ...
    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
    'BackgroundColor',[0.95 0.97 1.0], 'EdgeColor',[0.4 0.4 0.8], 'LineWidth', 1.0);

%% ===== 蒙特卡洛仿真 =====
snr_range       = -20:2:10;
monte_carlo_num = 500;

theta_fixed = sort(generate_angles(M, -50, 50, 10));

if isempty(gcp('nocreate'))
    parpool;
end

rmse_snr = zeros(length(snr_range), 1);
crb_snr  = zeros(length(snr_range), 1);

fprintf('\n开始蒙特卡洛仿真 (%d SNR点 × %d次)...\n', length(snr_range), monte_carlo_num);

for s = 1:length(snr_range)
    cur_SNR = snr_range(s);
    crb_vec    = compute_NC_CRB_internal(theta_fixed, N, cur_SNR, K);
    crb_snr(s) = mean(crb_vec);
    all_sq_errors = NaN(monte_carlo_num, M);

    parfor mc = 1:monte_carlo_num
        theta_mc  = theta_fixed;
        S_mc      = (randn(M, K) > 0) * 2 - 1;
        A_mc      = exp(1j * pi * (0:N-1)' * sin(deg2rad(theta_mc)));
        As_mc     = A_mc * S_mc;
        sigma2_mc = mean(abs(As_mc(:)).^2) / (10^(cur_SNR/10));
        V_mc      = sqrt(sigma2_mc/2) * (randn(N,K) + 1j*randn(N,K));
        Y_mc      = As_mc + V_mc;

        [~, est_mc] = NC_MUSIC(Y_mc, grid_theta, M);
        est_sorted  = sort(est_mc(:));
        true_sorted = sort(theta_mc(:));
        if length(est_sorted) == M
            all_sq_errors(mc, :) = (est_sorted - true_sorted).^2;
        else
            all_sq_errors(mc, :) = 90^2 * ones(1, M);
        end
    end

    rmse_snr(s) = sqrt(mean(all_sq_errors(:)));
    fprintf('SNR = %+3d dB  -->  RMSE = %.4f °,  CRB = %.4f °\n', ...
        cur_SNR, rmse_snr(s), crb_snr(s));
end

%% ===== 性能曲线绘图 =====
figure('Color','w','Position',[100 100 750 460]);

c_rmse = [0.15 0.45 0.78];
c_crb  = [0.80 0.10 0.10];

semilogy(snr_range, rmse_snr, '-s', 'Color',c_rmse, 'LineWidth',2, ...
    'MarkerSize',7, 'MarkerFaceColor',c_rmse, 'DisplayName','NC-MUSIC RMSE');
hold on;
semilogy(snr_range, crb_snr, '--^', 'Color',c_crb, 'LineWidth',2, ...
    'MarkerSize',7, 'MarkerFaceColor',c_crb, 'DisplayName','NC-CRB 理论下界');

% 每点数值标注
for i = 1:length(snr_range)
    text(snr_range(i), rmse_snr(i)*1.20, sprintf('%.3f', rmse_snr(i)), ...
        'FontSize',6.5, 'Color',c_rmse, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    text(snr_range(i), crb_snr(i)*0.75, sprintf('%.3f', crb_snr(i)), ...
        'FontSize',6.5, 'Color',c_crb, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top');
end

grid on;
set(gca, 'YMinorGrid','on');
xlabel('信噪比 SNR (dB)',          'FontSize',12, 'FontWeight','bold');
ylabel('均方根误差 RMSE (°)',      'FontSize',12, 'FontWeight','bold');
legend('Location','northeast',     'FontSize',11);
ylim([1e-2, 1e2]);
xticks(snr_range);

% 参数文字框
text(0.02, 0.18, ...
    sprintf('阵元数 N = %d\n信源数 M = %d\n快拍数 K = %d', N, M, K), ...
    'Units','normalized', 'FontSize',9, ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
    'BackgroundColor',[0.95 0.97 1.0], 'EdgeColor',[0.4 0.4 0.8], 'LineWidth',1.0);