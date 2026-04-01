%% ================== NC-ESPRIT 自由度压力测试 ==================
%% 固定 M=8，逐步增加 K，找崩溃点
clear; clc; close all;

%% ---------- 固定参数 ----------
M       = 8;
N_snap  = 100;
SNR     = 10;
K_list  = 1:15;          % NC-ESPRIT 理论极限 2M-1 = 15

angle_min_sep = 5;       % 最小角度间隔（度），防止信源过近

%% ---------- 逐K测试 ----------
n_cols = 3;
n_rows = ceil(length(K_list) / n_cols);
figure('Color','w','Position',[50 50 1400 900]);

for ki = 1:length(K_list)
    K = K_list(ki);

    %% 随机生成角度（保证最小间隔）
    max_try = 1000;
    success_gen = false;
    for t = 1:max_try
        theta_true = sort(-70 + 140 * rand(1, K));
        if K == 1 || min(diff(theta_true)) >= angle_min_sep
            success_gen = true;
            break;
        end
    end
    if ~success_gen
        fprintf('K=%d: 无法生成满足间隔要求的角度\n', K);
        continue;
    end

    %% 构造信号
    theta_rad = deg2rad(theta_true);
    A = exp(1j * pi * (0:M-1).' * sin(theta_rad));
    S = randn(K, N_snap);
    sigma2 = 10^(-SNR/10);
    noise = sqrt(sigma2/2) * (randn(M, N_snap) + 1j*randn(M, N_snap));
    Y = A * S + noise;

    %% NC-ESPRIT 估计
    try
        est_theta = NC_ESPRIT(Y, K);
        % 匹配真实角度计算误差
        theta_sorted = sort(theta_true);
        err = est_theta(:)' - theta_sorted;
        rmse = sqrt(mean(err.^2));
        if rmse < 3
            status = 'OK';
        else
            status = 'WARN';        % 估计偏差过大
        end
    catch ME
        est_theta = nan(K, 1);
        rmse = nan;
        status = 'FAIL';
    end

    %% 子图绘制
subplot(n_rows, n_cols, ki);

% 估计角度用竖线
if ~strcmp(status,'FAIL')
    stem(est_theta, 0.75*ones(K,1), 'b', 'LineWidth', 1.8, ...
        'MarkerFaceColor','b', 'MarkerSize', 5);
    hold on;
    % 估计角度数值标注
    for k = 1:K
        text(est_theta(k), 0.82, sprintf('%.1f°', est_theta(k)), ...
            'HorizontalAlignment','center', 'FontSize', 6, 'Color','b');
    end
end
hold on;

% 真实角度用 o 散点
plot(sort(theta_true), 0.3*ones(1,K), 'ro', ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'MarkerFaceColor','r');

% 真实角度数值标注
for k = 1:K
    text(sort(theta_true(k)), 0.37, sprintf('%.1f°', sort(theta_true(k))), ...
        'HorizontalAlignment','center', 'FontSize', 6, 'Color','r');
end

xlim([-90 90]); ylim([0 1.3]);
grid on;

% 标题颜色反映状态
switch status
    case 'OK'
        color = [0.0 0.55 0.27];
        tag   = sprintf('K=%d ', K);
    case 'WARN'
        color = [0.85 0.55 0.0];
        tag   = sprintf('K=%d ', K);
    case 'FAIL'
        color = [0.8 0.0 0.0];
        tag   = sprintf('K=%d ', K);
end
title(tag, 'Color', color, 'FontSize', 9, 'FontWeight','bold');

if mod(ki-1, n_cols) == 0
    ylabel('幅度','FontSize',8);
end
if ki > (n_rows-1)*n_cols
    xlabel('角度（°）','FontSize',8);
end

% 图例只在第一个子图显示
if ki == 1
    legend('估计角度','真实角度','Location','north','FontSize',7);
end

fprintf('K=%2d 状态: %-4s ', K, status);
if ~isnan(rmse), fprintf('%.3f°\n', rmse);
else,            fprintf('N/A\n'); end
end
