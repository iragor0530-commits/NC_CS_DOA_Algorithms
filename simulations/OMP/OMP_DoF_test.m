%% --- 基于非圆扩展 OMP 的自由度压力测试 (M = 3, 5, 8, 10, 12, 15) ---
clear; clc;

N = 8;                        % 物理阵元数
M_test_list = [3, 5, 8, 10, 12, 15];
SNR_fixed = 20;               % 固定信噪比 (dB)
K_snapshots = 500;            % 快拍数

% 构建角度搜索网格与字典矩阵
grid_theta = -90:0.1:90;
P = length(grid_theta);
A_dict = exp(1j * pi * (0:N-1)' * sin(deg2rad(grid_theta)));

% 创建图窗
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.05, 0.7, 0.85]);

for i = 1:length(M_test_list)
    M_curr = M_test_list(i);

    % 1. 生成随机角度（最小间隔 4.5°，范围 [-75, 75]）
    theta_curr = sort(generate_angles(M_curr, -65, 65, 12));

    % 2. 生成 BPSK 非圆信号
    S_curr = (randn(M_curr, K_snapshots) > 0) * 2 - 1;
    A_curr = exp(1j * pi * (0:N-1)' * sin(deg2rad(theta_curr)));

    % 3. 添加高斯白噪声
    As_curr = A_curr * S_curr;
    Pn_curr = mean(abs(As_curr(:)).^2) / (10^(SNR_fixed/10));
    V_curr  = sqrt(Pn_curr/2) * (randn(N, K_snapshots) + 1j*randn(N, K_snapshots));
    Y_curr  = As_curr + V_curr;

    % 4. 调用 NC_OMP 算法
    fprintf('正在计算 M = %2d 的空间谱 (OMP)...\n', M_curr);
    [P_db_curr, ~] = NC_OMP(Y_curr, A_dict, grid_theta, M_curr);

    % 5. 绘制子图
    subplot(3, 2, i);
    plot(grid_theta, P_db_curr, 'b', 'LineWidth', 1.2, 'DisplayName', '估计谱');
    hold on;
    stem(theta_curr, zeros(size(theta_curr)), 'r--', 'LineWidth', 1, ...
         'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', '真实值');

    grid on;
    ylim([-40, 5]);
    xlim([-90, 90]);

    % 标题
    if M_curr < N
        title(['M = ', num2str(M_curr), ' (常规工况)']);
    elseif M_curr < 2*N-1
        title(['M = ', num2str(M_curr), ' (突破物理极限 N-1)']);
    else
        title(['M = ', num2str(M_curr), ' (理论上限 2N-1)']);
    end

    if i >= 5, xlabel('入射角 (degree)'); end
    if mod(i,2) ~= 0, ylabel('归一化功率 (dB)'); end
    if i == 1, legend('Location', 'SouthWest', 'FontSize', 7); end
end

sgtitle(['基于非圆信号扩展的 NC-OMP 自由度测试 (物理阵元 N=8, SNR=', ...
         num2str(SNR_fixed), 'dB)'], 'FontSize', 14);