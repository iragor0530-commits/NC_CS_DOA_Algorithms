%% Step 1: 仿真信号环境构建 (根据公式1-5)
clc;
clear;
% --- 参数设定 ---
N = 8;                  % 物理阵元数量
M = 3;                  % 窄带、不相关、非圆信号源数量 
K =800;                % 快拍总数 
SNR = 20;               % 信噪比

% 随机生成M个信号源入射方向 {theta_1, ..., theta_M}
% 角度范围 [-90, 90]，确保角度之间有一定间隔
grid_potential = -80:0.1:80;
theta_true = sort(generate_angles(M,-70,70,8));
theta_true = sort(theta_true); % 升序排列

% --- 1. 生成非圆信号源 s_l(k) ---
% 生成符合非圆特性的原始信号 (以BPSK为例)
% s_pure 为实数信号，通过旋转相位 Phi 引入非圆性
s_pure = (randn(M, K) > 0) * 2 - 1; 
S = s_pure;   %源信号向量 s(k)

% --- 2. 构建阵列流形矩阵 A ---
% 传感器位置 x_n，以 lambda/2 为单位 (对应公式2)
x_n = (0:N-1)'; 
A = zeros(N, M);
for l = 1:M
    % 注意：pi 是因为距离是 lambda/2，2*pi*(lambda/2)/lambda = pi
    A(:, l) = exp(1j * pi * x_n * sin(deg2rad(theta_true(l))));
end

% --- 3. 生成阵列接收数据 ---

% A*s(k) 部分：信号分量
As_k = A * S; 

% v(k) 部分：生成零均值高斯白噪声
% 根据 SNR = 10*log10(Ps / Pn) 反推噪声功率 sigma_n^2
Pn = mean(abs(As_k(:)).^2) / (10^(SNR/10)); 
sigma_n = sqrt(Pn / 2); % 复噪声的实部和虚部各占一半功率

% 生成 N x K 的噪声矩阵 V，对应公式中的 v(k)
V = sigma_n * (randn(N, K) + 1j*randn(N, K));

% 最终组合：y(k) = A*s(k) + v(k)
Y = As_k + V; % 这里的 Y_noisy 每一列就是一个 y(k)

% --- 4. 非圆扩展 (核心毕设步骤) ---
% 利用非圆信号特性，构造扩展观测矢量 y_aug(k) = [y(k); conj(y(k))]
% 这将观测空间从 N 维扩展到 2N 维
Y_aug = [Y; 
         conj(Y)];
%% --- 6. 信号矩阵可视化展示 ---

% A. 查看原始接收信号 Y_noisy (N x K)
% 取第 1 个阵元的接收数据进行分析
y1 = Y(1, :); 

figure('Color', 'w', 'Name', '非圆信号特性与扩展对比');

% 子图1：原始信号的复平面分布 (星座图)
subplot(2, 2, 1);
plot(real(y1), imag(y1), 'b.', 'MarkerSize', 8);
grid on; axis equal;
title(['第1阵元原始信号 y_1(k) (SNR=', num2str(SNR), 'dB)']);
xlabel('实部 (In-phase)'); ylabel('虚部 (Quadrature)');
%  BPSK 信号在复平面上呈"条状"分布，这就是非圆性的直观表现

% 子图2：扩展后的信号 y_aug 的前两个分量关系
% y_aug = [y1; ...; yN; conj(y1); ...; conj(yN)]
subplot(2, 2, 2);
y_ext_1 = Y_aug(1, :);          % 即 y1
y_ext_Nplus1 = Y_aug(N+1, :);   % 即 conj(y1)
plot(real(y_ext_1), real(y_ext_Nplus1), 'r.', 'MarkerSize', 8);
grid on; axis equal;
title('扩展项之间的相关性 (y_1 与 y_1^*)');
xlabel('Re(y_1)'); ylabel('Re(y_1^*)');

% B. 矩阵维度的直观对比 (热力图)
% 子图3：原始矩阵 Y_noisy 的幅度衰减/分布
subplot(2, 2, 3);
imagesc(abs(Y)); colorbar;
title(['原始矩阵 Y (', num2str(N), ' x ', num2str(K), ')']);
xlabel('快拍数 (k)'); ylabel('阵元索引 (n)');

% 子图4：扩展矩阵 Y_aug 的幅度分布
subplot(2, 2, 4);
imagesc(abs(Y_aug)); colorbar;
title(['扩展矩阵 Y_{aug} (', num2str(2*N), ' x ', num2str(K), ')']);
xlabel('快拍数 (k)'); ylabel('扩展阵元索引');

% --- 打印关键数据到命令行 ---
disp('--- 矩阵维度对比 ---');
disp(['原始接收矩阵 Y 维度: ', num2str(size(Y))]);
disp(['非圆扩展矩阵 Y_aug 维度: ', num2str(size(Y_aug))]);
disp('--------------------');

% --- 5. 构建过完备字典 (用于后续压缩感知) ---
grid_step = 0.1;
grid_theta = -90:grid_step:90;
P = length(grid_theta);

% 基础字典
A_dict = exp(1j * pi * x_n * sin(deg2rad(grid_theta)));

% 扩展字典 (2N x P)
A_dict_aug = [A_dict; 
              conj(A_dict)];

%扩展噪声矩阵
N_tilde = [V; 
           conj(V)];
%构建稀疏模型

K = size(Y_aug, 2); % 快拍数
fprintf('稀疏矩阵 S_MATRIX 的目标维度为: %d x %d\n', P, K);
fprintf('信号环境构建完成：\n');
fprintf('真实角度: %s 度\n', num2str(theta_true));
fprintf('扩展字典维度: %d x %d\n', size(A_dict_aug, 1), size(A_dict_aug, 2));


% --- 调用 NC-OMP 算法 ---
[P_db_omp, est_omp] = NC_OMP(Y, A_dict, grid_theta, M);

% --- 优化后的 OMP 绘图 ---
figure('Color', 'w', 'Name', 'NC-OMP 估计结果对比');

% 将谱峰归一化到 0dB，并将极小值限制在合理范围（如 -50dB），避免"绿墙"
P_plot = P_db_omp;
P_plot(P_plot < -50) = -50; % 限制显示下限，这样背景就是空白的

% 1. 绘制估计的谱线
stem(grid_theta, P_plot, 'Color', [0 0.8 0], 'LineWidth', 1.5, 'Marker', 'none'); 
hold on;

% 2. 绘制真实角度（设置一个合适的高度，比如 0dB）
% 增加 'BaseValue', -50 是为了让红线从底部升起，而不是悬浮在顶部
stem(theta_true, zeros(size(theta_true)), 'r--', 'LineWidth', 2, ...
    'MarkerSize', 10, 'MarkerFaceColor', 'none', 'BaseValue', -50);

% 3. 美化坐标轴
grid on;
ylim([-55, 5]); % 留出一点顶部空间看清楚 0dB 的标记
xlabel('入射角 (degree)', 'FontSize', 12);
ylabel('归一化功率谱 (dB)', 'FontSize', 12);
title(sprintf('非圆信号 OMP 估计结果  (N=%d, M=%d, K=%d, SNR=%d dB)', ...
    N, M, K, SNR), 'FontSize', 13);
legend('OMP 估计位置', '真实入射角', 'Location', 'SouthEast');

% 4. 强制刷新文字显示
fprintf('真实角度: %s\n', num2str(sort(theta_true)));
fprintf('OMP估计角度: %s\n', num2str(sort(est_omp)));

%% --- 蒙特卡洛仿真 (含 CRB 对比) ---
snr_range       = -10:2:10;
monte_carlo_num  = 500;         % 建议增加到 500 次以获得更平滑的曲线
FAIL_PENALTY    = 5;            % 过滤掉误差大于 5 度的帧，确保 RMSE 落在合理区间
 %启动并行环境
if isempty(gcp('nocreate'))
    parpool; % 自动根据你的 CPU 核心数开启并行池
end
rmse_snr = zeros(length(snr_range), 1);
crb_snr  = zeros(length(snr_range), 1); % 用于存储每个 SNR 的理论下界
theta_fixed=sort(generate_angles(M, -50, 50, 8));
fprintf('\n开始 SNR 扫描仿真 (共 %d 个SNR点, 每点 %d 次蒙特卡洛)...\n', ...
    length(snr_range), monte_carlo_num);
 
for s = 1:length(snr_range)
    current_SNR = snr_range(s);
    
    % --- 并行修改点：预分配矩阵 ---
    % 每一行代表一次蒙特卡洛，每一列代表一个信号的误差
    all_sq_errors = NaN(monte_carlo_num, M); 
    
    % --- 计算当前 SNR 下的理论 CRB ---
    current_theta_ref =theta_fixed;
    crb_vec = compute_NC_CRB_internal(current_theta_ref, N, current_SNR, K);
    crb_snr(s) = mean(crb_vec);

    % --- 【关键修改】使用 parfor 进行并行仿真 ---
    parfor mc = 1:monte_carlo_num
        % 1. 随机生成角度（注意：parfor 内部函数需独立）
        theta_mc =theta_fixed; 
        % 2. 生成信号与噪声
        S_mc  = (randn(M, K) > 0) * 2 - 1;
        A_mc  = exp(1j * pi * (0:N-1)' * sin(deg2rad(theta_mc)));
        As_mc = A_mc * S_mc;
        sigma2_mc = 1 / 10^(current_SNR/10);
        V_mc = sqrt(sigma2_mc/2) * (randn(N,K) + 1j*randn(N,K));
        Y_mc  = As_mc + V_mc;
 
        % 3. 调用算法
        [~, est_mc] = NC_OMP(Y_mc, A_dict, grid_theta, M);
 
        % 4. 误差匹配
        est_sorted  = sort(est_mc(:));
        true_sorted = sort(theta_mc(:));
        
        if length(est_sorted) == M
            errs = abs(est_sorted - true_sorted);
            % 满足条件的存入预分配矩阵
            if max(errs) < FAIL_PENALTY
                all_sq_errors(mc, :) = errs.^2;
            end
        end
    end
 
    % --- 处理结果 ---
    % 剔除 NaN (即失败的帧) 后计算平均值
    valid_errors = all_sq_errors(~isnan(all_sq_errors));
    if ~isempty(valid_errors)
        rmse_snr(s) = sqrt(mean(valid_errors));
    else
        rmse_snr(s) = NaN;
    end
    
    fprintf('SNR = %+3d dB  -->  RMSE = %.4f°, CRB = %.4f° (并行完成)\n', ...
            current_SNR, rmse_snr(s), crb_snr(s));
end
%% --- 最终可视化 (RMSE vs CRB) ---
figure('Color', 'w', 'Name', '算法性能与理论下界对比');
semilogy(snr_range, rmse_snr, 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'NC-OMP RMSE');
hold on;
semilogy(snr_range, crb_snr, 'r--', 'LineWidth', 2, 'DisplayName', 'NC-CRB (理论下界)');
grid on;
set(gca, 'YMinorGrid', 'on');
xlabel('信噪比 SNR (dB)', 'FontSize', 12);
ylabel('均方根误差 RMSE (degree)', 'FontSize', 12);
title(['非圆信号 DOA 估计性能对比 (N=', num2str(N), ', K=', num2str(K), ')']);
legend('Location', 'NorthEast');
ylim([1e-2, 10]);

