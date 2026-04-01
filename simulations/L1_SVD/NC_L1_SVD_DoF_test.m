%% --- 非圆信号自由度压力测试 (M = 3, 5, 8, 10, 12, 15) ---
M_test_list = [3, 5, 8, 10, 12, 15]; 
SNR_fixed = 20;               % 较高信噪比以观察算法潜力
K_snapshots = 100;

% 创建大图窗口，调整尺寸以适应 6 个子图
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.05, 0.7, 0.85]);

for i = 1:length(M_test_list)
    M_curr = M_test_list(i);
    
    % 1. 生成当前信源数下的信号环境
    % 确保角度在 [-75, 75] 之间，且最小间隔为 4.5 度以保证基本可分性
    theta_curr = sort(generate_angles(M_curr, -75, 75, 4.5)); 
    
    % 生成 BPSK 非圆信号
    S_curr = (randn(M_curr, K_snapshots) > 0) * 2 - 1; 
    A_curr = exp(1j * pi * (0:N-1)' * sin(deg2rad(theta_curr)));
    
    % 添加高斯白噪声
    As_curr = A_curr * S_curr;
    Pn_curr = mean(abs(As_curr(:)).^2) / (10^(SNR_fixed/10));
    V_curr = sqrt(Pn_curr/2) * (randn(N, K_snapshots) + 1j*randn(N, K_snapshots));
    Y_curr = As_curr + V_curr;
    
    % 2. 调用 NC_L1_SVD 算法进行重构
    fprintf('正在计算 M = %2d 的空间谱...\n', M_curr);
    [P_db_curr, ~] = NC_L1_SVD(Y_curr, A_dict, grid_theta, M_curr);
    
    % 3. 绘制子图 (3行2列布局)
    subplot(3, 2, i);
    plot(grid_theta, P_db_curr, 'b', 'LineWidth', 1.2, 'DisplayName', '估计谱'); 
    hold on;
    
    % 绘制真实角度参考 (用红色虚线和圆圈标注)
    stem(theta_curr, zeros(size(theta_curr)), 'r--', 'LineWidth', 1, ...
         'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', '真实值');
    
    % 坐标轴与标题美化
    grid on;
    ylim([-40, 5]); 
    xlim([-90, 90]);
    
    % 标题显示关键信息：是否突破物理阵元数 N=8
    if M_curr < N
        title(['M = ', num2str(M_curr), ' (常规工况)']);
    elseif M_curr < 15
        title(['M = ', num2str(M_curr), ' (突破物理极限 N-1)']);
    else
        title(['M = ', num2str(M_curr), ' (理论上限 2N-1)']);
    end
    
    % 只有底部的图显示横坐标标签，最左侧显示纵坐标标签
    if i >= 5, xlabel('入射角 (degree)'); end
    if mod(i, 2) ~= 0, ylabel('归一化功率 (dB)'); end
    if i == 1, legend('Location', 'SouthWest', 'FontSize', 7); end
end

sgtitle(['基于非圆信号扩展的 L1-SVD 自由度测试 (物理阵元 N=8, SNR=', num2str(SNR_fixed), 'dB)'], 'FontSize', 14);