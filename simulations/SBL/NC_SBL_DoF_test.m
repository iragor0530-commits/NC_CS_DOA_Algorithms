%% =========================================================
%  SBL_test.m
%  NC-SBL 自由度压力测试：不同信源数 K 对估计性能的影响
%  依赖：NC_SBL.m、generate_angles.m
% =========================================================

%% 1. 参数配置
M          = 8;
N          = 500;
SNR_dB     = 10;
K_max      = 15;
min_sep    = 8;
grid_step  = 1;
grid_theta = -90:grid_step:90;

    %% 2. 压力测试主循环
    results = struct();

for K = 1:K_max

    % 2.1 随机生成真实 DOA（off-grid，不对齐网格）
    true_theta = sort(generate_angles(K, -60, 60, min_sep));
    true_theta = true_theta(:)';   % 强制行向量
    K_act      = length(true_theta);

    % 2.2 构造接收信号
    A       = exp(1j*pi*(0:M-1)' * sind(true_theta));  % [M × K_act]
    S       = 2*randi([0,1], K_act, N) - 1;            % [K_act × N]
    sig_pow = mean(abs(A*S).^2, 'all');
    noi_pow = sig_pow / 10^(SNR_dB/10);
    noise   = sqrt(noi_pow/2) * (randn(M,N) + 1j*randn(M,N));
    Y       = A*S + noise;

    % 2.3 调用 NC_SBL
    est_theta = NC_SBL(Y, grid_theta);

    % 2.4 存储结果
    results(K).true_theta = true_theta;
    results(K).est_theta  = est_theta;
    results(K).K_act      = K;
    results(K).success    = (length(est_theta) == K);
    if results(K).success
        results(K).rmse = sqrt(mean((sort(true_theta)-sort(est_theta)).^2));
    else
        results(K).rmse = NaN;
    end
end

%% 3. 茎图逐K对比
rows = ceil(K_max / 3);
cols = 3;

figure('Name','NC-SBL 自由度压力测试','NumberTitle','off', ...
       'Position',[100 50 1200 900]);

for K = 1:K_max
    true_theta = results(K).true_theta;
    est_theta  = results(K).est_theta;
    success    = results(K).success;

    subplot(rows, cols, K);
    plot(true_theta, ones(1,K), 'ro', ...
        'MarkerSize',7,'MarkerFaceColor','r','LineStyle','none');
    hold on;
    stem(est_theta, ones(1,length(est_theta)), 'b', ...
        'LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','b');

    % ★ 修正：强制用 %g 避免科学计数法
    for k = 1:K
        text(true_theta(k), 1.12, sprintf('%.1f°', true_theta(k)), ...
            'Color','r','FontSize',7,'HorizontalAlignment','center');
    end
    for k = 1:length(est_theta)
        text(est_theta(k), 0.62, sprintf('%.1f°', est_theta(k)), ...
            'Color','b','FontSize',7,'HorizontalAlignment','center');
    end

    xlim([-90 90]); ylim([0 1.3]); yticks([]);
    xlabel('角度 (°)','FontSize',9);

    if success
        title_str   = sprintf('K=%d | 估计K=%d | RMSE=%.2f°  ✓', ...
            K, length(est_theta), results(K).rmse);
        title_color = [0 0.5 0];
    else
        title_str   = sprintf('K=%d | 估计K=%d  ✗', K, length(est_theta));
        title_color = [0.8 0 0];
    end
    title(title_str,'FontSize',9,'Color',title_color);
    grid on; hold off;
end

sgtitle(sprintf('NC-SBL DOA 自由度压力测试  |  SNR=%ddB  |  M=%d  |  N=%d', ...
    SNR_dB, M, N), 'FontSize',13,'FontWeight','bold');