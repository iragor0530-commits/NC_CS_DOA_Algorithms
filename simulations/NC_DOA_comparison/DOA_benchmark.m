%% NC_DOA_benchmark.m
%  五种 NC 算法 RMSE vs SNR 对比 + 确定性模型 NC-CRB 下界
%  信源数 K=3，阵元数 M=8，快拍数 N=200，MC=500

clear; clc; close all;

if isempty(gcp('nocreate'))
    parpool('local', 4);
end

%% ================== 公共参数 ==================
M_arr  = 8;
K      = 3;
N_snap = 20;
n_MC   = 500;
SNR_vec = -20:4:10;
n_snr   = length(SNR_vec);

%% 各算法使用各自合适的网格
% NC-SBL：1° 足够，算法自带 off-grid 细化
grid_SBL   = -90:1:90;

% NC-MUSIC：谱峰搜索
grid_MUSIC = -90:0.5:90;

% NC-OMP / NC-L1-SVD：字典直接决定精度，用 0.1°
grid_sparse = -90:0.1:90;

% NC-ESPRIT：无需网格，直接解析求解

%% 各算法字典矩阵（仅稀疏类需要）
A_dict_sparse = exp(1j*pi*(0:M_arr-1).' * sind(grid_sparse));  % M×361，OMP/L1-SVD用

%% 预生成角度池（固定种子，可重复）
rng(42);
theta_pool = zeros(n_MC, K);
for mc = 1:n_MC
    theta_pool(mc, :) = generate_angles(K, -60, 60, 10);
end
theta_rep = mean(theta_pool, 1);   % 代表性角度，用于 CRB

%% ================== CRB 计算 ==================
crb_curve = zeros(1, n_snr);
for si = 1:n_snr
    crb_vec       = compute_NC_CRB_internal(theta_rep, M_arr, SNR_vec(si), N_snap);
    crb_curve(si) = sqrt(mean(crb_vec.^2));
end

%% ================== 算法仿真 ==================
methods   = {'NC-SBL','NC-MUSIC','NC-OMP','NC-ESPRIT','NC-L1-SVD'};
n_methods = length(methods);
RMSE      = nan(n_methods, n_snr);
RATE      = zeros(n_methods, n_snr);

for si = 1:n_snr
    snr_db = SNR_vec(si);
    sigma2 = 10^(-snr_db/10);

    % 修复1：预分配改为 n_MC × n_methods，parfor 按行切片更安全
    err = nan(n_MC, n_methods);
    ok  = false(n_MC, n_methods);

    parfor mc = 1:n_MC
        theta_true = theta_pool(mc, :);
        A     = exp(1j*pi*(0:M_arr-1).' * sind(theta_true));
        S     = sign(randn(K, N_snap));
        noise = sqrt(sigma2/2)*(randn(M_arr,N_snap)+1j*randn(M_arr,N_snap));
        Y     = A*S + noise;

        err_mc = nan(1, n_methods);
        ok_mc  = false(1, n_methods);

        t_sorted = sort(theta_true(:));   % 预排序真值，避免重复调用

        % 1. NC-SBL（1° 网格，自带细化，取前K个最强峰）
        try
            est = NC_SBL(Y, grid_SBL);
            est = sort(est(1:min(numel(est), K)));   % 修复：取前K个，不要求恰好等于K
            if numel(est) == K && all(isfinite(est))
                err_mc(1) = mean((t_sorted - est(:)).^2);
                ok_mc(1)  = true;
            end
        catch; end

        % 2. NC-MUSIC（0.1° 扫描网格）
        try
            [~, est] = NC_MUSIC(Y, grid_MUSIC, K);
            est = sort(est(:));
            if numel(est) == K && all(isfinite(est))
                err_mc(2) = mean((t_sorted - est).^2);
                ok_mc(2)  = true;
            end
        catch; end

        % 3. NC-OMP（0.5° 字典）
        try
            [~, est] = NC_OMP(Y, A_dict_sparse, grid_sparse, K);
            est = sort(est(:));
            if numel(est) == K && all(isfinite(est))
                err_mc(3) = mean((t_sorted - est).^2);
                ok_mc(3)  = true;
            end
        catch; end

        % 4. NC-ESPRIT（无需网格）
        try
            est = sort(NC_ESPRIT(Y, K));
            est = est(:);
            if numel(est) == K && all(isfinite(est))   % 修复：isfinite 替代仅长度检查
                err_mc(4) = mean((t_sorted - est).^2);
                ok_mc(4)  = true;
            end
        catch; end

        % 5. NC-L1-SVD（0.5° 字典）
        try
            [~, est] = NC_L1_SVD(Y, A_dict_sparse, grid_sparse, K, sigma2, []);
            est = sort(est(:));
            if numel(est) == K && all(isfinite(est))
                err_mc(5) = mean((t_sorted - est).^2);
                ok_mc(5)  = true;
            end
        catch; end

        % 修复2：按行赋值，parfor sliced variable 安全
        err(mc, :) = err_mc;
        ok(mc, :)  = ok_mc;
    end

    % 修复3：汇总时按列索引（对应新的 n_MC × n_methods 布局）
    for m = 1:n_methods
        valid = ok(:, m);
        RATE(m, si) = mean(valid);
        if any(valid)
            RMSE(m, si) = sqrt(mean(err(valid, m)));
        end
    end

    fprintf('SNR=%+3ddB  ', snr_db);
    for m = 1:n_methods
        fprintf('%s %.0f%%  ', methods{m}, RATE(m,si)*100);
    end
    fprintf('\n');
end

%% ================== 绘图 ==================
colors  = [0.15 0.45 0.78;
           0.80 0.10 0.10;
           0.13 0.70 0.35;
           0.85 0.55 0.00;
           0.55 0.20 0.75];
markers = {'o','^','s','d','v'};
lines   = {'-','--','-.',':','-'};

figure('Color','w','Position',[60 60 1300 520]);

%% 子图1：RMSE vs SNR
for m = 1:n_methods
    semilogy(SNR_vec, RMSE(m,:), ...
        [lines{m}, markers{m}], ...
        'Color', colors(m,:), 'LineWidth', 2, ...
        'MarkerFaceColor', colors(m,:), 'MarkerSize', 7, ...
        'DisplayName', methods{m});
    hold on;
end
semilogy(SNR_vec, crb_curve, 'k-', 'LineWidth', 2.5, 'DisplayName', 'NC-CRB');

grid on; box on;
ylim([1e-2, 1e2]); xlim([SNR_vec(1), SNR_vec(end)]);
xlabel('信噪比 (dB)', 'FontSize',11);
ylabel('RMSE (°)',    'FontSize',11);
title('RMSE vs SNR（含 NC-CRB 下界）', 'FontSize',12);
legend('Location','southwest', 'FontSize',9, 'NumColumns',2);
text(0.02, 0.98, ...
    sprintf('M=%d, K=%d, N=%d, MC=%d\n随机角度 ±60°, 最小间隔 10°\nSBL:1° / MUSIC:0.1° / OMP&L1:0.5°', ...
    M_arr, K, N_snap, n_MC), ...
    'Units','normalized','VerticalAlignment','top', ...
    'FontSize',8.5,'BackgroundColor',[0.97 0.97 0.97],'EdgeColor',[0.6 0.6 0.6]);


fprintf('\n仿真完成。\n');