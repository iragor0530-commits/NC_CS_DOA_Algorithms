%% DOA_benchmark_snapshots.m
%  五种 NC 算法 RMSE vs 快拍数 对比 + NC-CRB 下界
%  信源数 K=3，阵元数 M=8，固定 SNR=0dB，MC=500

clear; clc; close all;

if isempty(gcp('nocreate'))
    parpool('local', 4);
end

%% ================== 公共参数 ==================
M_arr    = 8;
K        = 3;
SNR_dB   = 0;          % 固定信噪比
sigma2   = 10^(-SNR_dB/10);
n_MC     = 500;

% 快拍数扫描
Nsnap_vec = [1,3,5,10, 20, 50, 100, 200, 500, 1000];
n_nsnap   = length(Nsnap_vec);

%% ================== 网格设置 ==================
grid_SBL    = -90:1:90;
grid_MUSIC  = -90:0.5:90;
grid_sparse = -90:0.1:90;

A_dict_sparse = exp(1j*pi*(0:M_arr-1).' * sind(grid_sparse));

%% ================== 预生成角度池 ==================
rng(42);
theta_pool = zeros(n_MC, K);
for mc = 1:n_MC
    theta_pool(mc, :) = generate_angles(K, -60, 60, 10);
end
theta_rep = mean(theta_pool, 1);

%% ================== CRB 计算 ==================
crb_curve = zeros(1, n_nsnap);
for ni = 1:n_nsnap
    crb_vec       = compute_NC_CRB_internal(theta_rep, M_arr, SNR_dB, Nsnap_vec(ni));
    crb_curve(ni) = sqrt(mean(crb_vec.^2));
end

%% ================== 算法仿真 ==================
methods   = {'NC-SBL','NC-MUSIC','NC-OMP','NC-ESPRIT','NC-L1-SVD'};
n_methods = length(methods);
RMSE      = nan(n_methods, n_nsnap);
RATE      = zeros(n_methods, n_nsnap);

for ni = 1:n_nsnap
    N_snap = Nsnap_vec(ni);

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
        t_sorted = sort(theta_true(:));

        % 1. NC-SBL
        try
            est = NC_SBL(Y, grid_SBL);
            est = sort(est(1:min(numel(est), K)));
            if numel(est) == K && all(isfinite(est))
                err_mc(1) = mean((t_sorted - est(:)).^2);
                ok_mc(1)  = true;
            end
        catch; end

        % 2. NC-MUSIC
        try
            [~, est] = NC_MUSIC(Y, grid_MUSIC, K);
            est = sort(est(:));
            if numel(est) == K && all(isfinite(est))
                err_mc(2) = mean((t_sorted - est).^2);
                ok_mc(2)  = true;
            end
        catch; end

        % 3. NC-OMP
        try
            [~, est] = NC_OMP(Y, A_dict_sparse, grid_sparse, K);
            est = sort(est(:));
            if numel(est) == K && all(isfinite(est))
                err_mc(3) = mean((t_sorted - est).^2);
                ok_mc(3)  = true;
            end
        catch; end

        % 4. NC-ESPRIT
        try
            est = sort(NC_ESPRIT(Y, K));
            est = est(:);
            if numel(est) == K && all(isfinite(est))
                err_mc(4) = mean((t_sorted - est).^2);
                ok_mc(4)  = true;
            end
        catch; end

        % 5. NC-L1-SVD
        try
            [~, est] = NC_L1_SVD(Y, A_dict_sparse, grid_sparse, K, sigma2, []);
            est = sort(est(:));
            if numel(est) == K && all(isfinite(est))
                err_mc(5) = mean((t_sorted - est).^2);
                ok_mc(5)  = true;
            end
        catch; end

        err(mc, :) = err_mc;
        ok(mc, :)  = ok_mc;
    end

    for m = 1:n_methods
        valid = ok(:, m);
        RATE(m, ni) = mean(valid);
        if any(valid)
            RMSE(m, ni) = sqrt(mean(err(valid, m)));
        end
    end

    fprintf('N_snap=%4d  ', N_snap);
    for m = 1:n_methods
        fprintf('%s %.0f%%  ', methods{m}, RATE(m,ni)*100);
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

figure('Color','w','Position',[100 100 800 520]);

for m = 1:n_methods
    semilogy(Nsnap_vec, RMSE(m,:), ...
        [lines{m}, markers{m}], ...
        'Color', colors(m,:), 'LineWidth', 2, ...
        'MarkerFaceColor', colors(m,:), 'MarkerSize', 7, ...
        'DisplayName', methods{m});
    hold on;
end
semilogy(Nsnap_vec, crb_curve, 'k-', 'LineWidth', 2.5, 'DisplayName', 'NC-CRB');

set(gca, 'XScale', 'log');
grid on; box on;
xlim([Nsnap_vec(1), Nsnap_vec(end)]);
xlabel('快拍数 N', 'FontSize', 11);
ylabel('RMSE (°)',  'FontSize', 11);
title(sprintf('RMSE vs 快拍数（SNR=%ddB，含 NC-CRB 下界）', SNR_dB), 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 9, 'NumColumns', 2);
text(0.02, 0.98, ...
    sprintf('M=%d, K=%d, SNR=%ddB, MC=%d\n随机角度 ±60°, 最小间隔 10°\nSBL:1° / MUSIC:0.5° / OMP&L1:0.1°', ...
    M_arr, K, SNR_dB, n_MC), ...
    'Units','normalized','VerticalAlignment','top', ...
    'FontSize', 8.5, 'BackgroundColor', [0.97 0.97 0.97], 'EdgeColor', [0.6 0.6 0.6]);

fprintf('\n仿真完成（快拍数扫描）。\n');