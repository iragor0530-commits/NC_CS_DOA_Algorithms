%% DOA_benchmark_coherent.m
%  五种 NC 算法在相干源条件下的 RMSE vs SNR
%  场景：K=3 信源，其中信源1&2完全相干，信源3独立
%  M=8, N=200, MC=300

clear; clc; close all;

if isempty(gcp('nocreate'))
    parpool('local', 4);
end

%% ========== 参数 ==========
M_arr   = 8;
K       = 3;
N_snap  = 200;
n_MC    = 300;
SNR_vec = -10:2.5:15;
n_snr   = length(SNR_vec);

%% ========== 网格 / 字典 ==========
grid_SBL    = -90:1:90;
grid_MUSIC  = -90:0.5:90;
grid_sparse = -90:0.1:90;
A_dict_sparse = exp(1j*pi*(0:M_arr-1).' * sind(grid_sparse));

%% ========== 预生成角度池 ==========
rng(42);
theta_pool = zeros(n_MC, K);
for mc = 1:n_MC
    theta_pool(mc,:) = generate_angles(K, -60, 60, 10);
end
theta_rep = mean(theta_pool, 1);

%% ========== NC-CRB（非相干，作参考下界） ==========
crb_curve = zeros(1, n_snr);
for si = 1:n_snr
    cv = compute_NC_CRB_internal(theta_rep, M_arr, SNR_vec(si), N_snap);
    crb_curve(si) = sqrt(mean(cv.^2));
end

%% ========== 仿真 ==========
methods   = {'NC-SBL','NC-MUSIC','NC-OMP','NC-ESPRIT','NC-L1-SVD'};
n_methods = 5;
RMSE = nan(n_methods, n_snr);
RATE = zeros(n_methods, n_snr);

for si = 1:n_snr
    snr_db = SNR_vec(si);
    sigma2 = 10^(-snr_db/10);

    err = nan(n_MC, n_methods);
    ok  = false(n_MC, n_methods);

    parfor mc = 1:n_MC
        theta_true = theta_pool(mc,:);
        A = exp(1j*pi*(0:M_arr-1).' * sind(theta_true));

        % 相干信源矩阵：信源1&2共享同一基带，信源3独立
        s_base = sign(randn(1, N_snap));
        S = [s_base; s_base; sign(randn(1, N_snap))];

        noise = sqrt(sigma2/2)*(randn(M_arr,N_snap)+1j*randn(M_arr,N_snap));
        Y = A*S + noise;

        err_mc = nan(1, n_methods);
        ok_mc  = false(1, n_methods);
        t_sorted = sort(theta_true(:));

        try
            est = sort(NC_SBL(Y, grid_SBL)); est = est(1:min(end,K));
            if numel(est)==K && all(isfinite(est))
                err_mc(1)=mean((t_sorted-est(:)).^2); ok_mc(1)=true; end
        catch; end
        try
            [~,est]=NC_MUSIC(Y,grid_MUSIC,K); est=sort(est(:));
            if numel(est)==K && all(isfinite(est))
                err_mc(2)=mean((t_sorted-est).^2); ok_mc(2)=true; end
        catch; end
        try
            [~,est]=NC_OMP(Y,A_dict_sparse,grid_sparse,K); est=sort(est(:));
            if numel(est)==K && all(isfinite(est))
                err_mc(3)=mean((t_sorted-est).^2); ok_mc(3)=true; end
        catch; end
        try
            est=sort(NC_ESPRIT(Y,K)); est=est(:);
            if numel(est)==K && all(isfinite(est))
                err_mc(4)=mean((t_sorted-est).^2); ok_mc(4)=true; end
        catch; end
        try
            [~,est]=NC_L1_SVD(Y,A_dict_sparse,grid_sparse,K,sigma2,[]); est=sort(est(:));
            if numel(est)==K && all(isfinite(est))
                err_mc(5)=mean((t_sorted-est).^2); ok_mc(5)=true; end
        catch; end

        err(mc,:)=err_mc; ok(mc,:)=ok_mc;
    end

    for m = 1:n_methods
        v = ok(:,m); RATE(m,si)=mean(v);
        if any(v), RMSE(m,si)=sqrt(mean(err(v,m))); end
    end

    fprintf('SNR=%+5.1fdB  ', snr_db);
    for m=1:n_methods, fprintf('%s %.0f%%  ',methods{m},RATE(m,si)*100); end
    fprintf('\n');
end

%% ========== 绘图 ==========
colors  = [0.15 0.45 0.78;
           0.80 0.10 0.10;
           0.13 0.70 0.35;
           0.85 0.55 0.00;
           0.55 0.20 0.75];
markers = {'o','^','s','d','v'};
lines   = {'-','--','-.',':','-'};

figure('Color','w','Position',[100 100 780 520]);
for m = 1:n_methods
    semilogy(SNR_vec, RMSE(m,:), [lines{m},markers{m}], ...
        'Color',colors(m,:),'LineWidth',2,...
        'MarkerFaceColor',colors(m,:),'MarkerSize',7,...
        'DisplayName',methods{m});
    hold on;
end
semilogy(SNR_vec, crb_curve, 'k-', 'LineWidth',2.5, 'DisplayName','NC-CRB（非相干参考）');

grid on; box on;
ylim([1e-2, 1e2]); xlim([SNR_vec(1), SNR_vec(end)]);
xlabel('信噪比 (dB)', 'FontSize',12);
ylabel('RMSE (°)',    'FontSize',12);
title('相干源场景下 RMSE vs SNR', 'FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',9,'NumColumns',2);
text(0.02,0.98, ...
    sprintf('M=%d, K=%d（2相干+1独立）, N=%d, MC=%d\n随机角度 ±60°, 最小间隔 10°', ...
    M_arr,K,N_snap,n_MC), ...
    'Units','normalized','VerticalAlignment','top',...
    'FontSize',9,'BackgroundColor',[0.97 0.97 0.97],'EdgeColor',[0.6 0.6 0.6]);

fprintf('\n仿真完成。\n');