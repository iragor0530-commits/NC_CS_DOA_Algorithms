%% NC-MUSIC 性能曲线仿真
clear; clc; close all;

%% ================== 公共参数 ==================
M          = 8;
K          = 3;
n_MC       = 500;
theta_true =sort(generate_angles(K,-60,60,8));
grid_theta = -90:0.5:90;       % NC_MUSIC 所需角度扫描网格


%% ================== (1) RMSE vs 快拍数 ==================
SNR2      = 5;
sigma2_2  = 10^(-SNR2/10);
Nsnap_vec = [1,2,5,8,10, 20, 50, 100, 200, 500, 1000];
rmse_snap = zeros(1, length(Nsnap_vec));

for ni = 1:length(Nsnap_vec)
    N    = Nsnap_vec(ni);
    err2 = zeros(1, n_MC);
    for mc = 1:n_MC
        A = exp(1j*pi*(0:M-1).' * sind(theta_true));
        S = randn(K, N);
        noise = sqrt(sigma2_2/2)*(randn(M,N)+1j*randn(M,N));
        Y = A*S + noise;
        try
            [~, est] = NC_MUSIC(Y, grid_theta, K);   % 使用 NC_MUSIC
            est = sort(est);
            err2(mc) = mean((sort(theta_true) - est).^2);
        catch
            err2(mc) = NaN;
        end
    end
    rmse_snap(ni) = sqrt(mean(err2));
end

%% ================== (2) 相干源 vs 非相干源 ==================
SNR_vec3   = -20:2:20;
N_snap3    = 200;
rho        = 1.0;
rmse_coh   = zeros(1, length(SNR_vec3));
rmse_incoh = zeros(1, length(SNR_vec3));

for si = 1:length(SNR_vec3)
    sigma2    = 10^(-SNR_vec3(si)/10);
    err_coh   = zeros(1, n_MC);
    err_incoh = zeros(1, n_MC);
    for mc = 1:n_MC
        A  = exp(1j*pi*(0:M-1).' * sind(theta_true));
        s1 = randn(1, N_snap3);
        S_coh   = [s1; rho*s1; randn(1,N_snap3)];
        S_incoh = randn(K, N_snap3);
        noise   = sqrt(sigma2/2)*(randn(M,N_snap3)+1j*randn(M,N_snap3));
        try
            [~, est] = NC_MUSIC(A*S_coh+noise, grid_theta, K);   % 使用 NC_MUSIC
            err_coh(mc) = mean((sort(theta_true) - sort(est)).^2);
        catch; err_coh(mc) = NaN; end
        try
            [~, est] = NC_MUSIC(A*S_incoh+noise, grid_theta, K); % 使用 NC_MUSIC
            err_incoh(mc) = mean((sort(theta_true) - sort(est)).^2);
        catch; err_incoh(mc) = NaN; end
    end
    rmse_coh(si)   = sqrt(mean(err_coh));
    rmse_incoh(si) = sqrt(mean(err_incoh));
end

%% ================== 绘图辅助函数 ==================
function label_points(xdata, ydata, color, ax)
    for i = 1:length(xdata)
        if ~isnan(ydata(i))
            text(ax, xdata(i), ydata(i)*1.15, ...
                sprintf('%.3f', ydata(i)), ...
                'FontSize', 6.5, 'Color', color, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'bottom');
        end
    end
end

%% ================== 绘图 ==================
rmse_incoh = max(rmse_incoh, 1e-3);   % clip 防止 log(0)
rmse_coh   = max(rmse_coh,   1e-3);
figure('Color','w','Position',[100 100 1350 430]);


% ---------- 子图1：RMSE vs 快拍数 ----------
ax2 = subplot(1,2,1);
c2 = [0.85 0.33 0.10];
semilogx(Nsnap_vec, rmse_snap, '-s', 'Color',c2, 'LineWidth',2, ...
    'MarkerFaceColor',c2, 'MarkerSize',6);
label_points(Nsnap_vec, rmse_snap, c2, ax2);
grid on; box on;
title('RMSE vs 快拍数', 'FontSize',11);
xlabel('快拍数 N', 'FontSize',10);
ylabel('均方根误差 RMSE (°)', 'FontSize',10);
text(0.04, 0.18, ...
    sprintf('NC-MUSIC\n信噪比 SNR = %d dB\n阵元数 M = %d\n信源数 K = %d\n非相干信号', SNR2, M, K), ...
    'Units','normalized', 'FontSize',8.5, ...
    'BackgroundColor',[1.0 0.97 0.93], 'EdgeColor',[0.8 0.5 0.2], ...
    'LineWidth',1.0, 'VerticalAlignment','bottom');

% ---------- 子图2：相干源 vs 非相干源 ----------
ax3 = subplot(1,2,2);
c3b = [0.15 0.45 0.78];
c3r = [0.80 0.10 0.10];
semilogy(SNR_vec3, rmse_incoh, '-o', 'Color',c3b, 'LineWidth',2, ...
    'MarkerFaceColor',c3b, 'MarkerSize',6); hold on;
semilogy(SNR_vec3, rmse_coh,   '--^', 'Color',c3r, 'LineWidth',2, ...
    'MarkerFaceColor',c3r, 'MarkerSize',6);
label_points(SNR_vec3, rmse_incoh, c3b, ax3);
label_points(SNR_vec3, rmse_coh,   c3r, ax3);
grid on; box on;
title('相干源 vs 非相干源', 'FontSize',11);
xlabel('信噪比 (dB)', 'FontSize',10);
ylabel('均方根误差 RMSE (°)', 'FontSize',10);
ylim([1e-2, 1e2]);
legend('非相干信号', sprintf('相干信号 (ρ=%.1f)', rho), ...
    'Location','northeast', 'FontSize',9);
text(0.04, 0.18, ...
    sprintf('NC-MUSIC\n快拍数 N = %d\n阵元数 M = %d\n信源数 K = %d', N_snap3, M, K), ...
    'Units','normalized', 'FontSize',8.5, ...
    'BackgroundColor',[1.0 0.97 0.97], 'EdgeColor',[0.8 0.3 0.3], ...
    'LineWidth',1.0, 'VerticalAlignment','bottom');

fprintf('\n仿真完成。\n');