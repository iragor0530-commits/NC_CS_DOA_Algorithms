%% NC-OMP 性能曲线仿真
clear; clc; close all;

%% ================== 公共参数 ==================
M_arr      = 8;                          % 物理阵元数
K          = 3;                          % 信源数
n_MC       = 100;                        % 蒙特卡洛次数
theta_true = sort(generate_angles(K,-70,70,8));  % 随机生成真实角度
grid_theta = -90:0.1:90;                 % 字典网格

%% ================== (1) RMSE vs 快拍数 ==================
SNR2      = 10;
sigma2_2  = 10^(-SNR2/10);
Nsnap_vec = [1,2,5,8,10,20,50,100,200,500,1000];
rmse_snap = zeros(1, length(Nsnap_vec));

for ni = 1:length(Nsnap_vec)
    N     = Nsnap_vec(ni);
    err2  = zeros(1, n_MC);
    valid = false(1, n_MC);

    parfor mc = 1:n_MC
        % 构建导向矩阵与字典（相同网格）
        A_dict = exp(1j*pi*(0:M_arr-1).' * sind(grid_theta));   % (M x P)
        A_true = exp(1j*pi*(0:M_arr-1).' * sind(theta_true));   % (M x K)

        S     = sign(randn(K, N));
        noise = sqrt(sigma2_2/2) * (randn(M_arr,N) + 1j*randn(M_arr,N));
        Y     = A_true * S + noise;

        try
            [~, est] = NC_OMP(Y, A_dict, grid_theta, K);        % ← 正确签名
            if length(est) == K
                err2(mc)  = mean((sort(theta_true) - sort(est)).^2);
                valid(mc) = true;
            end
        catch
            err2(mc) = NaN;
        end
    end

    if any(valid)
        rmse_snap(ni) = sqrt(sum(err2(valid)) / n_MC);
    else
        rmse_snap(ni) = NaN;
    end
end

%% ================== (2) 相干源 vs 非相干源 ==================
SNR_vec3   = -10:2:20;
N_snap3    = 200;
rho        = 0.9;
rmse_coh   = zeros(1, length(SNR_vec3));
rmse_incoh = zeros(1, length(SNR_vec3));

for si = 1:length(SNR_vec3)
    sigma2    = 10^(-SNR_vec3(si)/10);
    err_coh   = zeros(1, n_MC);
    err_incoh = zeros(1, n_MC);
    valid_coh = false(1, n_MC);
    valid_inc = false(1, n_MC);

    parfor mc = 1:n_MC
        A_dict = exp(1j*pi*(0:M_arr-1).' * sind(grid_theta));   % (M x P)
        A_true = exp(1j*pi*(0:M_arr-1).' * sind(theta_true));   % (M x K)

        s1      = sign(randn(1, N_snap3));
        S_coh   = [s1; rho*s1; sign(randn(1, N_snap3))];        % 相干信源
        S_incoh = sign(randn(K, N_snap3));                       % 非相干信源
        noise   = sqrt(sigma2/2) * (randn(M_arr,N_snap3) + 1j*randn(M_arr,N_snap3));

        % 相干源
        try
            [~, est] = NC_OMP(A_true*S_coh + noise, A_dict, grid_theta, K);
            if length(est) == K
                err_coh(mc)   = mean((sort(theta_true) - sort(est)).^2);
                valid_coh(mc) = true;
            end
        catch
            err_coh(mc) = NaN;
        end

        % 非相干源
        try
            [~, est] = NC_OMP(A_true*S_incoh + noise, A_dict, grid_theta, K);
            if length(est) == K
                err_incoh(mc)  = mean((sort(theta_true) - sort(est)).^2);
                valid_inc(mc)  = true;
            end
        catch
            err_incoh(mc) = NaN;
        end
    end

    rmse_coh(si)   = sqrt(sum(err_coh(valid_coh))   / n_MC);
    rmse_incoh(si) = sqrt(sum(err_incoh(valid_inc)) / n_MC);
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

% ---------- 图1：RMSE vs 快拍数 ----------
figure('Color','w','Name','RMSE vs 快拍数');
ax1 = axes;
c1  = [0.85 0.33 0.10];

semilogx(Nsnap_vec, rmse_snap, '-s', 'Color',c1, 'LineWidth',2, ...
    'MarkerFaceColor',c1, 'MarkerSize',6);
hold on;
label_points(Nsnap_vec, rmse_snap, c1, ax1);

grid on; box on;
title('RMSE vs 快拍数', 'FontSize',11);
xlabel('快拍数', 'FontSize',10);
ylabel('RMSE (°)', 'FontSize',10);

text(0.04, 0.18, ...
    sprintf('NC-OMP\nSNR = %d dB\nM = %d\nK = %d\nIncoherent Signals', SNR2, M_arr, K), ...
    'Units','normalized', 'FontSize',8.5, ...
    'BackgroundColor',[1.0 0.97 0.93], 'EdgeColor',[0.8 0.5 0.2], ...
    'LineWidth',1.0, 'VerticalAlignment','bottom');

% ---------- 图2：相干源 vs 非相干源 ----------
figure('Color','w','Name','相干源信号 vs 非相干源信号');
ax2 = axes;
c2b = [0.15 0.45 0.78];
c2r = [0.80 0.10 0.10];

semilogy(SNR_vec3, rmse_incoh, '-o', 'Color',c2b, 'LineWidth',2, ...
    'MarkerFaceColor',c2b, 'MarkerSize',6);
hold on;
semilogy(SNR_vec3, rmse_coh, '--^', 'Color',c2r, 'LineWidth',2, ...
    'MarkerFaceColor',c2r, 'MarkerSize',6);

label_points(SNR_vec3, rmse_incoh, c2b, ax2);
label_points(SNR_vec3, rmse_coh,   c2r, ax2);

grid on; box on;
title('相干源信号 vs 非相干源信号', 'FontSize',11);
xlabel('SNR (dB)', 'FontSize',10);
ylabel('RMSE (°)', 'FontSize',10);
ylim([1e-2, 1e2]);

legend('非相干源', sprintf('相干源(\\rho = %.1f)', rho), ...
    'Location','northeast', 'FontSize',9);

text(0.04, 0.18, ...
    sprintf('NC-OMP\nN = %d\nM = %d\nK = %d', N_snap3, M_arr, K), ...
    'Units','normalized', 'FontSize',8.5, ...
    'BackgroundColor',[1.0 0.97 0.97], 'EdgeColor',[0.8 0.3 0.3], ...
    'LineWidth',1.0, 'VerticalAlignment','bottom');

fprintf('\n仿真完成。\n');