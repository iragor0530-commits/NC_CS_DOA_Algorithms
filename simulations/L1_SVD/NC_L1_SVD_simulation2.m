%% NC-L1-SVD 性能曲线仿真：RMSE vs 快拍数 & 相干/非相干源对比
clc; clear; close all;

%% ===== 公共参数 =====
N           = 8;
M_source    = 3;
n_MC        = 200;
FAIL_THRESH = 15;
rho         = 0.9;

grid_step  = 0.1;
grid_theta = -90:grid_step:90;
x_n        = (0:N-1)';
A_dict     = exp(1j * pi * x_n * sin(deg2rad(grid_theta)));

theta_fixed = sort(generate_angles(M_source, -50, 50, 10));
fprintf('固定真实角度: %s 度\n', num2str(theta_fixed));

if isempty(gcp('nocreate')); parpool; end

%% ===== CVX 单次校准 lambda（公共，只跑一次）=====
fprintf('\n正在用CVX校准lambda...\n');
K_cal  = 200;
SNR_cal = 10;
S_cal  = (randn(M_source,K_cal)>0)*2-1;
A_cal  = exp(1j*pi*x_n*sin(deg2rad(theta_fixed)));
As_cal = A_cal*S_cal;
sig2_cal = mean(abs(As_cal(:)).^2)/(10^(SNR_cal/10));
V_cal  = sqrt(sig2_cal/2)*(randn(N,K_cal)+1j*randn(N,K_cal));
Y_cal  = As_cal + V_cal;

Y_aug_c = [Y_cal; conj(Y_cal)];
Psi_c   = [A_dict; conj(A_dict)];
L       = length(grid_theta);
[U_c, S_c, ~] = svd(Y_aug_c, 'econ');
s_c     = diag(S_c);
Y_sv_c  = U_c(:,1:M_source)*S_c(1:M_source,1:M_source);
noise_sv_c = s_c(M_source+1:end);
sigma2_c   = sum(noise_sv_c.^2)/(length(noise_sv_c)*K_cal);
epsilon_c  = 0.8*sqrt(sigma2_c*M_source*K_cal);

cvx_begin quiet
    variable S_cvx(L, M_source) complex;
    variable r_cvx(L) nonnegative;
    minimize(sum(r_cvx))
    subject to
        norm(Y_sv_c - Psi_c*S_cvx,'fro') <= epsilon_c;
        for i = 1:L
            norm(S_cvx(i,:),2) <= r_cvx(i);
        end
cvx_end

S_cvx_val  = double(S_cvx);
resid_cvx  = norm(Y_sv_c - Psi_c*S_cvx_val,'fro');
sparse_sum = sum(sqrt(sum(abs(S_cvx_val).^2,2)));
if sparse_sum > 1e-8
    lambda_ref = resid_cvx/(sparse_sum+1e-8)*10;
else
    lambda_ref = epsilon_c/(sqrt(M_source*K_cal)+1e-6)*0.3;
end
fprintf('CVX校准完成 → lambda_ref=%.6f\n', lambda_ref);

%% ================== (1) RMSE vs 快拍数 ==================
SNR_snap  = 10;
Nsnap_vec = [3,5,10,20,50,100,200,500,1000];  % 最小快拍数 >= M_source，避免SVD秩不足
rmse_snap = nan(1, length(Nsnap_vec));

fprintf('\n实验1：RMSE vs 快拍数  (SNR=%ddB)\n', SNR_snap);
for ni = 1:length(Nsnap_vec)
    K_s    = Nsnap_vec(ni);
    all_sq = nan(n_MC, M_source);
    lam    = lambda_ref;

    parfor mc = 1:n_MC
        S_mc      = (randn(M_source,K_s)>0)*2-1;
        A_mc      = exp(1j*pi*(0:N-1)'*sin(deg2rad(theta_fixed)));
        As_mc     = A_mc*S_mc;
        sigma2_mc = mean(abs(As_mc(:)).^2)/(10^(SNR_snap/10));
        V_mc      = sqrt(sigma2_mc/2)*(randn(N,K_s)+1j*randn(N,K_s));
        Y_mc      = As_mc + V_mc;

        try
            [~, est] = NC_L1_SVD(Y_mc, A_dict, grid_theta, M_source, sigma2_mc, lam);
            est_s  = sort(est(:));
            true_s = sort(theta_fixed(:));
            if length(est_s)==M_source && max(abs(est_s-true_s))<FAIL_THRESH
                all_sq(mc,:) = (est_s - true_s).^2;
            end
        catch; end
    end

    valid = all_sq(~isnan(all_sq(:)));
    if ~isempty(valid)
        rmse_snap(ni) = sqrt(mean(valid));
    end
    fprintf('  N=%4d  RMSE=%.4f°\n', K_s, rmse_snap(ni));
end

%% ================== (2) 相干源 vs 非相干源 ==================
SNR_vec    = -10:2:20;
K_coh      = 200;
rmse_incoh = nan(1, length(SNR_vec));
rmse_coh   = nan(1, length(SNR_vec));

fprintf('\n实验2：相干源 vs 非相干源  (K=%d)\n', K_coh);
for si = 1:length(SNR_vec)
    cur_SNR = SNR_vec(si);
    sigma2  = 10^(-cur_SNR/10);
    sq_incoh = nan(n_MC, M_source);
    sq_coh   = nan(n_MC, M_source);
    lam      = lambda_ref;

    parfor mc = 1:n_MC
        A_mc  = exp(1j*pi*(0:N-1)'*sin(deg2rad(theta_fixed)));
        noise = sqrt(sigma2/2)*(randn(N,K_coh)+1j*randn(N,K_coh));

        % 非相干
        S_incoh  = (randn(M_source,K_coh)>0)*2-1;
        As_incoh = A_mc*S_incoh;
        sig2_i   = mean(abs(As_incoh(:)).^2)/(10^(cur_SNR/10));
        Y_incoh  = As_incoh + sqrt(sig2_i/2)*(randn(N,K_coh)+1j*randn(N,K_coh));

        % 相干
        s1     = (randn(1,K_coh)>0)*2-1;
        S_coh  = [s1; rho*s1+sqrt(1-rho^2)*((randn(1,K_coh)>0)*2-1); (randn(1,K_coh)>0)*2-1];
        As_coh = A_mc*S_coh;
        sig2_c = mean(abs(As_coh(:)).^2)/(10^(cur_SNR/10));
        Y_coh  = As_coh + sqrt(sig2_c/2)*(randn(N,K_coh)+1j*randn(N,K_coh));

        % 非相干估计
        [~, est] = NC_L1_SVD(Y_incoh, A_dict, grid_theta, M_source, sig2_i, lam);
        est_s = sort(est(:)); true_s = sort(theta_fixed(:));
        if length(est_s)==M_source && max(abs(est_s-true_s))<FAIL_THRESH
            sq_incoh(mc,:) = (est_s-true_s).^2;
        end

        % 相干估计
        [~, est] = NC_L1_SVD(Y_coh, A_dict, grid_theta, M_source, sig2_c, lam);
        est_s = sort(est(:));
        if length(est_s)==M_source && max(abs(est_s-true_s))<FAIL_THRESH
            sq_coh(mc,:) = (est_s-true_s).^2;
        end
    end

    v_i = sq_incoh(~isnan(sq_incoh(:)));
    v_c = sq_coh(~isnan(sq_coh(:)));
    if ~isempty(v_i), rmse_incoh(si) = sqrt(mean(v_i)); end
    if ~isempty(v_c), rmse_coh(si)   = sqrt(mean(v_c)); end
    fprintf('  SNR=%+3ddB  非相干=%.4f°  相干=%.4f°\n', cur_SNR, rmse_incoh(si), rmse_coh(si));
end

%% ================== 绘图辅助 ==================
function label_points(xdata, ydata, color, ax)
    for i = 1:length(xdata)
        if ~isnan(ydata(i))
            text(ax, xdata(i), ydata(i)*1.12, sprintf('%.3f', ydata(i)), ...
                'FontSize',6.5, 'Color',color, ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom');
        end
    end
end

%% ================== 绘图 ==================
c_blue = [0.15 0.45 0.78];
c_red  = [0.80 0.10 0.10];

%% 图1：RMSE vs 快拍数
figure('Color','w','Name','RMSE vs 快拍数');
ax1 = axes;
semilogx(Nsnap_vec, rmse_snap, '-s', 'Color',c_blue, 'LineWidth',2, ...
    'MarkerFaceColor',c_blue, 'MarkerSize',6);
hold on;
label_points(Nsnap_vec, rmse_snap, c_blue, ax1);
grid on; box on;
xlim([Nsnap_vec(1), Nsnap_vec(end)]);
xlabel('快拍数', 'FontSize',11);
ylabel('RMSE (°)', 'FontSize',11);
title('RMSE vs 快拍数', 'FontSize',12);
text(0.04, 0.18, ...
    sprintf('NC-L1-SVD\nSNR=%ddB\nN=%d\nM=%d\nIncoherent Signals', SNR_snap, N, M_source), ...
    'Units','normalized','FontSize',9, ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[1.0 0.97 0.93], 'EdgeColor',[0.8 0.5 0.2], 'LineWidth',1.0);

%% 图2：相干源 vs 非相干源
figure('Color','w','Name','相干源 vs 非相干源');
ax2 = axes;
semilogy(SNR_vec, rmse_incoh, '-o', 'Color',c_blue, 'LineWidth',2, ...
    'MarkerFaceColor',c_blue, 'MarkerSize',6);
hold on;
semilogy(SNR_vec, rmse_coh, '--^', 'Color',c_red, 'LineWidth',2, ...
    'MarkerFaceColor',c_red, 'MarkerSize',6);
label_points(SNR_vec, rmse_incoh, c_blue, ax2);
label_points(SNR_vec, rmse_coh,   c_red,  ax2);
grid on; box on;
ylim([1e-2, 1e2]); xlim([SNR_vec(1), SNR_vec(end)]);
xlabel('SNR (dB)', 'FontSize',11);
ylabel('RMSE (°)', 'FontSize',11);
title('相干源信号 vs 非相干源信号', 'FontSize',12);
legend('非相干源', sprintf('相干源(\\rho=%.1f)', rho), ...
    'Location','northeast', 'FontSize',10);
text(0.04, 0.18, ...
    sprintf('NC-L1-SVD\nK=%d\nN=%d\nM=%d', K_coh, N, M_source), ...
    'Units','normalized','FontSize',9, ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[1.0 0.97 0.97], 'EdgeColor',[0.8 0.3 0.3], 'LineWidth',1.0);

fprintf('\n仿真完成。\n');