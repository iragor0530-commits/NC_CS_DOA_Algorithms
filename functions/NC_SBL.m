function est_theta = NC_SBL(Y, grid_theta)
    [M, N] = size(Y);
    G  = length(grid_theta);
    M2 = 2*M;

    Phi     = exp(1j*pi*(0:M-1)' * sin(deg2rad(grid_theta)));
    Phi_aug = [Phi; conj(Phi)];
    Y_aug   = [Y;   conj(Y)];

    %% SBL 迭代
    gamma  = ones(G,1);
    sigma2 = 0.1*var(Y_aug(:));

    for iter = 1:200
        gamma_old = gamma;
        Sigma_y = Phi_aug*diag(gamma)*Phi_aug' + sigma2*eye(M2);
        K_post  = diag(gamma)*Phi_aug'/Sigma_y;
        Mu      = K_post*Y_aug;
        Sigma   = diag(gamma) - K_post*Phi_aug*diag(gamma);
        gamma   = max(real(sum(abs(Mu).^2,2)/N + diag(Sigma)), 1e-10);
        residual = Y_aug - Phi_aug*Mu;
        dof    = M2 - sum(1 - diag(Sigma)./gamma);
        sigma2 = max(real(norm(residual,'fro')^2/(N*max(dof,1))), 1e-10);
        if norm(gamma-gamma_old)/norm(gamma_old) < 1e-4
            break;
        end
    end

    R_aug = Sigma_y \ eye(M2); 
    %% 峰值检测：加 prominence 阈值过滤假峰
    grid_step    = mean(diff(grid_theta));
    min_peak_idx = max(1, round(3/grid_step));
    gamma_norm   = gamma / max(gamma);
    
    % 两端补 0，解决 findpeaks 边界盲区
    gamma_pad = [0; gamma_norm; 0];
    
    [peak_vals, locs] = findpeaks(gamma_pad, ...
        'MinPeakDistance',   min_peak_idx, ...
        'MinPeakProminence', 0.05, ...
        'MinPeakHeight',     0.03);
    
    % 补偿 padding 引入的索引偏移，并做边界保护
    locs = locs - 1;
    locs = max(1, min(G, locs));
    peak_vals = gamma_norm(locs);  % 用原始 gamma_norm 重新取值
    
    if isempty(locs)
        [~, locs] = max(gamma_norm);
        peak_vals = gamma_norm(locs);
    end
    
    % 按 gamma 从大到小排序
    [~, si] = sort(peak_vals, 'descend');
    locs     = locs(si);    

    %% 抛物线插值细化
est_theta  = zeros(1, length(locs));
coarse_est = grid_theta(locs);

for p = 1:length(locs)
    th = coarse_est(p);
    
    for iter = 1:20
        th_rad  = th * pi/180;
        m_vec   = (0:M-1)';
        a       = exp(1j*pi*m_vec*sin(th_rad));
        a_aug   = [a; conj(a)];
        da      = 1j*pi*m_vec*cos(th_rad)*(pi/180) .* a;
        da_aug  = [da; conj(da)];
        d2a     = -(pi*cos(th_rad)*(pi/180))^2 * m_vec.^2 .* a;
        d2a_aug = [d2a; conj(d2a)];

        g = 2*real(da_aug' * R_aug * a_aug);
        h = 2*real(da_aug' * R_aug * da_aug + d2a_aug' * R_aug * a_aug);

        if abs(h) < 1e-10, break; end
        step = g / h;
        
        % 信任域：步长限制在粗网格步长的一半
        step = sign(step) * min(abs(step), 0.5);
        th   = th - step;
        th   = max(-90, min(90, th));
        if abs(step) < 1e-6, break; end
    end
    
    % 只有在合理范围内才接受细化结果
    if abs(th - coarse_est(p)) <= 1.5  % 放宽到1.5°
        est_theta(p) = th;
    else
        est_theta(p) = coarse_est(p);
    end
end
end  % 函数结束