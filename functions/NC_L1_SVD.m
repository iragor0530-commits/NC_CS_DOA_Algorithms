function [P_db, est_theta] = NC_L1_SVD(Y, A_dict, grid_theta, M_source, sigma2, lambda_in)
    % NC_L1_SVD 基于非圆扩展模型的 L1-SVD DOA 估计 (ADMM求解器)
    %
    % 输入:
    %   Y          - 接收数据矩阵 [N x K]
    %   A_dict     - 字典矩阵 [N x P]
    %   grid_theta - 角度网格 [1 x P]
    %   M_source   - 信源数量
    %   sigma2     - 噪声功率（可选，已知时epsilon更准）
    %   lambda_in  - ADMM惩罚参数（可选，由CVX校准传入）

    %% 1. 非圆扩展
    Y_aug = [Y; conj(Y)];           % [2N x K]
    Psi   = [A_dict; conj(A_dict)]; % [2N x P]

    [M2, K] = size(Y_aug);
    L = length(grid_theta);

    %% 2. SVD 降维
    [U, S_val, ~] = svd(Y_aug, 'econ');
    s_vals = diag(S_val);
    Y_sv   = U(:, 1:M_source) * S_val(1:M_source, 1:M_source); % [2N x M_source]

    %% 3. epsilon 计算（用于 lambda 推导）
    if nargin >= 5 && ~isempty(sigma2)
        % 已知噪声功率：理论公式，最准确
        epsilon = 0.8 * sqrt(sigma2 * M_source * K);
    else
        % 盲估：从噪声子空间奇异值反推
        if length(s_vals) > M_source
            noise_sv   = s_vals(M_source+1:end);
            sigma2_est = sum(noise_sv.^2) / (length(noise_sv) * K);
        else
            sigma2_est = norm(Y_sv,'fro')^2 / (M2 * K) * 0.1;
        end
        epsilon = 0.8 * sqrt(sigma2_est * M_source * K);
    end

    %% 4. ADMM 稀疏重构
    % 求解: min sum_i‖S(i,:)‖₂  s.t. ‖Y_sv - Psi·S‖_F ≤ epsilon
    % 等价惩罚形式: min (1/2)‖Y_sv - Psi·S‖_F² + lambda·sum_i‖S(i,:)‖₂

    % lambda：优先用CVX校准值，否则从epsilon推导
    if nargin >= 6 && ~isempty(lambda_in)
        lambda = lambda_in;
    else
        lambda = epsilon / (sqrt(M_source * K) + 1e-6) * 0.3;
    end

    % ADMM 参数
    rho    = 1.0;
    max_it = 200;
    tol    = 1e-5;

    % 预计算（循环外做一次，速度关键）
    PsiTPsi = Psi' * Psi;                           % [P x P]
    PsiTY   = Psi' * Y_sv;                          % [P x M_source]
    A_mat   = PsiTPsi + rho * eye(L);
    A_chol  = chol(A_mat, 'upper');                 % Cholesky 预分解

    % 初始化
    S_hat = zeros(L, M_source);
    Z     = zeros(L, M_source);
    Uu    = zeros(L, M_source);

    for it = 1:max_it
        % S 更新：线性系统（用预分解加速）
        rhs   = PsiTY + rho * (Z - Uu);
        S_hat = A_chol \ (A_chol' \ rhs);

        % Z 更新：逐行 Group Lasso 近端算子（块软阈值）
        V_z       = S_hat + Uu;
        row_norms = sqrt(sum(abs(V_z).^2, 2)) + 1e-12;
        scale     = max(1 - (lambda/rho) ./ row_norms, 0);
        Z_new     = bsxfun(@times, scale, V_z);

        % U 更新：对偶变量
        Uu = Uu + S_hat - Z_new;

        % 收敛判断
        pri_res = norm(S_hat - Z_new, 'fro') / (norm(Z_new,'fro') + 1e-12);
        if pri_res < tol
            break;
        end
        Z = Z_new;
    end
    S_hat = Z_new;  % 取 Z 作为最终解（更稀疏）

    %% 5. 空间谱计算
    P_spatial = sum(abs(S_hat).^2, 2);

    if max(P_spatial) == 0
        P_db      = -60 * ones(L, 1);
        est_theta = grid_theta(1:M_source);
        return;
    end

    P_db = 10 * log10(P_spatial / max(P_spatial));
    P_db = max(P_db, -60);

    %% 6. 峰值提取
    [~, locs] = findpeaks(P_db, 'MinPeakProminence', 6);

    if isempty(locs)
        [~, idx] = sort(P_db, 'descend');
        est_theta = grid_theta(idx(1:M_source));
    else
        [~, sort_idx] = sort(P_db(locs), 'descend');
        actual_k  = min(length(locs), M_source);
        est_theta = grid_theta(locs(sort_idx(1:actual_k)));
        if actual_k < M_source
            [~, idx] = sort(P_db, 'descend');
            for ii = 1:length(idx)
                if length(est_theta) >= M_source, break; end
                if ~ismember(grid_theta(idx(ii)), est_theta)
                    est_theta(end+1) = grid_theta(idx(ii)); %#ok
                end
            end
        end
    end

    est_theta = sort(est_theta);
end