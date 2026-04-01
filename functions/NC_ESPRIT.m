function est_theta = NC_ESPRIT(Y, K_source)
    [M, N] = size(Y);

    %% Step 1: 增广协方差矩阵
    % 流形约定 exp(+jπ sinθ)，所以 J*conj(Y) 的下半阵
    % 等效流形为 exp(-jπ sinθ)，联合后覆盖正负频
    J = fliplr(eye(M));
    Y_aug = [Y; J * conj(Y)];          % 2M × N
    R = (Y_aug * Y_aug') / N;

    %% Step 2: 信号子空间，取 K 个（非 2K）
    [U, ~] = eig(R);
    % eig 返回升序，取最大 K 个
    Us = U(:, end-K_source+1:end);     % 2M × K

    %% Step 3: 子阵构造
    E1 = Us([1:M-1,   M+1:2*M-1], :); % 2(M-1) × K
    E2 = Us([2:M,     M+2:2*M  ], :); % 2(M-1) × K

    %% Step 4: 旋转算子（最小二乘）
    Phi = E1 \ E2;                     % K × K

    %% Step 5: 特征值（K个，全部使用）
    lambda = eig(Phi);
    % 数值保护：取模最近单位圆的 K 个
    [~, si] = sort(abs(abs(lambda) - 1));
    lambda = lambda(si(1:K_source));

    %% Step 6: 角度提取
    % 仿真脚本流形为 exp(+jπ sinθ)
    % 旋转算子特征值 λ = exp(+jπ sinθ)
    % 因此 sinθ = angle(λ) / π
    sin_val   = real(angle(lambda) / pi);
    sin_val   = max(-1, min(1, sin_val));
    est_theta = sort(rad2deg(asin(sin_val)));
end