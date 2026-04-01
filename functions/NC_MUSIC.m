function [P_db, est_theta] = NC_MUSIC(Y, grid_theta, M_source)
    % NC_MUSIC 基于非圆扩展模型的 MUSIC DOA 估计
    %
    % 输入:
    %   Y          - 接收数据矩阵 [N x K]
    %   grid_theta - 角度网格 [1 x P]
    %   M_source   - 信源数量
    %
    % 输出:
    %   P_db       - 归一化空间谱 (dB)
    %   est_theta  - 估计角度 (升序)

    [N, K] = size(Y);

    %% 1. 非圆扩展协方差矩阵
    Y_aug = [Y; conj(Y)];               % [2N x K]
    R_aug = (Y_aug * Y_aug') / K;       % [2N x 2N]

    %% 2. 特征分解，提取噪声子空间
    [E, D]   = eig(R_aug);
    [~, idx] = sort(real(diag(D)), 'ascend');
    E_n      = E(:, idx(1:2*N-M_source));  % 噪声子空间 [2N x (2N-M)]
    En_En    = E_n * E_n';

    %% 3. 扫描空间谱
    x_n     = (0:N-1)';
    P_music = zeros(1, length(grid_theta));

    for i = 1:length(grid_theta)
        a          = exp(1j * pi * x_n * sin(deg2rad(grid_theta(i))));
        a_aug      = [a; conj(a)];
        P_music(i) = 1 / real(a_aug' * En_En * a_aug);
    end

    %% 4. 归一化转 dB
    P_db = 10 * log10(P_music / max(P_music));
    P_db = max(P_db(:), -60);

    %% 5. 峰值提取
    [~, locs] = findpeaks(P_db, 'MinPeakProminence', 6);

    if isempty(locs)
        [~, idx2] = sort(P_db, 'descend');
        est_theta  = grid_theta(idx2(1:M_source));
    else
        [~, sort_idx] = sort(P_db(locs), 'descend');
        actual_k      = min(length(locs), M_source);
        est_theta     = grid_theta(locs(sort_idx(1:actual_k)));
        if actual_k < M_source
            [~, idx2] = sort(P_db, 'descend');
            for ii = 1:length(idx2)
                if length(est_theta) >= M_source, break; end
                if ~ismember(grid_theta(idx2(ii)), est_theta)
                    est_theta(end+1) = grid_theta(idx2(ii)); %#ok
                end
            end
        end
    end

    est_theta = sort(est_theta);
end