function [P_db, est_angles] = NC_OMP(Y, A_dict, grid_theta, M)
% NC_OMP: 基于非圆扩展的正交匹配追踪算法
% 输入:
%   Y         : 原始接收信号矩阵 (N x K)
%   A_dict    : 原始字典矩阵 (N x P)，未归一化
%   grid_theta: 角度搜索网格 (1 x P)
%   M         : 信号源数量（迭代次数）
% 输出:
%   P_db      : 空间谱（归一化 dB），非支撑集位置为 -100 dB
%   est_angles: 估计的 DOA 角度值

P = length(grid_theta);
grid_step=grid_theta(2)-grid_theta(1);
% -----------------------------------------------------------------------
% 1. 非圆扩展
% -----------------------------------------------------------------------
Y_aug      = [Y; conj(Y)];               % 扩展观测  (2N x K)
A_dict_aug = [A_dict; conj(A_dict)];     % 扩展字典  (2N x P)

% 归一化字典原子，确保投影公平
atom_norms = vecnorm(A_dict_aug, 2, 1);  % (1 x P)
A_dict_aug = A_dict_aug ./ atom_norms;   % 每列单位化

% -----------------------------------------------------------------------
% 2. 初始化
% -----------------------------------------------------------------------
R           = Y_aug;           % 残差矩阵
Support_Set = zeros(1, M);     % 支撑集（预分配）
P_spectrum  = zeros(P, 1);     % 空间谱强度

% -----------------------------------------------------------------------
% 3. 贪心迭代
% -----------------------------------------------------------------------
for m = 1:M
    % a. 计算所有原子与残差的投影能量（多快拍求和）
    correlation = sum(abs(A_dict_aug' * R).^2, 2);   % (P x 1)

    % b. 屏蔽已选原子，避免重复选取
    correlation(Support_Set(1:m-1)) = 0;

    % c. 选取相关性最强的原子
    [~, best_idx] = max(correlation);
    Support_Set(m) = best_idx;

    % d. 在当前支撑集上做最小二乘（用 \ 代替 pinv，更稳定高效）
    A_sel = A_dict_aug(:, Support_Set(1:m));
    S_hat = A_sel \ Y_aug;                           % (m x K)

    % e. 更新残差
    R = Y_aug - A_sel * S_hat;

    % f. 记录每个已选信源的功率（取各快拍的 Frobenius 范数）
    P_spectrum(Support_Set(1:m)) = vecnorm(S_hat, 2, 2);  % (m x 1)
end
% -----------------------------------------------------------------------
% 4. 【新增】离网格精化：抛物线插值
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% 4. 【离网格精化：仅用抛物线插值，稳定可靠】
% -----------------------------------------------------------------------
N2 = size(Y_aug, 1) / 2;
x_n = (0:N2-1)';
est_angles_refined = zeros(1, M);

for m = 1:M
    idx = Support_Set(m);
    theta_coarse = grid_theta(idx);
    
    % 计算相邻三点的投影能量
    if idx > 1 && idx < P
        c_prev = sum(abs(A_dict_aug(:, idx-1)' * Y_aug).^2, 'all');
        c_curr = sum(abs(A_dict_aug(:, idx  )' * Y_aug).^2, 'all');
        c_next = sum(abs(A_dict_aug(:, idx+1)' * Y_aug).^2, 'all');
        
        denom = c_prev - 2*c_curr + c_next;
        if abs(denom) > 1e-10
            delta = 0.5 * (c_prev - c_next) / denom;
            delta = max(min(delta, 1), -1); % 限制在[-1,1]格内
        else
            delta = 0;
        end
        est_angles_refined(m) = theta_coarse + delta * grid_step;
    else
        est_angles_refined(m) = theta_coarse;
    end
end
% -----------------------------------------------------------------------
% 5. 结果输出
% -----------------------------------------------------------------------
est_angles = grid_theta(Support_Set);

% 空间谱输出（仍基于网格，用于可视化）
P_db = -100 * ones(P, 1);
peak_vals = P_spectrum(Support_Set);
if max(peak_vals) > 0
    P_db(Support_Set) = 20 * log10(peak_vals / max(peak_vals));
end

end