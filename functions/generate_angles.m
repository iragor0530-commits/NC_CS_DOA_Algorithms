function theta = generate_angles(M, theta_min, theta_max, min_sep)
    % 生成 M 个在 [theta_min, theta_max] 范围内、
    % 且相邻角度间隔不小于 min_sep 的随机角度
    max_try = 1000;
    for i = 1:max_try
        candidates = theta_min + (theta_max - theta_min) * rand(1, M);
        candidates = sort(candidates);
        if all(diff(candidates) >= min_sep)
            theta = candidates;
            return;
        end
    end
    % 若随机失败，均匀分配
    theta = linspace(theta_min + min_sep, theta_max - min_sep, M);
end