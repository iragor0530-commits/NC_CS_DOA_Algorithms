function crb_rmse = compute_NC_CRB(theta_deg, N, K, sigma2)
%随机模型的CRB计算，适用于NC_MUSIC和NC_ESPRIT 算法比对。
%theta_deg 1×M，M为信源个数
%N：阵元数
%K：快拍数
%sigma2：噪声功率

    M         = length(theta_deg);   % 信源个数
    theta_rad = deg2rad(theta_deg(:)');
    m         = (0:N-1)';

    %% 流形矩阵和导数矩阵
    A = zeros(N, M);
    D = zeros(N, M);
    for k = 1:M
        a_k    = exp(1j * pi * m * sin(theta_rad(k)));
        A(:,k) = a_k;
        D(:,k) = 1j * pi * m * cos(theta_rad(k)) * (pi/180) .* a_k;
    end

    %% NC 增广（BPSK: ψ=0）
    B  = [A;      conj(A)];    % 2N×M
    Db = [D;      conj(D)];    % 2N×M

    %% 增广协方差矩阵
    R_aug = B * B' + sigma2 * eye(2*N);
    R_inv = inv(R_aug);

    %% Fisher 信息矩阵
    FIM = zeros(M, M);
    for i = 1:M
        dRi = Db(:,i)*B(:,i)' + B(:,i)*Db(:,i)';
        for j = 1:M
            dRj = Db(:,j)*B(:,j)' + B(:,j)*Db(:,j)';
            FIM(i,j) = 2 * K * real(trace(R_inv * dRi * R_inv * dRj));
        end
    end

    %% 输出每个信源的 RMSE 下界（度）
    crb_rmse = sqrt(diag(inv(FIM)))';    % 1×M
end