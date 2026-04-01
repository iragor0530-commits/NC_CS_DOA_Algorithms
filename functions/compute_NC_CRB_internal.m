function crb_deg = compute_NC_CRB_internal(theta_deg, N, SNR_dB, K)
    
    %确定性模型，适用于NC_L1_SVD,NC_OMP,NC_SBL算法。
    %theta_deg:信源角度（1✖M）
    %N：阵元数
    %K：快排数
    %SNR_dB：信噪比

    M   = length(theta_deg);
    rad = deg2rad(theta_deg(:)');
    x_n = (0:N-1)';

    A = exp(1j * pi * x_n * sin(rad));           % N x M
    D = 1j * pi * x_n .* cos(rad) .* A;          % N x M，逐列导数

    % NC 扩展流形
    A_tilde = [A;       conj(A)];                 % 2N x M
    D_tilde = [D;       conj(D)];                 % 2N x M

    % 扩展投影矩阵
    P_perp = eye(2*N) - A_tilde * ((A_tilde' * A_tilde) \ A_tilde');

    snr_lin = 10^(SNR_dB / 10);
    sigma2  = 1 / snr_lin;

    % Fisher 信息矩阵（正确，保留完整矩阵）
    FIM = (2 * K / sigma2) * real(D_tilde' * P_perp * D_tilde);

    crb_deg = (180/pi) * sqrt(diag(inv(FIM)));
end