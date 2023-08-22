function [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_without_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, ak_start, ak_end, err)
    % 生成在谷附近的 k-points
    akxs = linspace(ak_start(1), ak_end(1), Nx);
    akys = linspace(ak_start(2), ak_end(2), Ny);
    
    eig_enes_m_K = zeros(Nx, Ny, 2);
    eig_enes_m_Kp = zeros(Nx, Ny, 2);
    eig_enes_b_K = zeros(Nx, Ny, 4);
    eig_enes_b_Kp = zeros(Nx, Ny, 4);
    
    for i = 1:Nx
        akx = akxs(i);
        aky = akys(1); % 将第一行的第一个点作为参考
        [H_m_K, H_b_K, H_m_Kp, H_b_Kp] = construct_trilayer_effective_Ham_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, akx, aky);

        [eigvec_H_m_K_last,D] = eig(H_m_K);
        eigval_H_m_K_diag_last = diag(D);
        eig_enes_m_K(i,1,:)=diag(D);
        
        [eigvec_H_b_K_last,D] = eig(H_b_K);
        eigval_H_b_K_diag_last = diag(D);
        eig_enes_b_K(i,1,:)=diag(D);

        [eigvec_H_m_Kp_last,D] = eig(H_m_Kp);
        eigval_H_m_Kp_diag_last = diag(D);
        eig_enes_m_Kp(i,1,:)=diag(D);
        
        [eigvec_H_b_Kp_last,D] = eig(H_b_Kp);
        eigval_H_b_Kp_diag_last = diag(D);
        eig_enes_b_Kp(i,1,:)=diag(D);
        
        for j = 2:Ny
            aky = akys(j);
            [H_m_K, H_b_K, H_m_Kp, H_b_Kp] = construct_trilayer_effective_Ham_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, akx, aky);
            
            [eigvec_H_m_K_now,D] = eig(H_m_K);
            eigval_H_m_K_diag_now = diag(D);

            [eigvec_H_b_K_now,D] = eig(H_b_K);
            eigval_H_b_K_diag_now = diag(D);

            [eigvec_H_m_Kp_now,D] = eig(H_m_Kp);
            eigval_H_m_Kp_diag_now = diag(D);

            [eigvec_H_b_Kp_now,D] = eig(H_b_Kp);
            eigval_H_b_Kp_diag_now = diag(D);
        end
    end
end