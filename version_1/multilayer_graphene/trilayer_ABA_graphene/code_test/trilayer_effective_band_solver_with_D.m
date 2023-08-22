function [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, ak_start, ak_end)
    [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D_no_reorder(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, ak_start, ak_end);
end

function [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D_no_reorder(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, ak_start, ak_end)
% 生成在谷附近的 k-points
    dims = 6;
    akxs = linspace(ak_start(1), ak_end(1), Nx);
    akys = linspace(ak_start(2), ak_end(2), Ny);
    
    eig_enes_K = zeros(Nx, Ny, dims);
    eig_enes_Kp = zeros(Nx, Ny, dims);
    
     for i = 1:Nx
        akx = akxs(i);
        aky = akys(1);
        [H_t_K, H_t_Kp] = construct_trilayer_effective_Ham_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, akx, aky);

        [~, D] = eig(H_t_K);
        eigval_H_t_K_diag_last = diag(D);

        [~, D] = eig(H_t_Kp);
        eigval_H_t_Kp_diag_last = diag(D);

        eig_enes_K(i,1,:)=eigval_H_t_K_diag_last;
        eig_enes_Kp(i,1,:)=eigval_H_t_Kp_diag_last;

        for j = 2:Ny
            aky = akys(j);
            [H_t_K, H_t_Kp] = construct_trilayer_effective_Ham_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, akx, aky);

            [~, D] = eig(H_t_K);
            eigval_H_t_K_diag_now = diag(D);

            [~, D] = eig(H_t_Kp);
            eigval_H_t_Kp_diag_now = diag(D);

            eig_enes_K(i, j, :)=eigval_H_t_K_diag_now;
            eig_enes_Kp(i, j, :)=eigval_H_t_Kp_diag_now;
        end
    end
end