function [eig_enes_K, eig_enes_Kp] = bilayer_effective_four_band_model_solver(gamma0, gamma1, gamma3, gamma4, delta, ext_mat_K, ext_mat_Kp, akx_mesh, aky_mesh)
    
    dim_x = size(akx_mesh, 2);
    dim_y = size(akx_mesh, 1);
    
    eig_enes_K = zeros(dim_y, dim_x, 4);
    eig_enes_Kp = zeros(dim_y, dim_x, 4);
    
    for i = 1:dim_x
        for j = 1:dim_y
            akx = akx_mesh(j, i);
            aky = aky_mesh(j, i);
            
            [HK_ham, HKp_ham] = construct_bilayer_continuum_model(gamma0, gamma1, gamma3, gamma4, delta, ext_mat_K, ext_mat_Kp, akx, aky);
            
            [~, D] = eig(HK_ham);
            eig_enes_K(j, i) = diag(D);
            
            [~, D] = eig(HKp_ham);
            eig_enes_Kp(j, i) = diag(D);
        end
    end
    
end