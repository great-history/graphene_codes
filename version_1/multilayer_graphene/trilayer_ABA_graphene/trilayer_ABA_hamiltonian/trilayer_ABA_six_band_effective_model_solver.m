function [eig_enes_K, eig_enes_Kp] = trilayer_ABA_six_band_effective_model_solver(hopping_params, Delta1, akx_mesh, aky_mesh)
    gamma0 = hopping_params(1);
    gamma1 = hopping_params(2);
    gamma2 = hopping_params(3);
    gamma3 = hopping_params(4);  
    gamma4 = hopping_params(5);
    gamma5 = hopping_params(6);
    delta = hopping_params(7);
    Delta2 = hopping_params(8);
    
    dim_x = size(akx_mesh, 2);
    dim_y = size(akx_mesh, 1);
    
    eig_enes_K = zeros(dim_y, dim_x, 6);
    eig_enes_Kp = zeros(dim_y, dim_x, 6);
    
    for i = 1:dim_x
        for j = 1:dim_y
            akx = akx_mesh(j, i);
            aky = aky_mesh(j, i);
            
            [HK_m_ham, HKp_m_ham, HK_b_ham, HKp_b_ham] = construct_trilayer_ABA_six_bands_continuum_model(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky);
            if Delta1 == 0
                [~, D] = eig(HK_m_ham);
                eig_enes_K(j, i, 1:2) = diag(D);

                [~, D] = eig(HK_b_ham);
                eig_enes_K(j, i, 3:6) = diag(D);
                
                [~, D] = eig(HKp_m_ham);
                eig_enes_Kp(j, i, 1:2) = diag(D);

                [~, D] = eig(HKp_b_ham);
                eig_enes_Kp(j, i, 3:6) = diag(D);
                
            else
                Delta1_mat = zeros(2, 4);
                Delta1_mat(1, 1) = Delta1;
                Delta1_mat(2, 2) = Delta1;
                
                HK_ham = [HK_m_ham, Delta1_mat; Delta1_mat', HK_b_ham];
                HKp_ham = [HKp_m_ham, Delta1_mat; Delta1_mat', HKp_b_ham];
                
                [~, D] = eig(HK_ham);
                eig_enes_K(j, i) = diag(D);

                [~, D] = eig(HKp_ham);
                eig_enes_Kp(j, i) = diag(D);
            end
            
        end
    end
    
end