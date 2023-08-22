function [H_tmbg_K, H_tmbg_Kp] = construct_tmbg_continuum_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, Delta_asymm, w_aa, w_ab, ...
                                                                akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn)
   % dims = 3 * 2 * num_q_couple : 3是layer, 2是sublattice, N是耦合在一起的波矢数目
   dims_no_layer = 2 * num_q_couple; % 不考虑layer后的维数
   dims = 3 * 2 * num_q_couple;
   
   % transfer momentum
   q_norm = 8 * pi * abs(sin(theta / 2)) / 3; % 乘以晶格常数a，保证波矢无量纲
   q_b = q_norm * [0; -1];
   % q_tr = q_norm * [sqrt(3) / 2; 1 / 2];
   % q_tl = q_norm * [-sqrt(3) / 2; 1 / 2];
   
   % moire lattice vectors
   % bm_r = q_norm * sqrt(3) * [1 / 2; sqrt(3) / 2];
   % bm_l = q_norm * sqrt(3) * [-1 / 2; sqrt(3) / 2];
   % bm_sup_vecs = 0.246 * bm_sup_vecs; % optional, 保证波矢无量纲
   
   %% construct the 1st & 2nd layer ham @ valley K
   H_K_layer1 = zeros(dims_no_layer);
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢
       akx_temp = akx + bm_sup_vecs(1, i);
       aky_temp = aky + bm_sup_vecs(2, i);
       
       H_m_K = construct_HmK_with_rotate(gamma0, akx_temp, aky_temp, -theta / 2); % 第一层顺时针转theta/2
       
       H_m_K = H_m_K + Delta_asymm * [1/2, 0; 0, 1/2];
       H_K_layer1((2 * i - 1):2*i, (2 * i - 1):2*i) = H_m_K;
   end
   
   %% construct the 3rd & 4th layer ham @ valley K
   H_K_layer23 = zeros(dims_no_layer * 2);
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢
       akx_temp = akx + q_b(1) + bm_sup_vecs(1, i);
       aky_temp = aky + q_b(2) + bm_sup_vecs(2, i);
       
       H_b_K = construct_HbK_with_rotate(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx_temp, aky_temp, theta / 2, true); % 第二三层逆时针转theta/2
       H_b_K = H_b_K + Delta_asymm * [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, -1/2, 0; 0, 0, 0 -1/2];
       H_K_layer23((4 * i - 3):4*i, (4 * i - 3):4*i) = H_b_K;
   end
   
   %% construct the 3rd & 4th layer ham @ valley Kp
   H_Kp_layer23 = zeros(dims_no_layer * 2);
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢
       akx_temp = akx + bm_sup_vecs(1, i);
       aky_temp = aky + bm_sup_vecs(2, i);
       
       H_b_Kp = construct_HbKp_with_rotate(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx_temp, aky_temp, theta / 2, true); % 第二三层逆时针转theta/2
       H_b_Kp = H_b_Kp + Delta_asymm * [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, -1/2, 0; 0, 0, 0 -1/2];
       H_Kp_layer23((4 * i - 3):4*i, (4 * i - 3):4*i) = H_b_Kp;
   end
   
   %% construct the 1st & 2nd layer ham @ valley Kp
   H_Kp_layer1 = zeros(dims_no_layer);
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢
       akx_temp = akx + q_b(1) + bm_sup_vecs(1, i);
       aky_temp = aky + q_b(2) + bm_sup_vecs(2, i);
       
       H_m_Kp = construct_HmKp_with_rotate(gamma0, akx_temp, aky_temp, -theta / 2); % 第一层顺时针转theta/2       
       H_m_Kp = H_m_Kp + Delta_asymm * [1/2, 0; 0, 1/2];

       H_Kp_layer1((2 * i - 1):2*i, (2 * i - 1):2*i) = H_m_Kp;
   end
   
   %% construct the interlayer ham
   phi = 2 * pi / 3;
   %% valley K
   % from AA stacking
   % T matrix at valley K
   T_b_K = [w_aa, w_ab; w_ab, w_aa];
   T_tr_K = [w_aa, w_ab * exp(1j*phi); w_ab * exp(-1j*phi), w_aa];
   T_tl_K = [w_aa, w_ab * exp(-1j*phi); w_ab * exp(1j*phi), w_aa];
   
   H_K_layer2_to_layer1 = zeros(dims_no_layer, dims_no_layer * 2);
   for i = 1:num_q_couple % layer 1
       % 第i个参与耦合的波矢
       m_temp = bm_sup_mn(1, i); % bm_r前的系数
       n_temp = bm_sup_mn(2, i); % bm_l前的系数
       
       % from layer2 to layer1 @ K with T_b
       H_K_layer2_to_layer1((2 * i - 1):(2*i), (4 * i - 3):(4*i - 2)) = T_b_K;
       
       % from layer2 to layer1 @ K with T_tr
       count = 0;
       for j = 1:num_q_couple % layer 2
            if count == 2
                break
            end
            
            if (bm_sup_mn(2, j) == n_temp)
               if (bm_sup_mn(1, j) == m_temp + 1)
                    H_K_layer2_to_layer1((2 * i - 1):2*i, (4 * j - 3):(4*j - 2)) = T_tr_K;
                    count = count + 1;
               end
            end
            
            if (bm_sup_mn(1, j) == m_temp)
                if (bm_sup_mn(2, j) == n_temp + 1)
                    H_K_layer2_to_layer1((2 * i - 1):2*i, (4 * j - 3):(4*j - 2)) = T_tl_K;
                    count = count + 1;
               end
            end
       end
   end
   
   %% valley Kp
   % T matrix at valley Kp   
   T_b_Kp = [w_aa, w_ab; w_ab, w_aa];
   T_tr_Kp = [w_aa, w_ab * exp(-1j*phi); w_ab * exp(1j*phi), w_aa];
   T_tl_Kp = [w_aa, w_ab * exp(1j*phi); w_ab * exp(-1j*phi), w_aa];
   
   H_Kp_layer1_to_layer2 = zeros(dims_no_layer * 2, dims_no_layer);
   for i = 1:num_q_couple % layer 2
       % 第i个参与耦合的波矢
       m_temp = bm_sup_mn(1, i); % bm_r前的系数
       n_temp = bm_sup_mn(2, i); % bm_l前的系数
       
       % from layer1 to layer2 @ K with T_b
       H_Kp_layer1_to_layer2((4 * i - 3):(4*i - 2), (2 * i - 1):2*i) = T_b_Kp;
       
       % from layer2 to layer1 @ K with T_tr
       count = 0;
       for j = 1:num_q_couple % layer 1
            if count == 2
                break
            end
            
            if (bm_sup_mn(2, j) == n_temp)
               if (bm_sup_mn(1, j) == m_temp + 1)
                    H_Kp_layer1_to_layer2((4 * i - 3):(4*i - 2), (2 * j - 1):2*j) = T_tr_Kp;
                    count = count + 1;
               end
            end
            
            if (bm_sup_mn(1, j) == m_temp)
                if (bm_sup_mn(2, j) == n_temp + 1)
                    H_Kp_layer1_to_layer2((4 * i - 3):(4*i - 2), (2 * j - 1):2*j) = T_tl_Kp;
                    count = count + 1;
               end
            end
       end
   end
   
   %% get the whole ham
   H_tmbg_K = [H_K_layer1, H_K_layer2_to_layer1; H_K_layer2_to_layer1', H_K_layer23];
   H_tmbg_Kp = [H_Kp_layer1, H_Kp_layer1_to_layer2'; H_Kp_layer1_to_layer2, H_Kp_layer23];
end