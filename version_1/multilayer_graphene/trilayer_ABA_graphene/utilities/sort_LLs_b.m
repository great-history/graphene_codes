function [eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eigval_LL_K_b0_list, eigval_LL_Kp_b0_list] = sort_LLs_b(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_b_0, LL_Kp_b_0, ...
                                                                              LL_index_cutoff, B_steps, ene_eps, weight_K_b0, weight_Kp_b0)
    num_vec_K_b_min = eig_info_HK_b_select_cell{B_steps, 3}; % 最大磁场下的矢量数目是最少的
    num_vec_Kp_b_min = eig_info_HKp_b_select_cell{B_steps, 3}; % 最大磁场下的矢量数目是最少的
    eigval_LL_K_b0_list = zeros(1, B_steps);
    eigval_LL_Kp_b0_list = zeros(1, B_steps);
    
    % 至少要先保证在两个B_field下LL_K_b0, LL_Kp_b0, LL_K_b1, LL_Kp_b1是确定的
    for B_index = B_steps:-1:(B_steps - 1)
        %% 双层K LL0（只有在高磁场下才行）
        num_vec_K_b = eig_info_HK_b_select_cell{B_index, 3};
        eigvec_HK_b_now = eig_info_HK_b_select_cell{B_index, 1};
        eigval_HK_b_now = eig_info_HK_b_select_cell{B_index, 2};
        
        eigvec_HK_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HK_b_now, LL_index_cutoff, num_vec_K_b, + 1);
        weight_temp = 0.8;
        LL_K_b0_index_list = find_trilayer_ABA_LLL(eigvec_HK_b_phi_k_component_array, eigval_HK_b_now, LL_K_b_0, "bK_0", ene_eps, weight_temp);
        while isempty(LL_K_b0_index_list) % 没有找到，那么把weight稍微调小点
            weight_temp = weight_temp - 0.025;
            LL_K_b0_index_list = find_trilayer_ABA_LLL(eigvec_HK_b_phi_k_component_array, eigval_HK_b_now, LL_K_b_0, "bK_0", ene_eps, weight_K_b0);
        end
        
        if length(LL_K_b0_index_list) > 1% 没有找到或者找到两个，通过斜率来确定
            LL_K_b_0 = eigval_HK_b_now(LL_K_b0_index_list(1));
            LL_K_b0_index = LL_K_b0_index_list(1);
            for jj = 2:length(LL_K_b0_index_list)
                if eigval_HK_b_now(LL_K_b0_index_list(jj)) < LL_K_b_0
                    LL_K_b_0 = eigval_HK_b_now(LL_K_b0_index_list(jj));
                    LL_K_b0_index = LL_K_b0_index_list(jj);
                end
            end
        else
            LL_K_b0_index = find_closest_eigval(eigval_HK_b_now, LL_K_b_0);
        end
        
        LL_K_b_index_list = zeros(num_vec_K_b, 1);
        for ii = 1:num_vec_K_b
            LL_K_b_index_list(ii) = ii - LL_K_b0_index;
        end
        eig_info_HK_b_select_cell{B_index, 4} = LL_K_b_index_list;
        eigval_LL_K_b0_list(B_index) = eig_info_HK_b_select_cell{B_index, 2}(LL_K_b0_index);
        
         %% 双层Kp LL0（只有在高磁场下才行）
         num_vec_Kp_b = eig_info_HKp_b_select_cell{B_index, 3};
        eigvec_HKp_b_now = eig_info_HKp_b_select_cell{B_index, 1};
        eigval_HKp_b_now = eig_info_HKp_b_select_cell{B_index, 2};
        
        eigvec_HKp_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HKp_b_now, LL_index_cutoff, num_vec_Kp_b, - 1);
        weight_temp = 0.8;
        LL_Kp_b0_index_list = find_trilayer_ABA_LLL(eigvec_HKp_b_phi_k_component_array, eigval_HKp_b_now, LL_Kp_b_0, "bKp_0", ene_eps, weight_temp);
        while isempty(LL_Kp_b0_index_list) % 没有找到，那么把weight稍微调小点
            weight_temp = weight_temp - 0.025;
            LL_Kp_b0_index_list = find_trilayer_ABA_LLL(eigvec_HKp_b_phi_k_component_array, eigval_HKp_b_now, LL_Kp_b_0, "bKp_0", ene_eps, weight_K_b0);
        end
        
        if length(LL_Kp_b0_index_list) > 1% 没有找到或者找到两个，通过斜率来确定
            LL_Kp_b_0 = eigval_HKp_b_now(LL_Kp_b0_index_list(1));
            LL_Kp_b0_index = LL_Kp_b0_index_list(1);
            for jj = 2:length(LL_Kp_b0_index_list)
                if eigval_HKp_b_now(LL_K_b0_index_list(jj)) < LL_Kp_b_0
                    LL_Kp_b_0 = eigval_HKp_b_now(LL_Kp_b0_index_list(jj));
                    LL_Kp_b0_index = LL_Kp_b0_index_list(jj);
                end
            end
        else
            LL_Kp_b0_index = find_closest_eigval(eigval_HKp_b_now, LL_Kp_b_0);
        end
        
        LL_Kp_b_index_list = zeros(num_vec_Kp_b, 1);
        for ii = 1:num_vec_Kp_b
            LL_Kp_b_index_list(ii) = ii - LL_Kp_b0_index;
        end
        eig_info_HKp_b_select_cell{B_index, 4} = LL_Kp_b_index_list;
        eigval_LL_Kp_b0_list(B_index) = eig_info_HKp_b_select_cell{B_index, 2}(LL_Kp_b0_index);
    end
    
    for B_index = (B_steps - 2):-1:1
        %% 双层K LL0（只有在高磁场下才行）
        num_vec_K_b = eig_info_HK_b_select_cell{B_index, 3};
        eigvec_HK_b_now = eig_info_HK_b_select_cell{B_index, 1};
        eigval_HK_b_now = eig_info_HK_b_select_cell{B_index, 2};
        if (num_vec_K_b == num_vec_K_b_min) || (num_vec_K_b == num_vec_K_b_min - 1)
            eigvec_HK_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HK_b_now, LL_index_cutoff, num_vec_K_b, + 1);
            LL_K_b0_index = find_trilayer_ABA_LLL(eigvec_HK_b_phi_k_component_array, eigval_HK_b_now, LL_K_b_0, "bK_0", ene_eps, weight_K_b0);
        else % 通过斜率来确定
            LL_K_b_0 = 2 * eigval_LL_K_b0_list(B_index + 1) - eigval_LL_K_b0_list(B_index + 2);
            LL_K_b0_index = find_closest_eigval(eigval_HK_b_now, LL_K_b_0);
        end
        
        LL_K_b_index_list = zeros(num_vec_K_b, 1);
        for ii = 1:num_vec_K_b
            LL_K_b_index_list(ii) = ii - LL_K_b0_index;
        end
        eig_info_HK_b_select_cell{B_index, 4} = LL_K_b_index_list;
        eigval_LL_K_b0_list(B_index) = eig_info_HK_b_select_cell{B_index, 2}(LL_K_b0_index);
        
    %     if (~flag_slope_K_b0) && (num_vec_K_b == num_vec_K_b_min) % 通过LL权重来判定
    %         eigvec_HK_b_now = eig_info_HK_b_select_cell{B_index, 1};
    %         eigval_HK_b_now = eig_info_HK_b_select_cell{B_index, 2};
    %         
    %         eigvec_HK_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HK_b_now, LL_index_cutoff, num_vec_K_b, + 1);
    %         LL_K_b0_index = find_trilayer_ABA_LLL(eigvec_HK_b_phi_k_component_array, eigval_HK_b_now, LL_K_b_0, "bK_0", ene_eps, weight_K_b0);
    %         while isempty(LL_K_b0_index) % 没有找到，那么把weight稍微调小点
    %             disp("降低weight_K_b0")
    %             weight_K_b0
    %             weight_K_b0 = weight_K_b0 - 0.025;
    %             LL_K_b0_index = find_trilayer_ABA_LLL(eigvec_HK_b_phi_k_component_array, eigval_HK_b_now, LL_K_b_0, "bK_0", ene_eps, weight_K_b0);
    %         end
    %         
    %         % if weight_K_b0 < 0
    %         %     disp("weight_K_b0小于0")
    %         %     break
    %         % end
    %         
    %         if isempty(LL_K_b0_index) || length(LL_K_b0_index) > 1% 没有找到或者找到两个，通过斜率来确定
    %             LL_K_b_0 = 2 * eigval_LL_K_b0_list(B_index + 1) - eigval_LL_K_b0_list(B_index + 2);
    %             flag_slope_K_b0 = true;
    %             LL_K_b0_index = find_closest_eigval(eigval_HK_b_now, LL_K_b_0);
    %         end
    %         
    %     else % 通过斜率来确定
    %         LL_K_b_0 = 2 * eigval_LL_K_b0_list(B_index + 1) - eigval_LL_K_b0_list(B_index + 2);
    %         LL_K_b0_index = find_closest_eigval(eigval_HK_b_now, LL_K_b_0);
    %     end

        %% 双层Kp LL0（只有在高磁场下才行）
        num_vec_Kp_b = eig_info_HKp_b_select_cell{B_index, 3};
        eigvec_HKp_b_now = eig_info_HKp_b_select_cell{B_index, 1};
        eigval_HKp_b_now = eig_info_HKp_b_select_cell{B_index, 2};
        if (num_vec_Kp_b == num_vec_Kp_b_min) || (num_vec_Kp_b == num_vec_Kp_b_min - 1)
            eigvec_HKp_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HKp_b_now, LL_index_cutoff, num_vec_Kp_b, - 1);
            LL_Kp_b0_index = find_trilayer_ABA_LLL(eigvec_HKp_b_phi_k_component_array, eigval_HKp_b_now, LL_Kp_b_0, "bKp_0", ene_eps, weight_Kp_b0);
        else % 通过斜率来确定
            LL_Kp_b_0 = 2 * eigval_LL_Kp_b0_list(B_index + 1) - eigval_LL_Kp_b0_list(B_index + 2);
            LL_Kp_b0_index = find_closest_eigval(eigval_HKp_b_now, LL_Kp_b_0);
        end
        
        LL_Kp_b_index_list = zeros(num_vec_Kp_b, 1);
        for ii = 1:num_vec_Kp_b
            LL_Kp_b_index_list(ii) = ii - LL_Kp_b0_index;
        end
        eig_info_HKp_b_select_cell{B_index, 4} = LL_Kp_b_index_list;
        eigval_LL_Kp_b0_list(B_index) = eig_info_HKp_b_select_cell{B_index, 2}(LL_Kp_b0_index);
        
    %     if ~flag_slope_Kp_b0 % 通过LL权重来判定
    %         eigvec_HKp_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HKp_b_now, LL_index_cutoff, num_vec_Kp_b, - 1);
    %         LL_Kp_b0_index = find_trilayer_ABA_LLL(eigvec_HKp_b_phi_k_component_array, eigval_HKp_b_now, LL_Kp_b_0, "bKp_0", ene_eps, weight_Kp_b0);
    %         while isempty(LL_Kp_b0_index) % 没有找到，那么把weight稍微调小点
    %             weight_Kp_b0 = weight_Kp_b0 - 0.025;
    %             LL_Kp_b0_index = find_trilayer_ABA_LLL(eigvec_HKp_b_phi_k_component_array, eigval_HKp_b_now, LL_Kp_b_0, "bKp_0", ene_eps, weight_Kp_b0);
    %         end
    %         
    %         if isempty(LL_Kp_b0_index) || length(LL_Kp_b0_index) > 1 % 没有找到或者找到两个，通过斜率来确定
    %             LL_Kp_b_0 = 2 * eigval_LL_Kp_b0_list(B_index + 1) - eigval_LL_Kp_b0_list(B_index + 2);
    %             flag_slope_Kp_b0 = true;
    %             LL_Kp_b0_index = find_closest_eigval(eigval_HKp_b_now, LL_Kp_b_0);
    %         end
    %     else
    %         LL_Kp_b_0 = 2 * eigval_LL_Kp_b0_list(B_index + 1) - eigval_LL_Kp_b0_list(B_index + 2);
    %         LL_Kp_b0_index = find_closest_eigval(eigval_HKp_b_now, LL_Kp_b_0);
    %     end

    end
end