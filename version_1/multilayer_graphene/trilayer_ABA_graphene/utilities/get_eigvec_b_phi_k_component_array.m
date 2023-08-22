function eigvec_b_phi_k_component_array = get_eigvec_b_phi_k_component_array(eigvec_HK_b_now, LL_index_cutoff, num_vec, valley)
    % 得到valley处的bilayer-like中各个phi_k下的LL分量的分布
    
    eigvec_b_phi_k_component_array = zeros(num_vec, LL_index_cutoff + 1, 4);
    for vec_index = 1:num_vec
        eigvec_b_phi_k_component_temp = get_LL_b_components_each_phi_k(eigvec_HK_b_now(:, vec_index), LL_index_cutoff, valley);
        %         if abs(eigvec_b_phi_k_component_temp(1, 4)) > 0.8
        %             disp("this is bilayer LL0+")
        %             fprintf("%d,%d\n",B_index, vec_index)
        %             LLK_b0_index_list(B_index) = vec_index;
        %             count_LLK_b0 = count_LLK_b0 + 1;
        %         elseif abs(eigvec_b_phi_k_component_temp(2, 4)) > 0.8
        %             disp("this is bilayer LL1+")
        %             fprintf("%d,%d\n",B_index, vec_index)
        %             LLK_b1_index_list(B_index) = vec_index;
        %             count_LLK_b1 = count_LLK_b1 + 1;
        %         end
        eigvec_b_phi_k_component_array(vec_index, :, :) = eigvec_b_phi_k_component_temp;
    end
end