function eigvec_m_phi_k_component_array = get_eigvec_m_phi_k_component_array(eigvec_HK_m_now, LL_index_cutoff, num_vec, valley)
    % 得到valley处的bilayer-like中各个phi_k下的LL分量的分布
    
    eigvec_m_phi_k_component_array = zeros(num_vec, LL_index_cutoff + 1, 2);
    for vec_index = 1:num_vec
        eigvec_m_phi_k_component_temp = get_LL_m_components_each_phi_k(eigvec_HK_m_now(:, vec_index), LL_index_cutoff, valley);
        eigvec_m_phi_k_component_array(vec_index, :, :) = eigvec_m_phi_k_component_temp;
    end
end