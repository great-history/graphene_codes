function LL_K_b0_eigvec_alpha_component_array = get_LL_b_components_each_alpha(T_mat, LL_K_b0_eigvec_phi_component_array, LL_index_cutoff)
    LL_K_b0_eigvec_alpha_component_array = zeros(LL_index_cutoff + 1, 6);
    eigvec_phi_component = zeros(1,6);
    for n = 0:LL_index_cutoff
        eigvec_phi_component(3:6) = LL_K_b0_eigvec_phi_component_array(n + 1, :);
        for alpha = 1:6
            LL_K_b0_eigvec_alpha_component_array(n + 1, alpha) = dot(T_mat(alpha, :), eigvec_phi_component);
        end
    end
end